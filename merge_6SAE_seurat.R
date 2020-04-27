library('Seurat')
library(Matrix)
options(stringsAsFactors=FALSE)

samples <- data.frame(sample=c('DGM.00384.SAE', 'DGM.13427.SAE', 'DGM.13434.SAE', 'DGM.13451.SAE', 'DGM.13460.SAE', 'DGM.13471.SAE'),
                      pheno=c('S', 'NS', 'S', 'NS', 'NS', 'S'))


make.dds <- function (sample.id) {
  
  dge <- readMM(paste0(sample.id,".txt"))
  # colnames(dd$dge) <- paste(sample.id, colnames(dd$dge), sep='_')
  
  dds <- CreateSeuratObject(raw.data = dge, min.cells = 10, min.genes = 200,
                            project = sample.id)
  
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = dds@data), value = TRUE)
  percent.mito <- Matrix::colSums(dds@raw.data[mito.genes, ])/Matrix::colSums(dds@raw.data)
  dds <- AddMetaData(dds, metadata = percent.mito, col.name = "percent.mito")
  
  
  ribo.genes <- grep(pattern = "^RP[LS]", x = rownames(x = dds@data), value = TRUE)
  percent.ribo <- Matrix::colSums(dds@raw.data[ribo.genes, ])/Matrix::colSums(dds@raw.data)
  dds <- AddMetaData(dds, metadata = percent.ribo, col.name = "percent.ribo")
  
  dds <- FilterCells(object = dds,
                     subset.names = c("percent.mito", "nUMI"),
                     high.thresholds = c(0.25, 10000))
  dds <- NormalizeData(dds)
  dds <- ScaleData(dds, vars.to.regress="percent.mito")
  dds <- FindVariableGenes(dds)
  
  return (dds)
}

samples <- samples[order(samples$pheno), ]

dds.list <- lapply(samples$sample, make.dds)
names(dds.list) <- samples$sample

dds.list <- lapply(samples$sample, function (x) {
  dds.list[[x]]@meta.data$Pheno = samples$pheno[match(x, samples$sample)]
  dds.list[[x]]
})
names(dds.list) <- samples$sample

cell.ids <- Reduce(c, lapply(dds.list, function (x) colnames(x@data)))
duplicated.ids <- cell.ids[duplicated(cell.ids)]
dds.list <- lapply(dds.list, function (x) SubsetData(x, cells.use=colnames(x@data)[!colnames(x@data) %in% duplicated.ids]))

n.merge.genes <- 500
common.genes <- Reduce(intersect, lapply(dds.list, function (x) rownames(x@data)))
merge.genes <- Reduce(union, lapply(dds.list, function (x) head(rownames(x@hvg.info)[rownames(x@hvg.info) %in% common.genes], n.merge.genes)))

dds <- Reduce(MergeSeurat, dds.list)
dds <- ScaleData(dds)
dds@scale.data[is.na(x = dds@scale.data)] <- 0
dds@var.genes <- merge.genes


dds <- RunPCA(dds, pcs.print = 1:20)

dims.use <- c(1, 2, 5, 7, 8, 9, 11)

dds <- RunTSNE(dds, dims.use = dims.use)
dds <- FindClusters(dds, dims.use = dims.use, resolution=0.5, save.SNN = TRUE)


aggfun <- function (x) -1*(mean(x)+log(sum(x)))
renumber.factor <- function (x) {
  old.levels <- levels(x)
  new.levels <- as.character(1:length(old.levels))
  newx <- ordered(new.levels[match(as.character(x), old.levels)], levels=new.levels)
  names(newx) <- names(x)
  newx
}

dds <- ReorderIdent(dds, feature = 'KRT5', use.scaled=FALSE,
                    aggregate.fxn = aggfun)
dds@ident <- renumber.factor(dds@ident)

dds <- RunMagic(dds)


