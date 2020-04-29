
library('Seurat')
library(dplyr)
library(ggplot2)

# tsne plot
colors.use <- scales::hue_pal()(length(levels(dds@ident)))
p1 <- TSNEPlot(dds, do.label = FALSE, do.return=TRUE, colors.use = colors.use)
p1$data$pheno <- dds@meta.data$Pheno
p1 + facet_grid(. ~ pheno)


# violin plot for gene
VlnPlot(dds,feature,cols.use = cols.use ,point.size.use = NA, do.return = TRUE,use.imputed=TRUE)+ 
geom_boxplot(aes(fill=vp$data$ident),width=0.2,outlier.colour=NA)
  
## dot plot
print(DotPlot(object = dds, features =c() , plot.legend = TRUE,dot.scale = 20 ,do.return = TRUE,x.lab.rot = TRUE)+
        geom_point(aes(size = pct.exp),shape=21)+
        scale_colour_gradientn(colours=rev(rainbow(4)))+
        theme(legend.position = "bottom")+
        theme(axis.text.x=element_text(size=14,family = "Arial", face = "bold",angle = 45, vjust = 1, hjust=1),
              axis.text.y=element_text(size=14,family = "Arial", face = "bold"))+
        scale_y_discrete(limits = rev(levels(dds@ident)))
)


# heatmap plot
markers <- FindAllMarkers(dds,genes.use = genes.use ,only.pos = TRUE)
top10 <- markers %>% group_by(cluster) %>% top_n(10,avg_logFC)
DoHeatmap(dds, genes.use = top10$gene, group.by = "ident", slim.col.label = TRUE, 
          remove.key = FALSE,rotate.key=TRUE, group.label.rot = TRUE,cex.row=7)+
  theme(legend.position = "top")
