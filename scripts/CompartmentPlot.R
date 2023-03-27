library(pheatmap)
library(RColorBrewer)

# oemat = read.table("chr1oe.txt",sep = "\t")
coemat = read.table(snakemake@input[[1]],sep = "\t")

pdf(file = snakemake@output[[1]],width = 20,height = 20)
pheatmap(coemat,
         cluster_rows = F,
         cluster_cols = F,
         border = F,
         show_rownames = F,
         show_colnames = T,
         labels_col = 1:ncol(coemat),
         angle_col = 45,
         fontsize_col = 6)
 dev.off()
 
