library("ape")
library("ggplot2")
library("ggtree")
library("openxlsx")
library("phangorn")



dist_matrix <- read.table("dist3.txt", sep = "\t", header = TRUE, row.names = 1, dec = ",") 

d <- as.dist(dist_matrix, diag = TRUE, upper = TRUE)


#hclust(d)
#plot(hclust(d))

t <- upgma(d)
write.tree(t, "tree.nwk")

