library("ggtree")
library("ggimage")
library("ape")


setwd("C:/b/ss")

t <- read.tree("tree3_clust.nwk")

imgdir <- "C:/b/ss/tree_lab"

ggtree(t) + xlim(0, 9) + geom_tiplab(geom = "label") +
  geom_tiplab(aes(image=paste0(imgdir, '/', label, '.png')), 
                         geom="image", offset=.7, size = 0.2) 
 