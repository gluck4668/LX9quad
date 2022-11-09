\name{LX9quad}
\alias{LX9quad}
\title{LX9quad: the Nine-Quadrantal Diagram}
\usage{

LX9quad_p1(gene_file,meta_file,group1,group2)

LX9quad_p2(xmin,xmax,ymin,ymax)

}
\description{
The Nine Quadrantal Diagram of transcriptome and metabolome}
\examples{

install.packages("devtools")
library(devtools)

install_github("gluck4668/LX9quad")

library(LX9quad)

??LX9quad

#------------------------------------
data(gene_data_example)
data(meta_data_example)

#------------------------------------

#rm(list=ls())
setwd("D:/Desktop/R-example/LX9quad_example")

gene_file="genes_data.xlsx"
meta_file="Meta_data.xlsx"

group1="Didymin"
group2="Model"


# To prepare the orginal graphic

LX9quad_p1(gene_file,meta_file,group1,group2)


# The nine-quadrantal diagram were obtained by adjusting the x and y axis based on the orginal graphic above.

xmin=-4
xmax=2.5
ymin=-8
ymax=8

LX9quad_p2(xmin,xmax,ymin,ymax)

}



#-------check the conflict packages-----------------
conflict_scout()
conflict_prefer("cor", "WGCNA")


