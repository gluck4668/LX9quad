
install.packages("devtools")
library(devtools)

install_github("gluck4668/LX9quad")

library(LX9quad)
??LX9quad
#---------------------------------------
data(gene_data_example)
data(meta_data_example) # 可以包含Hightet.Mean
data(meta_example) # 也可以没有Hightet.Mean
#---------------------------------------

rm(list=ls())

devtools::load_all()

  gene_file="gene_data.xlsx"
# meta_file="meta_data.xlsx"
  meta_file="meta_data.xlsx"

# 第一行代码是绘出原始图（x,y轴未经调整）
LX9quad_p1(gene_file,meta_file)


# 第二代码，是根据上图来修改x，y轴，显示比较好看的九象图
xmin=-4
xmax=2
ymin=-8
ymax=8

# mycolor <- c("#3B5FCD","#00CE00","#CD3333","#778899")
mycolor <- c("#3B5FCD","#1edeaa","#ff6666","#9370db") # 自定义颜色

LX9quad_p2(xmin,xmax,ymin,ymax)




