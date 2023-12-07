
LX9quad_p1 <- function(gene_file,meta_file){

#-------installing the necessary R packages---------------------

insted.all <- installed.packages()[,1]

com.pack <- c("psych","pheatmap","ggplot2","openxlsx","ggrepel","dplyr","magrittr","htmltools")
com.no <- com.pack[!com.pack %in% insted.all]
com.fun <- function(i){install.packages(i)}
if(length(com.no)>0)
  sapply(com.no,com.fun,simplify = T)

bio.pack <- c("WGCNA","GO.db")
bio.no <- bio.pack[!bio.pack %in% insted.all]
bio.fun <- function(i){BiocManager::install(i)}
if(length(bio.no)>0)
  sapply(bio.no,bio.fun,simplify = T)

lib.fun <- function(i){library(i,character.only = T)}
sapply(c(com.pack,bio.pack),lib.fun)

rm(insted.all,com.pack,com.no,com.fun,bio.pack,bio.no,bio.fun,lib.fun)

#-------reading the gene data --------------------------
gene_all <- read.xlsx(gene_file)
gene_all$sum <- rowSums(gene_all[,3:8])
gene_all <- dplyr::filter(gene_all,sum>0)
gene_all <- gene_all[,-9]
colnames(gene_all)[1:2]=c("gene_id","log2FC")

gene_group <- gsub("\\d+$","", names(gene_all)) %>% .[duplicated(.)]%>% .[!duplicated(.)]
gene_group1 <- gene_group[1]
gene_group2 <- gene_group[2]
print("-----------------------------------------------------------")
group_g <- paste("The groups in the gene data file are",gene_group1,"and",gene_group2)
print(group_g)

gene_df <- data.frame(distinct(gene_all, gene_id, .keep_all = TRUE))
gene_df$log2FC <- as.numeric(gene_df$log2FC)
gene_df <- na.omit(gene_df)
gene_df <- dplyr::filter(gene_df,log2FC!=Inf)
gene_data <- gene_df[,-2]
rownames(gene_data)=gene_data[,1]
genes <- t(gene_data[,-1])
gene_n <- ncol(genes)
genes_FC <- gene_df[,1:2]

#-------reading the metabolite data -------------------
meta_all <- read.xlsx(meta_file)
colnames(meta_all)[1:2]=c("meta_id","FC")

meta_group <- gsub("\\d+$","", names(meta_all)) %>% .[duplicated(.)]%>% .[!duplicated(.)]
group1 <- meta_group[1]
group2 <- meta_group[2]
group1_n <- grep(group1,names(meta_all),ignore.case = T)
group2_n <- grep(group2,names(meta_all),ignore.case = T)

#meta_all$sum <- rowSums(meta_all[,c(group1_n,group2_n)])
#meta_all <- dplyr::filter(meta_all,sum>0)
#meta_all <- meta_all[,-ncol(meta_all)]


group_m <- paste("The groups in the metabolite data file are",group1,"and",group2)
print(group_m)
print("-----------------------------------------------------------")

#------------------------------
meta_df <- data.frame(distinct(meta_all, meta_id, .keep_all = TRUE))
meta_data <- meta_df[,c(group1_n,group2_n)]
rownames(meta_data)=meta_df[,1]
meta <- t(meta_data)
meta_n <- ncol(meta)

FC_df <- meta_df

FC_n <-grep("FC",names(FC_df),ignore.case = T)

FC_df[grepl("inf",FC_df$FC,ignore.case = T),FC_n] <- 100
FC_df$FC <- as.numeric(FC_df$FC)

#FC_df$Highest.Mean <- ifelse(rowMeans(FC_df[group2_n])/rowMeans(FC_df[group1_n])>1, group_meta[2], group_meta[1])
FC_df$FC <- ifelse(rowMeans(FC_df[group2_n])/rowMeans(FC_df[group1_n])>1, log2(FC_df$FC), log2(1/FC_df$FC))
names(FC_df) [c(1,2)]<- c("meta_id","log2FC")
meta_FC <- FC_df[,c(1,2)]


#--------calculating the correlation coefficient------------

# r_only <- cor(genes,meta,method="spearman")
if(nrow(genes)!=nrow(meta))
  stop("The group numbers between the genes and metabolites data are unequal. Please check it" )

cor_r_p <- WGCNA::corAndPvalue(genes,meta,use = "pairwise.complete.obs",method="spearman")
cor_r <- data.frame(cor_r_p[["cor"]])
cor_p <- data.frame(cor_r_p[["p"]])

r_m <- as.matrix(cor_r)
p_m <- as.matrix(cor_p)

dim(p_m)

x <- matrix(NA,ncol = 4,nrow = gene_n*meta_n)

  for(i in 1:gene_n){
    j <- 1
    while(j<=meta_n)
    {x[(i-1)*meta_n+j,1]  <- row.names(r_m)[i]
    x[(i-1)*meta_n+j,2]  <- colnames(r_m)[j]
    j  <- j+1
    }
  }


  for(i in 1:dim(x)[1]){
    x[i,3]  <- r_m[x[i,1],x[i,2]]
    x[i,4]  <- p_m[x[i,1],x[i,2]]
  }


colnames(x) <- c("gene","meta","r-value","p-value")

colnames(genes_FC) <- c("gene","Log2FC_genes")
genes_FC_m <- as.matrix(genes_FC)

colnames(meta_FC) <- c("meta","Log2FC_metabolites")
meta_FC_m <- as.matrix(meta_FC)

genes_meta_0 <- merge(x,genes_FC_m,by="gene")
genes_meta <- merge(genes_meta_0,meta_FC_m,by="meta")

genes_meta[,3] <- as.numeric(genes_meta [,3])
genes_meta[,4] <- as.numeric(genes_meta[,4])
genes_meta[,5] <- as.numeric(genes_meta[,5])
genes_meta[,6] <- as.numeric(genes_meta[,6])

new_genes_meta_0 <- genes_meta[which(abs(genes_meta$`r-value`)>=0.8),]
new_genes_meta_1 <- new_genes_meta_0[which(genes_meta$`p-value` <0.05),]

new_genes_meta <- new_genes_meta_1[ complete.cases(new_genes_meta_1[,5:6]),]
data <- new_genes_meta_1[ complete.cases(new_genes_meta_1[,5:6]),]

data$part <- case_when(abs(data$Log2FC_genes) >= 1 & abs(data$Log2FC_meta) >= 1 ~ "part1379",
                         abs(data$Log2FC_genes) < 1 & abs(data$Log2FC_meta) > 1 ~ "part28",
                         abs(data$Log2FC_genes) > 1 & abs(data$Log2FC_meta) < 1 ~ "part46",
                         abs(data$Log2FC_genes) < 1 & abs(data$Log2FC_meta) < 1 ~ "part5")


sig_cor <- dplyr::filter(data,part=="part1379")

sig_cor$correlation <- case_when(sig_cor$Log2FC_genes >= 1 & sig_cor$Log2FC_meta >= 1 ~ "positive",
                                  sig_cor$Log2FC_genes <= 1 & sig_cor$Log2FC_meta <= 1 ~ "positive",
                                  sig_cor$Log2FC_genes >= 1 & sig_cor$Log2FC_meta <= 1 ~ "negative",
                                  sig_cor$Log2FC_genes <= 1 & sig_cor$Log2FC_meta >= 1 ~ "negative")

postive_cor <- dplyr::filter(sig_cor,correlation=="positive")
postive_cor <- postive_cor[,-7]

negative_cor <- dplyr::filter(sig_cor,correlation=="negative")
negative_cor <-negative_cor[,-7]

if(!dir.exists("analysis result"))
      dir.create("analysis result")

write.xlsx(data,"analysis result/data_cor.xlsx")
write.xlsx(postive_cor,"analysis result/postive_cor.xlsx")
write.xlsx(negative_cor,"analysis result/negative_cor.xlsx")

p0 <-ggplot(data,aes(Log2FC_genes,Log2FC_metabolites,color=part))

p1 <- p0+geom_point(size = 2)+guides(color="none")

mycolor <- c("#3B5FCD","#00CE00","#CD3333","#778899")
p2 <- p1 + scale_colour_manual(name="",values=mycolor)

titile_text <- paste('Nine Quadrantal Diagram',"(",group2,"VS",group1,")")

p3 <- p2+geom_hline(yintercept = c(-1,1),
                      linewidth = 0.6,
                      color = "gray30",
                      lty = "dashed")+
    geom_vline(xintercept = c(-1,1),
               linewidth = 0.6,
               color = "gray30",
               lty = "dashed")+
    labs(title =titile_text) +
    theme(plot.title = element_text(hjust = 0.5))

p3

  }

#--------------------------------------------------------------------#

LX9quad_p2 <- function(xmin,xmax,ymin,ymax){

meta_group <- read.xlsx(meta_file)
meta_group <- gsub("\\d+$","", names(meta_group)) %>% .[duplicated(.)]%>% .[!duplicated(.)]

group1 <- meta_group[2]
group2 <- meta_group[1]

data <- read.xlsx("analysis result/data_cor.xlsx")

p0 <-ggplot(data,aes(Log2FC_genes,Log2FC_metabolites,color=part))

p1 <- p0+geom_point(size = 2)+guides(color="none")

if(!exists("mycolor"))
mycolor <- c("#3B5FCD","#00CE00","#CD3333","#778899")

p2 <- p1 + scale_colour_manual(name="",values=mycolor)

titile_text <- paste('Nine Quadrantal Diagram',"(",group1,"VS",group2,")")

p3 <- p2+geom_hline(yintercept = c(-1,1),
                        linewidth = 0.6,
                        color = "gray30",
                        lty = "dashed")+
      geom_vline(xintercept = c(-1,1),
                 linewidth = 0.6,
                 color = "gray30",
                 lty = "dashed")+
      labs(title =titile_text) +
      theme(plot.title = element_text(hjust = 0.5))


  xmin <- as.numeric(xmin)
  xmax <- as.numeric(xmax)
  ymin <- as.numeric(ymin)
  ymax <- as.numeric(ymax)

p4<-p3+
    scale_y_continuous(expand=expansion(add = c(0.01, 0.01)),
                       limits = c(ymin, ymax),
                       breaks = c(ymin,ymin/2,0,ymax/2,ymax),
                       label = c(as.character(ymin),as.character(ymin/2),"0",
                                 as.character(ymax/2),as.character(ymax)))+
    scale_x_continuous(expand=expansion(add = c(0.01, 0.01)),
                       limits = c(xmin, xmax),
                       breaks = c(xmin,xmin/2,0,xmax/2,xmax),
                       label = c(as.character(xmin),as.character(xmin/2),"0",
                                 as.character(xmax/2),as.character(xmax)))+
  theme_bw()+theme(panel.grid = element_blank())


  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",face="bold",colour ="black",size =14),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=18))+
    theme(axis.text.y = element_text(face="bold",color="black",size=18))

  cor  <- WGCNA::cor(data$Log2FC_genes,data$Log2FC_meta,use = "complete.obs",method ="spearman")
  cor <- round(cor,4)

  lab = paste("correlation=",cor,sep="")
  lab

  cor_value <- geom_text(x=xmin+1,y=ymax-0.5,label = lab, size=4,color="black")

p5 <- p4+mytheme+xytheme+
    annotate("text", x=-1+(xmin+1)/2, y=1+(ymax-1)/2, label= "①",size=12)+
    annotate("text", x=-1+(xmin+1)/2, y=0.1, label= "④", size=12)+
    annotate("text", x=-1+(xmin+1)/2, y=-1+(ymin+1)/2, label= "⑦",size=12)+

    annotate("text", x=0, y=1+(ymax-1)/2, label= "②",size=12)+
    annotate("text", x=0, y=0.1, label= "⑤",size=14)+
    annotate("text", x=0, y=-1+(ymin+1)/2, label= "⑧",size=12)+

    annotate("text", x=1+(xmax-1)/2, y=1+(ymax-1)/2, label= "③",size=12)+
    annotate("text", x=1+(xmax-1)/2, y=0.1, label= "⑥",size=14)+
    annotate("text", x=1+(xmax-1)/2, y=-1+(ymin+1)/2, label= "⑨",size=12)


ggsave("analysis result/The Nine Quadrantal Diagram.png", p5,width=1200, height =900, dpi=150,units = "px")

print("Please see the results in the folder of <analysis result>")

p5

}



