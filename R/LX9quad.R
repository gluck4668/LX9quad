
LX9quad_p1 <- function(gene_file,meta_file,group1,group2){


   R_packs_install <- function(){

     R_packs <- c("devtools", "psych","pheatmap","ggplot2","openxlsx","ggrepel","dplyr","magrittr")

    list_installed <- installed.packages()

    new_R_packs <- subset(R_packs, !(R_packs %in% list_installed[, "Package"]))

    if(length(new_R_packs)!=0){install.packages(new_R_packs,force=TRUE,quietly = TRUE)
      print(c(new_R_packs, " packages added..."))}

    if((length(new_R_packs)<1)){print("No new dependency packages added...")}

  }

  R_packs_install()


  #---------------------------------


  Bio_packages <- function(){
    Bio_pkgs <- c("WGCNA","GO.db")
    list_installed <- installed.packages()
    new_pkgs <- subset(Bio_pkgs, !(Bio_pkgs %in% list_installed[, "Package"]))
    if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      library(BiocManager)
      BiocManager::install(new_pkgs,force=TRUE,quietly = TRUE)
      print(c(new_pkgs, " packages added..."))
    }

    if((length(new_pkgs)<1)){
      print("No new BiocManager packages added...")
    }
  }


  Bio_packages()


  #---------------------------------

  packages <- c("devtools", "psych","pheatmap","ggplot2","openxlsx","ggrepel",
                "dplyr","magrittr","WGCNA","GO.db")

  for(i in packages){
    library(i, character.only = T)
  }

  rm(i)

 #---------------------------------


  group1 <- trimws(group1)
  group2 <- trimws(group2)

  group <- data.frame(group1,group2)

  if(dir.exists("temporary files")==FALSE)
    dir.create("temporary files")

  write.xlsx(group,"temporary files/group.xlsx")

  gene_all <- read.xlsx(gene_file)
  colnames(gene_all)[1]="gene_id"

  gene_df <- data.frame(distinct(gene_all, gene_id, .keep_all = TRUE))

  gene_data <- gene_df[,-2]
  rownames(gene_data)=gene_data[,1]
  genes <- t(gene_data[,-1])

  gene_n <- ncol(genes)

  genes_FC <- gene_df[,1:2]

  meta_all <- read.xlsx(meta_file)
  colnames(meta_all)[1]="meta_id"

  meta_df <- data.frame(distinct(meta_all, meta_id, .keep_all = TRUE))

  meta_data <- meta_df[,-2:-3]
  rownames(meta_data)=meta_data[,1]
  meta <- t(meta_data[,-1])

  meta_n <- ncol(meta)

  meta_FC_data <- meta_df[,1:3]
  colnames(meta_FC_data) <- c("meta_id","FC","group")
  meta_FC_data$log2FC <- NA
  meta_FC_data$FC <- as.numeric(meta_FC_data[,2])

  #--------------------------

  for(i in 1:meta_n){
    if(substring(tolower(meta_FC_data[i,2]),1,3)=="inf" & tolower(meta_FC_data[i,3])==tolower(group1) )
      meta_FC_data[i,2] <- 100 }
  rm(i)

  for(i in 1:meta_n){
    if(substring(tolower(meta_FC_data[i,2]),1,3)=="inf" & tolower(meta_FC_data[i,3])==tolower(group2))
      meta_FC_data[i,2] <- 100}
  rm(i)

  #--------------------------

  for(i in 1:meta_n){
    if(tolower(trimws(meta_FC_data[i,3]))==tolower(group2))
      meta_FC_data[i,4] <- log2(1/meta_FC_data[i,2])
    else
      meta_FC_data[i,4] <- log2(meta_FC_data[i,2])
  }

  meta_FC <- meta_FC_data[,-2]
  meta_FC <- meta_FC[,-2]

  r_only <- cor(genes,meta,method="spearman")

  cor_r_p <- corAndPvalue(genes,meta,use = "pairwise.complete.obs",method="spearman")

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

  write.xlsx(data,"temporary files/data.xlsx")

  p0 <-ggplot(data,aes(Log2FC_genes,Log2FC_metabolites,color=part))

  p1 <- p0+geom_point(size = 2)+guides(color="none")

  mycolor <- c("#3B5FCD","#00CE00","#CD3333","#778899")
  p2 <- p1 + scale_colour_manual(name="",values=mycolor)

  titile_text <- paste('Nine Quadrantal Diagram',"(",group1,"VS",group2,")")

  p3 <- p2+geom_hline(yintercept = c(-1,1),
                      size = 0.6,
                      color = "gray30",
                      lty = "dashed")+
    geom_vline(xintercept = c(-1,1),
               size = 0.6,
               color = "gray30",
               lty = "dashed")+
    labs(title =titile_text) +
    theme(plot.title = element_text(hjust = 0.5))


#--------------------------
gene_group <- colnames(gene_all)[c(3:8)]
gene_group <- gsub("\\d+$","",gene_group)

if(tolower(group1) %in% tolower(gene_group)==FALSE)
{g01 <- paste("'",group1,"'","is not found in the gene data file. Please check it")
print(g01)}

if(tolower(group2) %in% tolower(gene_group)==FALSE)
{g02 <- paste("'",group2,"'","is not found in the gene data file. Please check it")
print(g02)}



meta_group <- colnames(meta_all)[c(4:9)]
meta_group <- gsub("\\d+$","",meta_group)

if(tolower(group1) %in% tolower(meta_group)==FALSE)
{m01 <- paste("'",group1,"'","is not found in the metabolite data file. Please check it")
print(m01)}

if(tolower(group2) %in% tolower(meta_group)==FALSE)
{m02 <- paste("'",group2,"'","is not found in the metabolite data file. Please check it")
print(m02)}

#--------------------------


  p3


  }

#--------------------------------------------------------------------#

LX9quad_p2 <- function(xmin,xmax,ymin,ymax){

    data <- read.xlsx("temporary files/data.xlsx")
    group <- read.xlsx("temporary files/group.xlsx",colNames=TRUE)

    #file.remove(c("data.xlsx","group.xlsx"))

    p0 <-ggplot(data,aes(Log2FC_genes,Log2FC_metabolites,color=part))

    p1 <- p0+geom_point(size = 2)+guides(color="none")

    mycolor <- c("#3B5FCD","#00CE00","#CD3333","#778899")
    p2 <- p1 + scale_colour_manual(name="",values=mycolor)


    titile_text <- paste('Nine Quadrantal Diagram',"(",group$group1,"VS",group$group2,")")

    p3 <- p2+geom_hline(yintercept = c(-1,1),
                        size = 0.6,
                        color = "gray30",
                        lty = "dashed")+
      geom_vline(xintercept = c(-1,1),
                 size = 0.6,
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
          panel.border = element_rect(size = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(size = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=18))+
    theme(axis.text.y = element_text(face="bold",color="black",size=18))

  cor  <- round(cor(data$Log2FC_genes,data$Log2FC_meta,
                    use = "complete.obs",method ="spearman"),4)
  cor

  lab = paste("correlation=",cor,sep="")
  lab

  cor_value <- geom_text(x=xmin+1,y=ymax-0.5,label = lab, size=4,color="black")

  p5 <- p4+mytheme+xytheme+
        annotate("text", x=-1+(xmin+1)/2, y=ymax/2, label= "①",size=14)+
        annotate("text", x=-1+(xmin+1)/2, y=0.1, label= "④", size=14)+
        annotate("text", x=-1+(xmin+1)/2, y=ymin/2, label= "⑦",size=14)+

        annotate("text", x=0, y=ymax/2, label= "②",size=14)+
        annotate("text", x=0, y=0.1, label= "⑤",size=14)+
        annotate("text", x=0, y=ymin/2, label= "⑧",size=14)+

        annotate("text", x=1+(xmax-1)/2, y=ymax/2, label= "③",size=14)+
        annotate("text", x=1+(xmax-1)/2, y=0.1, label= "⑥",size=14)+
        annotate("text", x=1+(xmax-1)/2, y=ymin/2, label= "⑨",size=14)


  p5

}



