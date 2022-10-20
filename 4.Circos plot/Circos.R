#Import libraries
library(RCircos)
library(tidyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)

#Import reference dataset
data("UCSC.HG38.Human.CytoBandIdeogram")
cyto.info = UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
cyto.info = UCSC.HG38.Human.CytoBandIdeogram
cyto.info$Name = NA
cyto.info$Stain = NA
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)
chr_order = unique(cyto.info$Chromosome)

RCircos.Set.Plot.Area()
ideo = RCircos.Get.Plot.Ideogram()
ideo$BandColor = 'Yellow'
num = which(ideo$Chromosome == 'chrX')
ideo[num, 'BandColor'] = 'chartreuse'

num = which(ideo$Chromosome == 'chrY')
ideo[num, 'BandColor'] = 'Red'


RCircos.Reset.Plot.Ideogram(ideo)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#Import the Logcpm file
mat = read.csv("C:/Users/joash/OneDrive/Desktop/Anna/New.csv")

#Convert ensemblid to gene symbols
ensl=mat[,1]
mat$gene<-mapIds(org.Hs.eg.db,keys=ensl,keytype = "ENSEMBL",column = "SYMBOL")
mat = mat %>% drop_na()
rownames(mat) =  make.names(mat[,6],unique = T)
mat2 = mat[,-1]
mat2 = mat2[,-5]

#Run mart 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

coords = getBM(attributes=c('chromosome_name', 'start_position', 
                            'end_position', 'hgnc_symbol'),
               filters = c('hgnc_symbol'),
               values = list(rownames(mat2)),
               mart = mart)

write.csv(coords, file = 'coords.csv')
coords$chromosome_name = paste0('chr', coords$chromosome_name)
coords$chromosome_name = factor(coords$chromosome_name,levels = chr_order)

num = which(is.na(coords$chromosome_name))
coords = coords[-num, ]

#Show selected values in up
up = which((mat2$Pvalue < 0.01) &
             (mat2$Log2FC > 1))

#Create a dataset of selected values in "up" variables
upmat = mat2[up, ]

#Select gene symbols common in coords and upmat
num = which(coords$hgnc_symbol %in% rownames(upmat))

#Create another dataset that has the common symbols only
coords1 = coords[num, ]


#Plot RCircos plot 
RCircos.Gene.Name.Plot(coords1, name.col=4, track.num = 2, side = "in",
                       is.sorted = F)


#Common genes from coords and our dataset
genes = intersect(rownames(mat2), coords$hgnc_symbol)

mat1 = mat2[genes, ]

#Append selected columns from mat1 
df = cbind.data.frame(rownames(mat1), mat1[, c(1,2,4)])

#Change the name of the column
colnames(df)[1] = 'hgnc_symbol'

#Merge coords and df
data = merge(coords, df, by = 'hgnc_symbol',all.x = T)

#Drop na values
data = data %>% drop_na()

#Rearrange data columns
colnames(data)
data = data[,c(2,3,4,1,5,6,7)]
data = data[, c('chromosome_name', 'start_position',
                'end_position','hgnc_symbol','meanTumor', 'meanControl', 'Log2FC')]
  
#Plot Circos Heatmap
RCircos.Heatmap.Plot(data, data.col = 7, track.num = 6, side = "in",
                     min.value = -2, max.value =2, genomic.columns = 3,
                     is.sorted = F)
RC.param = RCircos.Get.Plot.Parameters()
RC.param['heatmap.color'] = "BlueWhiteRed"
RCircos.Reset.Plot.Parameters(RC.param)

RCircos.Heatmap.Plot(data, data.col = 7, track.num = 10, side = "in",
                     min.value = -2, max.value = 2,
                     is.sorted = F)