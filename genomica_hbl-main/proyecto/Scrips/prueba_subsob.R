
###################

###base

citrus<-("data_table/citrus.tsv")

#Sobreexpresados

sobre<-citrus[which(citrus$adj.P.Val < .05 & citrus$logFC > 0), ]

sobre

#Subexpresados

sube<-citrus[which(citrus$adj.P.Val < .05 & citrus$logFC < 0), ]
sube

######Ver si podemos sacar los nombre convencionales con biomart

###De acuerdo con DAVID y su GO ver si se sub o sobreexpresan

any("AJ012696"==sobre$GB_ACC)##true se sobreexpresa
any("CX070042"==sobre$GB_ACC)##true
any("AY243478"==sobre$GB_ACC)##true
any("AY36204"==sobre$GB_ACC)##false
any("CB292663"==sobre$GB_ACC)##falso
any("AF255013"==sobre$GB_ACC)##falso
any("CB292529"==sobre$GB_ACC)##true


library(biomaRt)


head(biomaRt::listMarts(host = "https://www.ensembl.org"), 10)# ver los marts/fuentes disponibles

##usamos la especializada en genomas vegetales
head(biomaRt::listDatasets(useMart("plants_mart", host = "https://plants.ensembl.org")), )


##atributos
head(listAttributes(biomaRt::useDataset(
  dataset = "cclementina_eg_gene",         
  mart    = useMart("plants_mart", host = "https://plants.ensembl.org"))), 10)        

###filtros
head(listFilters(biomaRt::useDataset(
  dataset = "cclementina_eg_gene",         
  mart    = useMart("plants_mart", host = "https://plants.ensembl.org"))), 10) 


###Establecer la base con la vamos a trabajar
datos<-useDataset(
  dataset = "cclementina_eg_gene",         
  mart    = useMart("plants_mart", host = "https://plants.ensembl.org"))

###
sobre2<-citrus[which(citrus$adj.P.Val < .05 & citrus$logFC > 0), 5]
sube2<-citrus[which(citrus$adj.P.Val < .05 & citrus$logFC < 0), 5]

conjunto<-c(sobre2,sube2)
conjunto
write.csv(conjunto,"resultados_fr/GB_ACC.csv")

identificadores<-getBM(attributes = c("ensembl_gene_id", "go_id"),       
      values     = conjunto,         
      mart       = datos)
identificadores
