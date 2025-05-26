############################################################################
#                             ONTOLOGIA DE GENES

library(BiocManager)
BiocManager::install("GEOquery")

citrus <- GSE33003.top.table
View(citrus)
citrus[ ,1]

library(GEOquery)
library(DT) 

annot <- getGEO("GPL5731") 
annot <- Table(annot) 
View(annot) #Tabla de ID con anotaciones y GO

DT::datatable(annot)

