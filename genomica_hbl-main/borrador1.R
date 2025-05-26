

#DADA2
library(Rcpp)
library (dada2)

##### Cargar biblioteca y datos
path_sano <- "archivos_fastq"
list.files(path_sano)

##### Extraer los nombres de las muestras
# Forward and reverse fastq tienen el formato: SAMPLENAME_1.fastq y SAMPLENAME_1.fastq
s1F_sano <- sort(list.files(path_sano, pattern="SRR17677269_1.fastq", full.names = TRUE))
s1R_sano <- sort(list.files(path_sano, pattern="SRR17677269_2.fastq", full.names = TRUE))
# Extraer nombres de las muestras, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename (s1F_sano), "_"), `[`, 1)
sample.names

#Para que corra: readRDS("dada2/seqtab.RDS")


## Verifcar calidad PHRED
length(s1F_sano)
length(s1R_sano)

##Forward
pdf("resultados_fr/Q_Forward_sano.pdf", width=13, height = 8)
plotQualityProfile(s1F_sano[1])
dev.off()

##Reverse
pdf("resultados_fr/Q_Reverse_sano.pdf", width=13, height = 8)
plotQualityProfile(s1R_sano[1])
dev.off()


##### Filtrar y Cortar
# Colocar las SRR filtradas en un subdirectorio 
# Cree el forder manualmente para que no se lleve más tiempo 

s1_filtF_sano<- file.path(path_sano, "filtered", paste0(sample.names, "_s1F_filt_sano.fastq.gz"))
s1_filtR_sano<- file.path(path_sano, "filtered", paste0(sample.names, "_s1R_filt_sano.fastq.gz"))
names(s1_filtF_sano) <- sample.names
names(s1_filtR_sano) <- sample.names

##### Filtrar por calidad 
s1_out_sano<- filterAndTrim(s1F_sano, s1_filtF_sano, s1R_sano, s1_filtR_sano, truncLen = c(120,120),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread= TRUE) # Este ultimo FALSE porque está en windows


##### Aprender de los errores 
s1_errF_sano<- learnErrors(s1_filtF_sano, multithread=TRUE)
s1_errR_sano<- learnErrors(s1_filtR_sano, multithread=TRUE)

#Gráficas de tasas de error
png("resultados_fr/Errores_Forward_sano.png")
plotErrors(s1_errF_sano, nominalQ=TRUE)
dev.off()

png("resultados_fr/Errores_Reverse_sano.png")
plotErrors(s1_errR_sano, nominalQ=TRUE)
dev.off()

##inferencia
s1_dadaF_sano <- dada(s1_filtF_sano, err=s1_errF_sano, multithread=TRUE)
s1_dadaR_sano<- dada(s1_filtR_sano, err=s1_errR_sano, multithread=TRUE)


##### Mezclar, hacer de las 2 direcciones
#| echo: true
s1_mergers_sano <- mergePairs(s1_dadaF_sano, s1_filtF_sano, s1_dadaR_sano, s1_filtR_sano, verbose=TRUE)
saveRDS(s1_mergers_sano, file="resultados_fr/s1_mergers_sano.RDS")

s1_mergers_sano
s1_seqtab_sano <- makeSequenceTable (s1_mergers_sano)
rownames(s1_seqtab_sano) <- sample.names##asingnar nombre a las filas sino no corre después
dim (s1_seqtab_sano)
saveRDS(s1_seqtab_sano, file="resultados_fr/seqtab_sano.RDS")


##### Remover quimeras
s1_seqtab_nochim_sano <- removeBimeraDenovo(s1_seqtab_sano, method="consensus", multithread=TRUE, verbose=TRUE)

rownames(s1_seqtab_nochim_sano)###revisando el error de las filas

sum (s1_seqtab_nochim_sano)/sum(s1_seqtab_sano)

#Save RDS --> se refiere a los objetos que ya están en la carpeta resultados

##### Quitar el ruido 
s1_getN_sano <- function(x) sum(getUniques(x))
s1_track_sano <- cbind(s1_out_sano,
                  s1_getN_sano(s1_dadaF_sano),
                  s1_getN_sano(s1_dadaR_sano),
                  s1_getN_sano(s1_mergers_sano),
                  rowSums(s1_seqtab_nochim_sano))
colnames(s1_track_sano) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(s1_track_sano) <- sample.names
head(s1_track_sano)
#denoised, quitar el ruido de forwards y reverse
#filtered, filtrado 
#merged son las que están unidas, que se toman como una abundancia?

#### Asignación taxonómica


# Género 
s1_taxa_sano<- assignTaxonomy(s1_seqtab_nochim_sano, "resultados_fr/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
saveRDS (s1_taxa_sano, file="resultados_fr/taxa_sano.RDS")
readRDS (file="resultados_fr/taxa_sano.RDS")

s1_taxa_sano<- assignTaxonomy(s1_taxa_sano, "resultados_fr/silva_nr99_v138.2_toSpecies_trainset.fa.gz")
saveRDS (s1_taxa_sano, file="resultados_fr/taxa_sano.RDS")
readRDS ("resultados_fr/taxa_sano.RDS")

# Especies
s1_taxa_sano<- addSpecies(s1_taxa_sano, "resultados_fr/silva_v138.2_assignSpecies.fa.gz")
saveRDS (s1_taxa_sano, file="resultados_fr/taxa_sano.RDS")
readRDS ("resultados_fr/taxa_sano.RDS")


##### Objeto más simple es una tabla
s1_taxa_sano_print <- s1_taxa_sano
rownames(s1_taxa_sano_print) <- NULL
head(s1_taxa_sano_print)


#### Data frame con metadatos 
library(phyloseq)
library(Biostrings)
library(ggplot2)

#Hay que editarlo, pero aún no se las características con la que está organizado el objeto
theme_set(theme_bw())
s1_samples_out_sano <- rownames(s1_seqtab_nochim_sano)
class(s1_samples_out_sano)
sujeto_sano<- sapply(strsplit(s1_samples_out_sano, "D"), `[`, 1)#tenía subect como nombre pero eso confunde a R

gender <- substr(sujeto_sano,1,1)
sujeto_sano <- substr(sujeto_sano,2,999)
day <- as.integer(sapply(strsplit(s1_samples_out_sano, "D"), `[`, 2))
s1_df_sano <- data.frame(Subject=sujeto_sano, Gender=gender, Day=day)
s1_df_sano$When <- "Early"
s1_df_sano$When[s1_df_sano$Day>100] <- "Late"
rownames(s1_df_sano) <- s1_samples_out_sano


#### Phyloseq
ps1_sano <- phyloseq(otu_table(s1_seqtab_nochim_sano, taxa_are_rows=FALSE), 
                sample_data(s1_df_sano), 
                tax_table(s1_taxa_sano))
ps1_sano<- prune_samples(sample_names(ps1_sano) != "Mock", ps1_sano) # Remove mock sample
dna1_sano <- Biostrings::DNAStringSet(taxa_names(ps1_sano))

names(dna1_sano) <- taxa_names(ps1_sano)
ps1_sano <- merge_phyloseq(ps1_sano, dna1_sano)

taxa_names(ps1_sano) <- paste0("ASV", seq(ntaxa(ps1_sano)))
ps1_sano

save(ps1_sano,file="resultados_fr/ps_sano.RDS")


#####ANÁLISIS

otu<-otu_table(ps1_sano)###OTU table

tax_table(ps1_sano)###TAX table

sample_data(ps1_sano)

plot_richness(otu,measures = c("Chao1", "Simpson", "Shannon"))##tenemos una sola muestra







































































































































