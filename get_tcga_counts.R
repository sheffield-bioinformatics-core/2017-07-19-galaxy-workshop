source("http://www.bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")

library(TCGAbiolinks)

listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
                 "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
                 "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
                 "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
                 "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07")

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(BRCAMatrix),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(BRCAMatrix), 
                                   typesample = c("TP"))

out_mat <- BRCAMatrix %>% data.frame() %>% tibble::rownames_to_column("GeneID") %>% 
  separate(GeneID, into=c("Gene","Version")) %>% 
  dplyr::select(-Version) %>% 
  filter(!duplicated(Gene))
colnames(out_mat) <- c(paste(rep("Normal",5),1:5,sep="-"),paste(rep("Tumour",5),1:5,sep="-"))


library(limma)


raw <- read.csv("tcga_raw_counts.csv")
voom_mat <- voom(raw[,-1])$E
voom_mat[1:4,1:4]

## Try and filter to top 1000 genes by IQR

iqrs <- apply(voom_mat,1,IQR)
picked_rows <- order(iqrs, decreasing = TRUE)[1:3000]
colnames(voom_mat) <- gsub(".","_",colnames(voom_mat),fixed=TRUE)
voom_out <- t(voom_mat[picked_rows,])
colnames(voom_out) <- raw[picked_rows,1]
voom_out <- data.frame(Subject = colnames(voom_mat), voom_out)
write.csv(voom_out,file="tcga_wgcna_counts.csv",row.names = FALSE,quote=FALSE)

clin <- data.frame(Subject=colnames(voom_mat),Group = c(rep(1,5),rep(2,5)))
readr::write_tsv(clin, path = "tcga_clinical.tsv")
