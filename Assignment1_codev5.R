### LSM3241 Assignment 1
## data: GSE4498
## matric no: A0099537L

setwd("GSE4498")
## load libraries
library(affy) # to read affymetrix data
library(limma) # linear modelling
library(magrittr)
library(hgu133plus2.db) # database used in study
library(biomaRt) # biomart query for later

## Read the affymetrix data and input it into R for analysis!
affyData = ReadAffy()

## normalization of affymetrix data
RMAnorm = rma(affyData)
MASnorm = mas5(affyData)
# the following step is essential since FC for MAS isn't log transformed (unlike RMA)
MAS.log = exprs(MASnorm)
MAS.log = log(MAS.log, 2)

## creating the model matrix + contrast matrix
# this assigns the sample into groups + sets comparison
# note: 22 samples: 12 non-smokers, 10 smokers
smokingstatus = c(rep("nonsmoker", 12), rep("smoker", 10))
modelmatrix = model.matrix(~smokingstatus+0)
colnames(modelmatrix) = c("nonsmoker", "smoker")
contrastmatrix = makeContrasts(smoker-nonsmoker, levels = modelmatrix)

## fitting affyData into linear model
linfit1 = lmFit(RMAnorm, modelmatrix)
linfit2 = lmFit(MASnorm, modelmatrix)
contrastfit1 = contrasts.fit(linfit1, contrastmatrix) # makes use of contrastmatrix to compute contrasts from linear fit model of affyData
contrastfit1 = eBayes(contrastfit1) # computes empirical bayes statistics
contrastfit2 = contrasts.fit(linfit2, contrastmatrix) # repeat for MAS5 normalization
contrastfit2 = eBayes(contrastfit2)

## identifying DGE
significantRMA = topTable(contrastfit1, p.value = 0.05, lfc = 1.5, n = Inf) # making use of the criteria used by the authors to determine which genes had significant changes in expression
significantMAS = topTable(contrastfit2, p.value = 0.05, lfc = 1.5, n = Inf)

# recording DGE results
write.table(significantRMA, file = "significantRMA.csv", quote = F, sep = "\t") # to an excel-readable file.
write.table(significantMAS, file = "significantMAS.csv", quote = F, sep = "\t")

## gene annotations
probeMappingA = select(hgu133plus2.db, columns = c("SYMBOL", "GO"), keys = rownames(significantRMA))
probeMappingB = select(hgu133plus2.db, columns = c("SYMBOL", "GO"), keys = rownames(significantMAS))

# note: this is for exporting data
write.table(probeMappingA, file = "probeMappingRMA.csv", quote = F, sep = "\t")
write.table(probeMappingB, file = "probeMappingMAS.csv", quote = F, sep = "\t")

## to get more data regarding genes upregulated/downregulated, we can use biomaRt!
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "asia.ensembl.org")
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
RMAmart = getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position",
                                "end_position", "ensembl_gene_id", "percentage_gc_content"),
                 filters = "hgnc_symbol", values = probeMappingA$SYMBOL,
                 mart = ensembl)
MASmart = getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position",
                                "end_position", "ensembl_gene_id", "percentage_gc_content"),
                 filters = "hgnc_symbol", values = probeMappingB$SYMBOL,
                 mart = ensembl)

# recording biomart search query results
write.table(RMAmart, file = "biomart_RMA.csv", quote = F, sep = "\t")
write.table(MASmart, file = "biomart_MAS.csv", quote = F, sep = "\t")