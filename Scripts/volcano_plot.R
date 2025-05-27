#Created Date: 20220805
#Modified Date: 20250527
#Author: Hein Ko Oo

#Install and load packagegs in R studio
install.packages(ggpubr)
install.packages(ggplot2)
install.pacakges(ggrepel)
inatall.packages(writexl)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library("writexl")

#Import protein ratios and p values csv file into R studio
protein_expression=read.csv("C:/Users/Yourname/Documents/Protein_Ratios_AND_pValues.csv", header=T);head(data);dim(data)

#Create additional column for labeling significant proteins
protein_expression$Expression= "Unchaged"

#Define thresholds for significance:
# if logFC > 2 and pvalue < 0.05, set as "Upregulated"
protein_expression$Expression[protein_expression$Log2.Fold.Change > 1.0 & protein_expression$p.Value.Adj < 0.05]="Upregulated"
# if logFC < -2 and pvalue < 0.05, set as "Downregulated"
protein_expression$Expression[protein_expression$Log2.Fold.Change < -1.0 & protein_expression$p.Value.Adj < 0.05]="Downregulated"

#Convert protein "Uniprot ID" to ENSEMBL, ENTREZID, and SYMBOL ID 
colnames(protein_expression)[colnames(protein_expression) == "Accession"] = "UNIPROT"
sig_genes<- protein_expression$UNIPROT
sig_genes
protein_expression.df <- bitr(sig_genes, fromType = "UNIPROT", toType = c("ENSEMBL","ENTREZID", "SYMBOL"),
                             OrgDb = org.Mm.eg.db)

protein_expression_all_IDs=merge(protein_expression,proteome_expression.df, by=c("UNIPROT"))

# Create a new column "delabel" in protein_expression table, that will contain the gene symbol of proteins diferentially expressed (NA in case they are not).
protein_expression_all_IDs$delabel <- NA
protein_expression_all_IDs$delabel[protein_expression_all_IDs$Expression != "Unchaged"] <- protein_expression_all_IDs$Symbol[protein_expression_all_IDs$Expression !="Unchaged"]

#Save the protein expression table including multiple gene IDs
write_xlsx(protein_expression_all_IDs,"C:/Users/Yourname/protein_expression_all_IDs_merge.xlsx")

#Plot Volcano plot and highlight upregulate and downregulate proteins
ggplot(data=protein_expression_all_IDs, aes(x=Log2.Fold.Change, y=-log10(p.Value.Adj), col=Expression, label=delabel)) +
  geom_point() + 
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p.Adj"))+
  geom_text_repel() +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_vline(xintercept=c(-2.0, 2.0), col="black", lty=2) + #Change according to the logFC value defined
  geom_hline(yintercept=-log10(0.05), col="black",lty=2)

