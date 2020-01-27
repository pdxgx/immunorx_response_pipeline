# Load relevant libraries
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(pROC)
library(RColorBrewer)
library(scales)

# Set working directory
setwd("YOUR DIRECTORY HERE") ## UPDATE THIS DIRECTORY 

##### Data loading/pre-processing #####

# Load data
immunorx <- read.table('immunotherapy_data_table.tsv', sep='\t', header=T, row.names = NULL)

# Add up indels
immunorx$Total_indels <- immunorx$Inframe_insertions + immunorx$Inframe_deletions + immunorx$Frameshift_deletions + immunorx$Frameshift_insertions
immunorx$Frameshift_indels <- immunorx$Frameshift_deletions + immunorx$Frameshift_insertions
immunorx$Inframe_indels <- immunorx$Inframe_insertions + immunorx$Inframe_deletions
immunorx$Insertions <- immunorx$Inframe_insertions + immunorx$Frameshift_insertions
immunorx$Deletions <- immunorx$Inframe_deletions + immunorx$Frameshift_deletions

# Adjust mutational burdens by coverage
immunorx$Total_mutations_adj <- immunorx$Total_mutations/immunorx$Coverage
immunorx$SNVs_adj <- immunorx$SNVs/immunorx$Coverage
immunorx$Inframe_ins_adj <- immunorx$Inframe_insertions/immunorx$Coverage
immunorx$Inframe_del_adj <- immunorx$Inframe_deletions/immunorx$Coverage
immunorx$Frameshift_ins_adj <- immunorx$Frameshift_insertions/immunorx$Coverage
immunorx$Frameshift_del_adj <- immunorx$Frameshift_deletions/immunorx$Coverage
immunorx$Total_indels_adj <- immunorx$Total_indels/immunorx$Coverage
immunorx$Frameshift_indels_adj <- immunorx$Frameshift_indels/immunorx$Coverage
immunorx$Inframe_indels_adj <- immunorx$Inframe_indels/immunorx$Coverage
immunorx$Insertions_adj <- immunorx$Insertions/immunorx$Coverage
immunorx$Deletions_adj <- immunorx$Deletions/immunorx$Coverage
immunorx$Muse_variants_adj <- immunorx$Muse_variants/immunorx$Coverage
immunorx$Mutect_variants_adj <- immunorx$Mutect_variants/immunorx$Coverage
immunorx$Pindel_variants_adj <- immunorx$Pindel_variants/immunorx$Coverage
immunorx$Radia_variants_adj <- immunorx$Radia_variants/immunorx$Coverage
immunorx$Somaticsniper_variants_adj <- immunorx$Somaticsniper_variants/immunorx$Coverage
immunorx$Varscan_variants_adj <- immunorx$Varscan_variants/immunorx$Coverage
immunorx$AtoC_adj <- immunorx$AtoC/immunorx$Coverage
immunorx$AtoG_adj <- immunorx$AtoG/immunorx$Coverage
immunorx$AtoT_adj <- immunorx$AtoT/immunorx$Coverage
immunorx$CtoA_adj <- immunorx$CtoA/immunorx$Coverage
immunorx$CtoG_adj <- immunorx$CtoG/immunorx$Coverage
immunorx$CtoT_adj <- immunorx$CtoT/immunorx$Coverage
immunorx$Synonymous_SNVs_adj <- immunorx$Synonymous_SNVs/immunorx$Coverage
immunorx$Nonsynonymous_SNVs_adj <- immunorx$Nonsynonymous_SNVs/immunorx$Coverage
immunorx$Start_codon_mutations_adj <- immunorx$Start_codon_mutations/immunorx$Coverage
immunorx$Stop_codon_mutations_adj <- immunorx$Stop_codon_mutations/immunorx$Coverage

# Create data tables for each disease
melanoma <- immunorx[immunorx$Disease == 'melanoma',]
RCC <- immunorx[immunorx$Disease == 'RCC',]
NSCLC <- immunorx[immunorx$Disease == 'NSCLC',]
MMR <- subset(immunorx, immunorx$Disease == 'colon' | immunorx$Disease == 'endometrial' | immunorx$Disease == 'thyroid')
prostate <- immunorx[immunorx$Disease == 'prostate',]

##### Mutation/neoepitope burdens #####

# Find median mutation burdens
median(immunorx$Total_mutations)
median(immunorx$Total_mutations[immunorx$Disease == 'melanoma'])
median(immunorx$Total_mutations[immunorx$Disease == 'NSCLC'])
median(immunorx$Total_mutations[immunorx$Disease == 'RCC'])
median(immunorx$Total_mutations[immunorx$Disease == 'prostate'])
median(MMR$Total_mutations)
median(immunorx$Total_mutations_adj)
median(immunorx$Total_mutations_adj[immunorx$Disease == 'melanoma'])
median(immunorx$Total_mutations_adj[immunorx$Disease == 'NSCLC'])
median(immunorx$Total_mutations_adj[immunorx$Disease == 'RCC'])
median(immunorx$Total_mutations_adj[immunorx$Disease == 'prostate'])
median(MMR$Total_mutations_adj)

# Find average % of SNVs/indels
median(immunorx$SNVs/immunorx$Total_mutations)
median(immunorx$Total_indels[immunorx$Disease == 'melanoma']/immunorx$Total_mutations[immunorx$Disease == 'melanoma'])
median(immunorx$Total_indels[immunorx$Disease == 'NSCLC']/immunorx$Total_mutations[immunorx$Disease == 'NSCLC'])
median(immunorx$Total_indels[immunorx$Disease == 'RCC']/immunorx$Total_mutations[immunorx$Disease == 'RCC'])
median(immunorx$Total_indels[immunorx$Disease == 'prostate']/immunorx$Total_mutations[immunorx$Disease == 'prostate'])
median(MMR$Total_indels/MMR$Total_mutations) # 31.2

# Find median neoepitope burdens
median(immunorx$MHCnuggets_eps)
median(immunorx$MHCnuggets_ClassI_eps)
median(immunorx$MHCnuggets_ClassII_eps)

median(melanoma$MHCnuggets_eps)
median(melanoma$MHCnuggets_ClassI_eps)
median(melanoma$MHCnuggets_ClassII_eps)

median(MMR$MHCnuggets_eps)
median(MMR$MHCnuggets_ClassI_eps)
median(MMR$MHCnuggets_ClassII_eps)

median(NSCLC$MHCnuggets_eps)
median(NSCLC$MHCnuggets_ClassI_eps)
median(NSCLC$MHCnuggets_ClassII_eps)

median(RCC$MHCnuggets_eps)
median(RCC$MHCnuggets_ClassI_eps)
median(RCC$MHCnuggets_ClassII_eps)

median(prostate$MHCnuggets_eps)
median(prostate$MHCnuggets_ClassI_eps)
median(prostate$MHCnuggets_ClassII_eps)


##### Plotting burdens #####

# Function to plot raw mutational burden
plot_mutation_burden <- function(disease, start_pos, end_pos){
  for (i in 1:length(disease$Total_mutations)){
    if (disease$MSI_status[i] == 1){
      color <- 'red'
    } else {
      color <- 'black'
    }
    points(i*1.75/length(disease$Total_mutations)+start_pos+0.1,log10(disease$Total_mutations[i]),pch=20, col=color)
  }
  med_val <- median(log10(disease$Total_mutations))
  points(c(start_pos+0.2,end_pos-0.2),c(med_val,med_val),col="firebrick",type="l")  
}

# Function to plot coverage-adjusted mutational burden
plot_adj_burden <- function(disease, start_pos, end_pos){
  for (i in 1:length(disease$Total_mutations_adj)){
    if (disease$MSI_status[i] == 1){
      color <- 'red'
    } else {
      color <- 'black'
    }
    points(i*1.75/length(disease$Total_mutations_adj)+start_pos+0.1,log10(disease$Total_mutations_adj[i]),pch=20,col=color)
  }
  med_val <- median(log10(disease$Total_mutations_adj))
  points(c(start_pos+0.2,end_pos-0.2),c(med_val,med_val),col="firebrick",type="l")  
}

# Function to plot FS indel burden
plot_fs_indel_burden <- function(disease, start_pos, end_pos){
  for (i in 1:length(disease$Frameshift_indels)){
    if (disease$MSI_status[i] == 1){
      color <- 'red'
    } else {
      color <- 'black'
    }
    points(i*1.75/length(disease$Frameshift_indels)+start_pos+0.1,log10(disease$Frameshift_indels[i]),
           pch=20, col=color)
  }
  med_val <- median(log10(disease$Frameshift_indels))
  points(c(start_pos+0.2,end_pos-0.2),c(med_val,med_val),col="firebrick",type="l")  
}

# Function to plot deletion burden
plot_inframe_indel_burden <- function(disease, start_pos, end_pos){
  for (i in 1:length(disease$Inframe_indels)){
    if (disease$MSI_status[i] == 1){
      color <- 'red'
    } else {
      color <- 'black'
    }
    points(i*1.75/length(disease$Inframe_indels)+start_pos+0.1,log10(disease$Inframe_indels[i]),
           pch=20, col=color)
  }
  med_val <- median(log10(disease$Inframe_indels))
  points(c(start_pos+0.2,end_pos-0.2),c(med_val,med_val),col="firebrick",type="l")  
}

# Function to plot neoepitope burden
plot_binding_burden <- function(disease, start_pos, end_pos){
  for (i in 1:length(disease$MHCnuggets_eps)){
    if (disease$MSI_status[i] == 1){
      color <- 'red'
    } else {
      color <- 'black'
    }
    points(i*1.75/length(disease$MHCnuggets_eps)+start_pos+0.1,log10(disease$MHCnuggets_eps[i]),
           pch=20, col=color)
  }
  med_val <- median(log10(disease$MHCnuggets_eps))
  points(c(start_pos+0.2,end_pos-0.2),c(med_val,med_val),col="firebrick",type="l")  
}

# Layout combined adjusted mutation burden/neoepitope burden plot - Figure 1
layout(matrix(c(1,2), nrow = 1, byrow=TRUE))
par(mar=c(4,5,0.5,1))
# Coverage - adjusted mutations
MMR <- MMR[order(MMR$Total_mutations_adj),]
melanoma <- melanoma[order(melanoma$Total_mutations_adj),]
NSCLC <- NSCLC[order(NSCLC$Total_mutations_adj),]
RCC <- RCC[order(RCC$Total_mutations_adj),]
prostate <- prostate[order(prostate$Total_mutations_adj),]
plot(0,0,xlim=c(0.2,9.8),ylim=c(-1,3),type="n",axes=F,xlab="",
     ylab="Mutation burden (# somatic variants/Mbp)", cex.lab=1.75)
axis(1, at = seq(1, 9, 2), labels = c("Prostate\nn=10", "NSCLC\nn=34", "RCC\nn=57",  
                                      "Melanoma\nn=285", "MMR\ndeficient\nn=28"), 
     cex.axis=0.7, tick=F, padj=0.5)
axis(2, at = seq(-1,3,1), labels = c(0.1,1,10,100,1000), cex.axis=1)
rect(xleft=2, xright=4, ytop=6, ybottom=-10, col="lightgray", border=NA)
rect(xleft=6, xright=8, ytop=6, ybottom=-10, col="lightgray", border=NA)
box(lty="solid", col="black")
mtext("A", side=3, line=-1.2, at=-3, cex=2)
plot_adj_burden(prostate,0,2)
plot_adj_burden(NSCLC,2,4)
plot_adj_burden(RCC,4,6)
plot_adj_burden(melanoma,6,8)
plot_adj_burden(MMR,8,10)
# Neoepitope
MMR <- MMR[order(MMR$MHCnuggets_eps),]
melanoma <- melanoma[order(melanoma$MHCnuggets_eps),]
NSCLC <- NSCLC[order(NSCLC$MHCnuggets_eps),]
RCC <- RCC[order(RCC$MHCnuggets_eps),]
prostate <- prostate[order(prostate$MHCnuggets_eps),]
plot(0,0,xlim=c(0.2,9.8),ylim=c(0,6.5),type="n",axes=F,xlab="",
     ylab="Neoepitope burden (#)", cex.lab=1.75)
axis(1, at = seq(1, 9, 2), labels = c( "NSCLC\nn=34", "Prostate\nn=10", 
                                      "Melanoma\nn=285", "RCC\nn=57", "MMR\ndeficient\nn=28"),
     cex.axis=0.7, tick=F, padj=0.5)
axis(2, at = seq(1,6,1), labels = c("10","100","1000","10000","100000","1000000"), cex.axis=1)
rect(xleft=2, xright=4, ytop=8, ybottom=-1, col="lightgray", border=NA)
rect(xleft=6, xright=8, ytop=8, ybottom=-1, col="lightgray", border=NA)
box(lty="solid", col="black")
mtext("B", side=3, line=-1.2, at=-3, cex=2)
plot_binding_burden(NSCLC,0,2)
plot_binding_burden(prostate,2,4)
plot_binding_burden(melanoma,4,6)
plot_binding_burden(RCC,6,8)
plot_binding_burden(MMR,8,10)
layout(matrix(c(1), nrow = 1, byrow=TRUE))

# Plot raw mutational burden - Supplementary Figure 2
MMR <- MMR[order(MMR$Total_mutations),]
melanoma <- melanoma[order(melanoma$Total_mutations),]
NSCLC <- NSCLC[order(NSCLC$Total_mutations),]
RCC <- RCC[order(RCC$Total_mutations),]
prostate <- prostate[order(prostate$Total_mutations),]
par(mar=c(5,7,2,5))
plot(0,0,xlim=c(0.2,9.8),ylim=c(1.7,4.9),type="n",axes=F,xlab="",ylab="Total mutation burden (raw #)", 
     cex.lab=2)
axis(1, at = seq(1, 9, 2), labels = c("Prostate\nn=10", "NSCLC\nn=34", 
                                      "Melanoma\nn=302", "RCC\nn=57", "MMR\ndeficient\nn=28"),
     cex.axis=1.4, tick=F, padj=0.7)
axis(2, at = seq(0,5,1), labels = c("0","10","100","1000","10000","100000"), cex.axis=1.5)
rect(xleft=2, xright=4, ytop=6, ybottom=-1, col="lightgray", border=NA)
rect(xleft=6, xright=8, ytop=6, ybottom=-1, col="lightgray", border=NA)
box(lty="solid", col="black")
plot_mutation_burden(prostate,0,2)
plot_mutation_burden(NSCLC,2,4)
plot_mutation_burden(melanoma,4,6)
plot_mutation_burden(RCC,6,8)
plot_mutation_burden(MMR,8,10)

# Plot indel burdens - Supplementary Figure 3
# First FS
layout(matrix(c(1,2), nrow = 1, byrow=TRUE))
par(mar=c(4,4.5,0.5,1))
MMR <- MMR[order(MMR$Frameshift_indels),]
melanoma <- melanoma[order(melanoma$Frameshift_indels),]
NSCLC <- NSCLC[order(NSCLC$Frameshift_indels),]
RCC <- RCC[order(RCC$Frameshift_indels),]
prostate <- prostate[order(prostate$Frameshift_indels),]
plot(0,0,xlim=c(0.2,9.8),ylim=c(0,5),type="n",axes=F,xlab="",ylab="Total FS indel burden (raw #)", 
     cex.lab=1.75)
axis(1, at = seq(1, 9, 2), labels = c("Prostate\nn=10", "RCC\nn=57", "NSCLC\nn=34", 
                                      "Melanoma\nn=285", "MMR\ndeficient\nn=28"),
     cex.axis=0.7, tick=F, padj=0.7)
axis(2, at = seq(0,5,1), labels = seq(0,5,1), cex.axis=1)
rect(xleft=2, xright=4, ytop=6, ybottom=-1, col="lightgray", border=NA)
rect(xleft=6, xright=8, ytop=6, ybottom=-1, col="lightgray", border=NA)
box(lty="solid", col="black")
mtext("A", side=3, line=-1.2, at=-2.5, cex=2)
plot_fs_indel_burden(NSCLC,0,2)
plot_fs_indel_burden(prostate,2,4)
plot_fs_indel_burden(RCC,4,6)
plot_fs_indel_burden(melanoma,6,8)
plot_fs_indel_burden(MMR,8,10)
# Then in-frame
MMR <- MMR[order(MMR$Inframe_indels),]
melanoma <- melanoma[order(melanoma$Inframe_indels),]
NSCLC <- NSCLC[order(NSCLC$Inframe_indels),]
RCC <- RCC[order(RCC$Inframe_indels),]
prostate <- prostate[order(prostate$Inframe_indels),]
plot(0,0,xlim=c(0.2,9.8),ylim=c(0,5),type="n",axes=F,xlab="",ylab="Total in-frame indel burden (raw #)", 
     cex.lab=1.75)
axis(1, at = seq(1, 9, 2), labels = c("Prostate\nn=10", "RCC\nn=57", "NSCLC\nn=34", 
                                      "Melanoma\nn=285", "MMR\ndeficient\nn=28"),
     cex.axis=0.7, tick=F, padj=0.7)
axis(2, at = seq(0,5,1), labels = seq(0,5,1), cex.axis=1)
rect(xleft=2, xright=4, ytop=6, ybottom=-1, col="lightgray", border=NA)
rect(xleft=6, xright=8, ytop=6, ybottom=-1, col="lightgray", border=NA)
box(lty="solid", col="black")
mtext("B", side=3, line=-1.2, at=-2.5, cex=2)
plot_inframe_indel_burden(prostate,0,2)
plot_inframe_indel_burden(NSCLC,2,4)
plot_inframe_indel_burden(melanoma,4,6)
plot_inframe_indel_burden(RCC,6,8)
plot_inframe_indel_burden(MMR,8,10)
layout(matrix(c(1), nrow = 1, byrow=TRUE))


##### Comparing mutational burden and neoepitope burden #####

# Correlation between burdens - Supplementary Figure 4
x <- immunorx$Total_mutations
y <- immunorx$MHCnuggets_eps
summary(rlm(y~x)) # Get slope/intercept
par(mar=c(4.5,5,0.75,1))
plot(x, y, ylim=c(0, 2800000), xlim=c(0,73600), 
     xlab="Total mutation burden (raw #)", ylab="Neoepitope burden (raw #)", 
     cex.lab=2, cex.axis=1.75, pch=19, cex=0.5, axes=F)
axis(1, at = seq(0, 60000, 20000), labels = seq(0, 60000, 20000), cex.axis=1.75)
axis(2, at = seq(0,2500000,500000), 
     labels = c("0","500000","1000000","1500000", "2000000", "2500000"), cex.axis=1.75)
abline(rlm(y~x), col='red')
text(65000, -700, labels="y = 5.82x - 619.59",cex=1.5, col="black") # adjust for slope/intercept
box(lty="solid", col="black")

# Get Pearson correlation
cor.test(x,y)

# Set up panel.cor
panel.cor <- function(x, y, digits=2, prefix="") 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex=1.25) 
  text(.8, .8, Signif, cex=1.25, col=2) 
}

# Corplot for mutational burden across callers - Supplementary Figure 15
pairs(cbind(immunorx$Total_mutations, immunorx$Muse_variants, immunorx$Mutect_variants, immunorx$Radia_variants,
            immunorx$Somaticsniper_variants, immunorx$Varscan_variants, immunorx$Pindel_variants), 
      lower.panel=panel.smooth, upper.panel=panel.cor,
      labels = c('Consensus\ncalls', 'MuSE', 'MuTect', 'RADIA', 'SomaticSniper', 'VarScan', 'Pindel'),
      cex.labels = 1.25)

# Variance across callers
caller_variance <- rep(0, length(immunorx$Patient))
max_diff <- rep(0, length(immunorx$Patient))
for (i in 1:length(immunorx$Patient)){
  counts <- c(immunorx$Muse_variants[i], immunorx$Mutect_variants[i], immunorx$Pindel_variants[i], 
              immunorx$Radia_variants[i], immunorx$Somaticsniper_variants[i], immunorx$Varscan_variants[i])
  caller_variance[i] <- mad(counts)
  max_diff[i] <- max(counts) - min(counts)
}
immunorx$Caller_variation <- caller_variance
immunorx$Max_TMB_diff <- max_diff
immunorx$Percent_TMB_diff <- immunorx$Max_TMB_diff/immunorx$Total_mutations

median(immunorx$Caller_variation)
median(immunorx$Max_TMB_diff)
median(immunorx$Percent_TMB_diff)

# Plot variance across callers - Supplementary Figure 16
x <- immunorx$Total_mutations
y <- immunorx$Caller_variation
summary(lm(y~x)) # Get slope/intercept
par(mar=c(4.5,5,0.75,1))
plot(x, y, 
     xlab="Somatic mutations (consensus calls, raw #)", ylab="Variation in mutation count (med. absolute deviation)", 
     cex.lab=1.7, cex.axis=1.25, pch=19, cex=0.5)
abline(lm(y~x), col='red')
text(66000, -300, labels="y = 0.15x + 108.1",cex=1.5, col="black") # adjust for slope/intercept

cor.test(x,y)


##### Identity of neoepitopes across callers #####

# Load upset plot data
df <- read.table('epitope_upset.tsv', sep='\t', header=T)
df$overlap <- factor(df$overlap, levels=c("Consensus", "RADIA", "MuTect", "Pindel", "SomaticSniper", 
                                          "VarScan", "MuSE", "Pindel_Consensus", "MuTect_SomaticSniper", 
                                          "MuSE_MuTect_Consensus", "SomaticSniper_VarScan_Consensus",
                                          "MuSE_MuTect_RADIA_Consensus", "RADIA_SomaticSniper_VarScan_Consensus",
                                          "MuSE_MuTect_RADIA_VarScan_Consensus", 
                                          "MuSE_MuTect_SomaticSniper_VarScan_Consensus",
                                          "MuSE_MuTect_RADIA_SomaticSniper_VarScan_Consensus",
                                          "MuSE_MuTect_Pindel_RADIA_VarScan_Consensus",
                                          "MuSE_MuTect_Pindel_RADIA_SomaticSniper_VarScan_Consensus", "Other"
                                          )
)
df$tool <- factor(df$tool, levels=c('MuSE', 'MuTect', 'Pindel', 'RADIA', 'SomaticSniper', 'VarScan', 'Consensus'))
df$classification <- factor(df$classification, levels=c('present', 'absent', 'other'))
df <- df[df$overlap != 'Other',]
df$overlap <- factor(df$overlap, levels=c("Consensus", "RADIA", "MuTect", "Pindel", "SomaticSniper", 
                                          "VarScan", "MuSE", "Pindel_Consensus", "MuTect_SomaticSniper", 
                                          "MuSE_MuTect_Consensus", "SomaticSniper_VarScan_Consensus",
                                          "MuSE_MuTect_RADIA_Consensus", "RADIA_SomaticSniper_VarScan_Consensus",
                                          "MuSE_MuTect_RADIA_VarScan_Consensus", 
                                          "MuSE_MuTect_SomaticSniper_VarScan_Consensus",
                                          "MuSE_MuTect_RADIA_SomaticSniper_VarScan_Consensus",
                                          "MuSE_MuTect_Pindel_RADIA_VarScan_Consensus",
                                          "MuSE_MuTect_Pindel_RADIA_SomaticSniper_VarScan_Consensus", "Other"
                                          )
)

# Set up data for pie chart
p <- ggplot(df, aes(x="", y=proportion, fill=classification))+geom_bar(width = 1, stat = 'identity')+theme(axis.text = element_blank(),
                                                                                                           axis.ticks = element_blank(),
                                                                                                           panel.grid  = element_blank())
pie <- p + coord_polar('y',start = 0) + scale_fill_manual(values=c("black", "white", "white"))

# Facet with tools (rows) and columns (overlap in recognition)
final_graph <- pie+facet_grid(rows = vars(tool), labeller=label_wrap_gen(width=18), 
                              cols = vars(overlap), switch='y')+theme(axis.text = element_blank(),
                                                                      axis.ticks = element_blank(),
                                                                      panel.grid  = element_blank(),
                                                                      axis.title.x = element_blank(),
                                                                      axis.title.y = element_blank(),
                                                                      strip.text.x = element_blank(),
                                                                      strip.text.y = element_text(angle=180, 
                                                                                                  size = rel(1.5)),
                                                                      strip.background = element_rect(fill = 'white'),
                                                                      legend.position = "none",
                                                                      panel.background = element_rect(fill="white", 
                                                                                                      colour = 'white'),
                                                                      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                                                                      panel.spacing=unit(0, "lines"),
                                                                      plot.margin=unit(c(-5,1.2,0,0), 'cm'))

# Load count data for bar graph
df2 <- read.table('epitope_upset_counts.tsv', sep='\t', header=T)
df2$overlap <- factor(df2$overlap, levels=c("Consensus", "RADIA", "MuTect", "Pindel", "SomaticSniper", 
                                            "VarScan", "MuSE", "Pindel_Consensus", "MuTect_SomaticSniper", 
                                            "MuSE_MuTect_Consensus", "SomaticSniper_VarScan_Consensus",
                                            "MuSE_MuTect_RADIA_Consensus", "RADIA_SomaticSniper_VarScan_Consensus",
                                            "MuSE_MuTect_RADIA_VarScan_Consensus", 
                                            "MuSE_MuTect_SomaticSniper_VarScan_Consensus",
                                            "MuSE_MuTect_RADIA_SomaticSniper_VarScan_Consensus",
                                            "MuSE_MuTect_Pindel_RADIA_VarScan_Consensus",
                                            "MuSE_MuTect_Pindel_RADIA_SomaticSniper_VarScan_Consensus", "Other"))
df2$count_name <- as.character(df2$count)

# Create bar plot color palette
color_pallette <- brewer.pal(7, "Set2")
full_color_pallette <- c(rep("gray", 7))
for (i in 1:5){
  full_color_pallette <- append(full_color_pallette, rep(color_pallette[i], 2))
}
full_color_pallette <- append(full_color_pallette, c(color_pallette[6], color_pallette[7]))

# Create bar plot
p2 <- ggplot(df2, 
             aes(x=overlap, 
                 y=count))+geom_bar(stat = 'identity', 
                                    fill = full_color_pallette, 
                                    color='black')+geom_text(aes(label=count_name), 
                                                             hjust=0.5,
                                                             vjust=-0.2,
                                                             size = rel(2.5),
                                                             color="black")+theme(axis.text.x = element_blank(),
                                                                                  panel.grid  = element_blank(),
                                                                                  axis.title.x = element_blank(),
                                                                                  axis.text = element_text(size = rel(1.5)),
                                                                                  axis.title.y = element_text(size = rel(1.5),
                                                                                  margin = margin(t = 0, r = 15, b = 0, l = 0)),
                                                                                  panel.background = element_rect(fill="white", colour = 'white'),
                                                                                  axis.line.y.left = element_line(colour = 'black'),
                                                                                  axis.ticks.x = element_blank(),
                                                                                   plot.margin=unit(c(2,0.2,0,1.1), 'cm'))+scale_y_continuous(trans = log10_trans(),
                                                                                                                                              breaks = trans_breaks("log10", function(x) 10^x),
                                                                                                                                              labels = trans_format("log10", math_format(10^.x)))+annotate("text", 
                                                                                                                                                                                                           label = "*", 
                                                                                                                                                                                                           x = 19, 
                                                                                                                                                                                                           y = 8e5, 
                                                                                                                                                                                                           size=rel(10))

# Arrange plots - Supplementary Figure 19
grid.arrange(p2+labs(y = "Neopeptide sequences (#)"), final_graph, layout_matrix=rbind(c(1), c(2)))


##### Driver mutation analysis #####

drivers <- read.table('drivers.tsv', sep='\t', header=T)
drivers$Percent_presented <- drivers$Presented_clinvar_variants/drivers$Total_clinvar_variants
drivers$Percent_eps_presented <- drivers$Presented_clinvar_epitopes/drivers$Total_clinvar_neopeptides

length(drivers$Total_clinvar_variants[drivers$Total_clinvar_variants > 0])
mean(drivers$Percent_presented[drivers$Total_clinvar_variants != 0])
mean(drivers$Percent_eps_presented[drivers$Total_clinvar_neopeptides != 0])

mean(drivers$Total_clinvar_variants)
mean(drivers$Total_clinvar_variants[drivers$Total_clinvar_variants > 0])




##### HLA presentation trends #####

ep_binders <- read.table('ep_binders.tsv', sep='\t', header=T)
mut_binders <- read.table('mut_binders.tsv', sep='\t', header=T)

ep_drivers <- read.table('driver_ep_binders.tsv', sep='\t', header=T)
mut_drivers <- read.table('driver_mut_binders.tsv', sep='\t', header=T)

immunorx$Eps_presented <- immunorx$MHCnuggets_eps/immunorx$Total_comprehensive_neoepitopes
median(immunorx$Eps_presented)

length(ep_binders$Epitope[ep_binders$Allele_count > 0])/length(ep_binders$Epitope) # 0.1135

# Find HLA heterozygosity
allele_count <- rep(NA, length(mut_binders$Patient))
for (i in 1:length(mut_binders$Patient)){
  pat <- as.character(mut_binders$Patient[i])
  hla1 <- immunorx$Tumor_HLA1_count[immunorx$Patient == pat]
  hla2 <- immunorx$Tumor_HLA2_count[immunorx$Patient == pat]
  allele_count[i] <- hla1 + hla2
}
mut_binders$HLA_count <- as.factor(allele_count)

allele_count <- rep(NA, length(mut_drivers$Patient))
for (i in 1:length(mut_drivers$Patient)){
  pat <- as.character(mut_drivers$Patient[i])
  hla1 <- immunorx$Tumor_HLA1_count[immunorx$Patient == pat]
  hla2 <- immunorx$Tumor_HLA2_count[immunorx$Patient == pat]
  allele_count[i] <- hla1 + hla2
}
mut_drivers$HLA_count <- as.factor(allele_count)

# Correlation plots
muts <- hist(mut_binders$Allele_count, breaks=0:9, plot=F)$counts/length(mut_binders$Allele_count)
drivers <- hist(mut_drivers$Allele_count, breaks=0:9, plot=F)$counts/length(mut_drivers$Allele_count)
counter1 <- c(0:8, 0:8)
full_muts <- c(muts, drivers)
n <- c(hist(mut_binders$Allele_count, breaks=0:9, plot=F)$counts, hist(mut_drivers$Allele_count, breaks=0:9, plot=F)$counts)
categories <- c(rep('All', length(muts)), rep('Driver', length(drivers)))
df1 <- as.data.frame(cbind(counter1, full_muts, categories, n))
colnames(df1) <- c('Count', 'Alleles', 'Category', 'N')
df1$Alleles <- as.numeric(as.character(df1$Alleles))
df1$Count <- as.numeric(as.character(df1$Count))
df1$N <- as.numeric(as.character(df1$N))
df1$Category <- factor(df1$Category, levels = rev(levels(df1$Category)))

eps <- hist(ep_binders$Allele_count[ep_binders$Allele_count > 0], breaks=1:8, plot=F)$counts/length(ep_binders$Allele_count[ep_binders$Allele_count > 0])
driver_eps <- hist(ep_drivers$Allele_count[ep_drivers$Allele_count > 0], breaks=1:8, plot=F)$counts/length(ep_drivers$Allele_count[ep_drivers$Allele_count > 0])
full_eps <- c(eps, driver_eps)
n2 <- c(hist(ep_binders$Allele_count[ep_binders$Allele_count > 0], breaks=1:8, plot=F)$counts, hist(ep_drivers$Allele_count[ep_drivers$Allele_count > 0], breaks=1:8, plot=F)$counts)
counter2 <- c(1:7, 1:7)
categories <- c(rep('All', length(eps)), rep('Driver', length(driver_eps)))
df2 <- as.data.frame(cbind(counter2, full_eps, categories, n2))
colnames(df2) <- c('Count', 'Alleles', 'Category', 'N')
df2$Alleles <- as.numeric(as.character(df2$Alleles))
df2$Count <- as.numeric(as.character(df2$Count))
df2$N <- as.numeric(as.character(df2$N))
df2$Category <- factor(df2$Category, levels = rev(levels(df2$Category)))


p0 <- ggplot(df1, 
             aes(x=Count,
                 y=Alleles, shape=Category, colour=Category, fill=Category))+geom_point(aes(size = N), shape=1)+geom_smooth(method='loess',
                                                      formula=y~x)+scale_x_continuous(breaks=c(0:8))+scale_color_manual(values=c('blue', 'red'))+scale_fill_manual(values=c('#00BFC4', '#F8766D'))
p1 <- p0+labs(title="",
              x="HLA Alleles (#)",
              y = "Mutations (% of total)")+scale_y_continuous(labels=percent_format(accuracy = .0001, 
                                                                                     suffix = ''))+theme(panel.grid  = element_blank(),
                                                                                                         axis.title.x = element_text(size = rel(1.08)),
                                                                                                         axis.text = element_text(size = rel(1.08)),
                                                                                                         axis.title.y = element_text(size = rel(1.08)),
                                                                                                         panel.background = element_rect(fill="white", colour = 'white'),
                                                                                                         axis.line.y.left = element_line(colour = 'black'),
                                                                                                         axis.line.x.bottom = element_line(colour = 'black'),
                                                                                                         axis.line.x.top = element_line(colour = 'black'),
                                                                                                         axis.line.y.right = element_line(colour = 'black'),
                                                                                                         plot.tag = element_text(size = rel(1.5)),
                                                                                                         legend.position = "bottom",
                                                                                                         plot.margin=unit(c(0.1,0.1,0.2,0.1), 'cm'))+ labs(tag = "B")

p2 <- ggplot(df2, 
             aes(x=Count, 
                 y=Alleles, shape=Category, colour=Category, fill=Category))+geom_point(aes(size = N), shape=1)+scale_size(trans="log10", breaks=c(1,100,10000,1000000), labels=c("1","100","10000","1000000"))+geom_smooth(method='lm',
                                                      formula=y~x,
                                                      fullrange=T)+scale_color_manual(values=c('blue', 'red'))+scale_fill_manual(values=c('#00BFC4', '#F8766D'))
p3 <- p2+scale_y_continuous(trans = log10_trans(), 
                            labels=percent_format(accuracy = .0001, 
                                                  suffix = ''))+scale_x_continuous(breaks=c(1:7))+labs(title="",
                                                                                                       x="HLA Alleles (#)",
                                                                                                       y = "Neoepitopes (% of total)")+theme(panel.grid  = element_blank(),
                                                                                                                                             axis.title.x = element_text(size = rel(1.08)),
                                                                                                                                             axis.text = element_text(size = rel(1.08)),
                                                                                                                                             axis.title.y = element_text(size = rel(1.08)),
                                                                                                                                             panel.background = element_rect(fill="white", 
                                                                                                                                                                             colour = 'white'),
                                                                                                                                             axis.line.y.left = element_line(colour = 'black'),
                                                                                                                                             axis.line.x.bottom = element_line(colour = 'black'),
                                                                                                                                             axis.line.x.top = element_line(colour = 'black'),
                                                                                                                                             axis.line.y.right = element_line(colour = 'black'),
                                                                                                                                             plot.tag = element_text(size = rel(1.5)),
                                                                                                                                             plot.margin=unit(c(0.1,0.1,0.2,0.1), 'cm'))+ labs(tag = "A")

hla_levels <- levels(mut_binders$HLA_count)
patient_count <- rep(NA, length(hla_levels))
patient_count_drivers <- rep(NA, length(hla_levels))
med_alleles <- rep(NA, length(hla_levels))
med_alleles_driver <- rep(NA, length(hla_levels))
percent_presentable <- rep(NA, length(hla_levels))
percent_drivers <- rep(NA, length(hla_levels))
for (i in 1:length(hla_levels)){
  mut_df <- mut_binders[mut_binders$HLA_count == hla_levels[i],]
  med <- median(mut_df$Allele_count[mut_df$Allele_count > 0])
  perc <- length(mut_df$Allele_count[mut_df$Allele_count > 0])/length(mut_df$Allele_count)
  med_alleles[i] <- med
  percent_presentable[i] <- perc
  patient_count[i] <- length(unique(mut_df$Patient))
  driver_df <- mut_drivers[mut_drivers$HLA_count == hla_levels[i],]
  med <- median(driver_df$Allele_count[driver_df$Allele_count > 0])
  perc <- length(driver_df$Allele_count[driver_df$Allele_count > 0])/length(driver_df$Allele_count)
  med_alleles_driver[i] <- med
  percent_drivers[i] <- perc
  patient_count_drivers[i] <- length(unique(driver_df$Patient))
}

df3 <- as.data.frame(cbind(c(hla_levels, hla_levels), c(med_alleles, med_alleles_driver), 
                           c(percent_presentable, percent_drivers), 
                           c(rep('All', length(hla_levels)), rep('Driver', length(hla_levels))), 
                           c(patient_count, patient_count_drivers)))
colnames(df3) <- c('HLA_count', 'Median_number_presenting', 'Percent_variants_presentable', 'Category', 'N')
df3$HLA_count <- as.numeric(as.character(df3$HLA_count ))
df3$N <- as.numeric(as.character(df3$N ))
df3$Median_number_presenting <- as.numeric(as.character(df3$Median_number_presenting ))
df3$Percent_variants_presentable <- as.numeric(as.character(df3$Percent_variants_presentable))
df3$Category <- factor(df3$Category, levels = rev(levels(df3$Category)))

p4 <- ggplot(df3, 
             aes(x=HLA_count, 
                 y=Percent_variants_presentable, shape=Category, colour=Category, fill=Category))+geom_point(aes(size = N), shape=1)+geom_smooth(method='lm',
                                                                           formula=y~x, fullrange=T)+scale_color_manual(values=c('blue', 'red'))+scale_fill_manual(values=c('#00BFC4', '#F8766D'))
p5 <- p4+scale_y_continuous(labels=percent_format(accuracy = 0.1, 
                                                  suffix = ''))+scale_x_continuous(breaks=c(6:12))+labs(title="",
                                                                                                        x="Unique patient HLA Alleles (#)",
                                                                                                        y = "Presentable variants (% of total)")+theme(panel.grid  = element_blank(),
                                                                                                                                                       axis.title.x = element_text(size = rel(1.08)),
                                                                                                                                                       axis.text = element_text(size = rel(1.08)),
                                                                                                                                                       axis.title.y = element_text(size = rel(1.08)),
                                                                                                                                                       panel.background = element_rect(fill="white", 
                                                                                                                                                                                       colour = 'white'),
                                                                                                                                                       axis.line.y.left = element_line(colour = 'black'),
                                                                                                                                                       axis.line.x.bottom = element_line(colour = 'black'),
                                                                                                                                                       axis.line.x.top = element_line(colour = 'black'),
                                                                                                                                                       axis.line.y.right = element_line(colour = 'black'),
                                                                                                                                                       plot.tag = element_text(size = rel(1.5)),
                                                                                                                                                       plot.margin=unit(c(0.1,0.1,0.2,0.1), 'cm'))+ labs(tag = "C")

# Figure 3
grid.arrange(p3, p1, p5, layout_matrix=rbind(c(1), c(2), c(3)))

# Supplementary Figure 7
p6 <- ggplot(df3, 
             aes(x=HLA_count, 
                 y=Median_number_presenting))+geom_point()+geom_smooth(method='lm',
                                                                       formula=y~x, 
                                                                       color='red', fullrange=T)
p7 <- p6+scale_x_continuous(breaks=c(6:12))+labs(title="",
                                                 x="Unique patient HLA Alleles (#)",
                                                 y = "Median alleles presenting per variant (#)")+theme(panel.grid  = element_blank(),
                                                                                                        axis.title.x = element_text(size = rel(1.5)),
                                                                                                        axis.text = element_text(size = rel(1.25)),
                                                                                                        axis.title.y = element_text(size = rel(1.5)),
                                                                                                        panel.background = element_rect(fill="white", colour = 'white'),
                                                                                                        axis.line.y.left = element_line(colour = 'black'),
                                                                                                        axis.line.x.bottom = element_line(colour = 'black'),
                                                                                                        axis.line.x.top = element_line(colour = 'black'),
                                                                                                        axis.line.y.right = element_line(colour = 'black'),
                                                                                                        plot.tag = element_text(size = rel(1.5)),
                                                                                                        legend.position="bottom",
                                                                                                        legend.text = element_text(size = rel(1.25)),
                                                                                                        legend.title = element_text(size = rel(1.25)),
                                                                                                        plot.margin=unit(c(0.1,0.1,0.2,0.1), 'cm'))

p0 <- ggplot(mut_binders,
             aes(x=Disease,
                 y=Allele_count))+geom_violin(fill='gray')+labs(title="",
                                                                x="", 
                                                                y = "Alleles presenting per mutation (#)")+scale_x_discrete(limits=c("MMR-deficient",
                                                                                                                                     "melanoma",
                                                                                                                                     "prostate",
                                                                                                                                     "NSCLC",
                                                                                                                                     "RCC"))
p1 <- p0 + geom_boxplot(width=0.1, fill="white")+theme(panel.grid  = element_blank(),
                                                       axis.title.x = element_blank(),
                                                       axis.text = element_text(size = rel(1.25)),
                                                       axis.title.y = element_text(size = rel(1.5)),
                                                       panel.background = element_rect(fill="white", 
                                                                                       colour = 'white'),
                                                       axis.line.y.left = element_line(colour = 'black'),
                                                       axis.line.x.bottom = element_line(colour = 'black'),
                                                       axis.line.x.top = element_line(colour = 'black'),
                                                       axis.line.y.right = element_line(colour = 'black'),
                                                       legend.position = "none",
                                                       plot.tag = element_text(size = rel(1.5)),
                                                       plot.margin=unit(c(0.1,0.1,0.5,0.1), 'cm'))+ labs(tag = "B")

p2 <- ggplot(ep_binders[ep_binders$Allele_count > 0,],
             aes(x=Disease,
                 y=Allele_count))+geom_violin(fill='gray')+labs(title="",
                                                                x="", 
                                                                y = "Alleles presenting per epitope (#)")+scale_x_discrete(limits=c("MMR-deficient",
                                                                                                                                    "melanoma",
                                                                                                                                    "prostate",
                                                                                                                                    "NSCLC",
                                                                                                                                    "RCC"))
p3 <- p2 + geom_boxplot(width=0.1, fill="white")+theme(panel.grid  = element_blank(),
                                                       axis.title.x = element_blank(),
                                                       axis.text = element_text(size = rel(1.25)),
                                                       axis.title.y = element_text(size = rel(1.5)),
                                                       panel.background = element_rect(fill="white", 
                                                                                       colour = 'white'),
                                                       axis.line.y.left = element_line(colour = 'black'),
                                                       axis.line.x.bottom = element_line(colour = 'black'),
                                                       axis.line.x.top = element_line(colour = 'black'),
                                                       axis.line.y.right = element_line(colour = 'black'),
                                                       legend.position = "none",
                                                       plot.tag = element_text(size = rel(1.5)),
                                                       plot.margin=unit(c(0.1,0.1,0.5,0.1), 'cm'))+ labs(tag = "A")

# Supplementary Figure 6
grid.arrange(p3, p1, layout_matrix=rbind(c(1), c(2)))


##### Compare author-reported burdens with our burdens #####

# Mutational burden - Supplementary figure 17
x <- immunorx$Total_mutations
y <- immunorx$Original_total_mutations
summary(rlm(y~x)) # Get slope/intercept
par(mar=c(4.5,5,0.75,1))
plot(x, y, ylim=c(0, 26115), xlim=c(0,73160), 
     xlab="Total mutation burden (raw #)", ylab="Author-reported mutation burden (raw #)", 
     cex.lab=2, cex.axis=1.75, pch=19, cex=0.5)
abline(rlm(y~x), col='red')
text(65500, -500, labels="y = 0.04x + 473.8",cex=1.5, col="black") # adjust for slope/intercept

cor.test(x,y)

# Nonsynonymous mutational burden - supplementary figure 20
x <- immunorx$Nonsynonymous_SNVs
y <- immunorx$Original_nonsynonymous_mutations
summary(rlm(y~x)) # Get slope/intercept
par(mar=c(4.5,5,0.75,1))
plot(x, y, 
     xlab="Nonsynonymous mutation burden (raw #)", ylab="Author-reported nonsynonymous mutation burden (raw #)", 
     cex.lab=1.7, cex.axis=1.25, pch=19, cex=0.5)
abline(rlm(y~x), col='red')
text(12500, -100, labels="y = 0.90x - 3.40",cex=1.5, col="black") # adjust for slope/intercept

cor.test(x,y)

# Neoantigen burden - supplementary figure 21
x <- immunorx$Manuscript_binding_eps
y <- immunorx$Original_neoantigens
summary(rlm(y~x)) # Get slope/intercept
par(mar=c(4.5,5,0.75,1))
plot(x, y, 
     xlab="Neoepitope burden (# MHC-binding)", ylab="Author-reported neoepitope burden (# MHC-binding)", 
     cex.lab=1.7, cex.axis=1.25, pch=19, cex=0.5)
abline(rlm(y~x), col='red')
text(150000, 16500, labels="y = -0.0015x + 377.60",cex=1.5, col="black") # adjust for slope/intercept

cor.test(x,y)

# Disagreement in rank across all burdens - supplementary figure 19
par(mar=c(4.5,4.75,0.5,0.5))
plot(0,0,xlim=c(0,100),ylim=c(0,50),type="n",axes=T,xlab="Burden 'high' threshold (%ile)",
     ylab="Disagreement (% of patients)", cex.lab=2, cex.axis=1.5)
box(lty="solid", col="black")
# plot nonsynonymous mutations
temp <- function(x) {
  length(which(immunorx$Nonsynonymous_SNVs>=quantile(immunorx$Nonsynonymous_SNVs,x,na.rm=TRUE) & immunorx$Original_nonsynonymous_mutations < quantile(immunorx$Original_nonsynonymous_mutations,x,na.rm=TRUE))) + length(which(immunorx$Nonsynonymous_SNVs<quantile(immunorx$Nonsynonymous_SNVs,x,na.rm=TRUE) & immunorx$Original_nonsynonymous_mutations >=quantile(immunorx$Original_nonsynonymous_mutations,x,na.rm=TRUE)))
}
data <- cbind(1:99, 0)
for (i in 1:99) {
  data[i,2] <- 100*temp(i/100)/length(which(!is.na(immunorx$Original_nonsynonymous_mutations)))
}
points(data[,1], data[,2], type="l", col='gray')
# plot all mutations
temp <- function(x) {
  length(which(immunorx$Total_mutations>=quantile(immunorx$Total_mutations,x,na.rm=TRUE) & immunorx$Original_total_mutations < quantile(immunorx$Original_total_mutations,x,na.rm=TRUE))) + length(which(immunorx$Total_mutations<quantile(immunorx$Total_mutations,x,na.rm=TRUE) & immunorx$Original_total_mutations >=quantile(immunorx$Original_total_mutations,x,na.rm=TRUE)))
}
data <- cbind(1:99, 0)
for (i in 1:99) {
  data[i,2] <- 100*temp(i/100)/length(which(!is.na(immunorx$Original_total_mutations)))
}
points(data[,1], data[,2], type="l", col='black')
# plot epitopes
temp <- function(x) {
  length(which(immunorx$Manuscript_binding_eps>=quantile(immunorx$Manuscript_binding_eps,x,na.rm=TRUE) & immunorx$Original_neoantigens < quantile(immunorx$Original_neoantigens,x,na.rm=TRUE))) + length(which(immunorx$Manuscript_binding_eps<quantile(immunorx$Manuscript_binding_eps,x,na.rm=TRUE) & immunorx$Original_neoantigens >=quantile(immunorx$Original_neoantigens,x,na.rm=TRUE)))
}
data <- cbind(1:99, 0)
for (i in 1:99) {
  data[i,2] <- 100*temp(i/100)/length(which(!is.na(immunorx$Original_neoantigens)))
}
points(data[,1], data[,2], type="l", col='blue')
rect(xleft=66, xright=103, ytop=50, ybottom=42, col="white", border="black")
rect(xleft=67, xright=72, ytop=49.5, ybottom=47.5, col="black", border="black")
rect(xleft=67, xright=72, ytop=47, ybottom=45, col="gray", border="black")
rect(xleft=67, xright=72, ytop=44.5, ybottom=42.5, col="blue", border="black")
text(81.5, 48.5, 'Total mutations')
text(88, 46, 'Nonsynonymous mutations')
text(80.5, 43.5, 'Neoantigens')


##### Extended neoepitope burden ##### 

extended <- read.table('summarized_epitope_burden.tsv', sep='\t', header=T)
colnames(extended)[2] <- "Tumor_ID"
merged_extended <- merge(immunorx, extended, by=c("Patient", "Tumor_ID"))


##### Incorporate neojunction burden #####

jx_data <- read.table('patient_jx_burdens.tsv', sep='\t', header=T)
merged_jx <- merge(immunorx, jx_data, by=c("Patient", "Tumor_ID"))

median(jx_data$Jx_burden)
median(merged_jx$Jx_burden[merged_jx$Disease == 'melanoma'])
median(merged_jx$Jx_burden[merged_jx$Disease == 'NSCLC'])
median(merged_jx$Jx_burden[merged_jx$Disease == 'prostate'])
median(merged_jx$Jx_burden[merged_jx$Disease == 'RCC'])
median(merged_jx$Jx_burden[merged_jx$Disease %in% c("colon", "endometrial", "thyroid")])


##### Incorporate retained introns #####

ri <- read.table('full_intron_retention_burden.tsv', sep='\t', header=T)

# Merge data
full_merged_rna <- merge(merged_extended, jx_data, by=c("Patient", "Tumor_ID"))
fuller_merged_rna <- merge(full_merged_rna, ri, by=c("Patient", "Tumor_ID"))
merged_rna <- fuller_merged_rna[!is.na(fuller_merged_rna$Jx_burden)]

median(merged_rna$Intron_burden)
median(merged_rna$Intron_burden[merged_rna$Intron_burden > 0])
median(merged_rna$Binding_intron_epitope_burden)
median(merged_rna$Binding_intron_epitope_burden[merged_rna$Binding_intron_epitope_burden > 0])

complete_merged_rna <- merge(merged_extended, jx_data, by=c("Patient", "Tumor_ID"), all=T)
final_complete_merged_rna <- merge(complete_merged_rna, ri, by=c("Patient", "Tumor_ID"), all=T)
write.table(final_complete_merged_rna, file="merged_rna_data.tsv", quote=F, col.names = T, row.names = F)

##### Investigating TVB #####

merged_rna$TVB <- merged_rna$Total_mutations + merged_rna$Jx_burden + merged_rna$Intron_burden
merged_rna$Delta_burden <- merged_rna$TVB - merged_rna$Total_mutations
merged_rna$RNA_variants <- merged_rna$Jx_burden + merged_rna$Intron_burden

# Proportion of RNA variation
mean(merged_rna$RNA_variants/merged_rna$TVB)
max(merged_rna$RNA_variants/merged_rna$TVB)

median(merged_rna$Total_mutations)
median(merged_rna$TVB)
median(merged_rna$Delta_burden)

# Load mutations for bar plotting
rna_enumerated_mutations <- read.table('rna_enumerated_mutations.tsv', sep='\t', header=T)
rna_enumerated_mutations$Class <- factor(rna_enumerated_mutations$Class, levels = rev(levels(rna_enumerated_mutations$Class)))

# Figure 2
ggplot(rna_enumerated_mutations, 
       aes(reorder(Patient, 
                   -Mut_count)))+geom_bar(aes(fill=Class))+scale_fill_brewer(palette="Set1")+labs(title="",
                                                                                                   x="Patients", 
                                                                                                   y = "Tumor Variant Burden (#)")+theme(axis.text.x = element_blank(),
                                                                                                                                         panel.grid  = element_blank(),
                                                                                                                                         #axis.title.x = element_blank(),
                                                                                                                                         axis.text = element_text(size = rel(1.5)),
                                                                                                                                         axis.title.y = element_text(size = rel(1.5),
                                                                                                                                                                     margin = margin(t = 0, r = 15, b = 0, l = 0)),
                                                                                                                                         axis.title.x = element_text(size = rel(1.5)),
                                                                                                                                         panel.background = element_rect(fill="white", colour = 'white'),
                                                                                                                                         axis.line.y.left = element_line(colour = 'black'),
                                                                                                                                         axis.line.x.bottom = element_line(colour = 'black'),
                                                                                                                                         axis.ticks.x = element_blank(),
                                                                                                                                         legend.text = element_text(size = rel(1.5)),
                                                                                                                                         legend.title = element_text(size = rel(1.5)),
                                                                                                                                         legend.position = c(0.95, 0.5),
                                                                                                                                         plot.margin=unit(c(0,0.25,0,0.5), 'cm'))+scale_y_continuous(expand = c(0,0))+scale_x_discrete(expand = c(0.01,0))



# Load neoepitopes for plotting
rna_enumerated_neoepitopes <- read.table('rna_enumerated_neoepitopes.tsv', sep='\t', header=T)
rna_enumerated_neoepitopes$Class <- factor(rna_enumerated_neoepitopes$Class, levels = rev(levels(rna_enumerated_neoepitopes$Class)))

# Supplementary Figure 5                                                                                                                         plot.margin=unit(c(0,0.25,0,0.5), 'cm'))
color_palette <- brewer.pal(3, "Dark2")
ggplot(rna_enumerated_neoepitopes, 
       aes(reorder(Patient, 
                   -Mut_count)))+geom_bar(aes(fill=Class))+scale_fill_manual(values = c(color_palette[1], color_palette[3]))+labs(title="",
                                                                                                   x="Patients", 
                                                                                                   y = "Tumor Neoepitope Burden (#)")+theme(axis.text.x = element_blank(),
                                                                                                                                         panel.grid  = element_blank(),
                                                                                                                                         #axis.title.x = element_blank(),
                                                                                                                                         axis.text = element_text(size = rel(1.5)),
                                                                                                                                         axis.title.y = element_text(size = rel(1.5),
                                                                                                                                                                     margin = margin(t = 0, r = 15, b = 0, l = 0)),
                                                                                                                                         axis.title.x = element_text(size = rel(1.5)),
                                                                                                                                         panel.background = element_rect(fill="white", colour = 'white'),
                                                                                                                                         axis.line.y.left = element_line(colour = 'black'),
                                                                                                                                         axis.line.x.bottom = element_line(colour = 'black'),
                                                                                                                                         axis.ticks.x = element_blank(),
                                                                                                                                         legend.text = element_text(size = rel(1.5)),
                                                                                                                                         legend.title = element_text(size = rel(1.5)),
                                                                                                                                         legend.position = c(0.95, 0.5),
                                                                                                                                         plot.margin=unit(c(0,0.25,0,0.5), 'cm'))+scale_y_continuous(expand = c(0,0))+scale_x_discrete(expand = c(0.01,0))
##### NetCTLpan epitopes #####

netctlpan <- read.table('netctlpan_burdens.tsv', sep='\t', header=T, row.names = NULL)

immunorx_netctl <- merge(immunorx, netctlpan, by=c("Patient", "Tumor_ID"))
melanoma_ctl <- immunorx_netctl[immunorx_netctl$Disease == 'melanoma',]
RCC_ctl <- immunorx_netctl[immunorx_netctl$Disease == 'RCC',]
NSCLC_ctl <- immunorx_netctl[immunorx_netctl$Disease == 'NSCLC',]
MMR_ctl <- subset(immunorx_netctl, immunorx_netctl$Disease == 'colon' | immunorx_netctl$Disease == 'endometrial' | immunorx_netctl$Disease == 'thyroid')
prostate_ctl <- immunorx_netctl[immunorx_netctl$Disease == 'prostate',]

median(immunorx_netctl$NetCTLpan_epitopes)
median(melanoma_ctl$NetCTLpan_epitopes)
median(RCC_ctl$NetCTLpan_epitopes)
median(NSCLC_ctl$NetCTLpan_epitopes)
median(MMR_ctl$NetCTLpan_epitopes)
median(prostate_ctl$NetCTLpan_epitopes)


##### Predicting response to immunotherapy #####

# Set function for finding cells in tables
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)}


### MHCnuggets 1 vs 2 - supplementary figure 9 ####

roc.all.mhcnuggets1 <- roc(Combined_response~MHCnuggets_ClassI_eps, data=immunorx)
roc.all.mhcnuggets2 <- roc(Combined_response~MHCnuggets_ClassI_eps, data=immunorx)

roc.mel.mhcnuggets1 <- roc(Combined_response~MHCnuggets_ClassI_eps, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.mhcnuggets2 <- roc(Combined_response~MHCnuggets_ClassI_eps, data=immunorx[immunorx$Disease == 'melanoma',])

roc.rcc.mhcnuggets1 <- roc(Combined_response~MHCnuggets_ClassI_eps, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.mhcnuggets2 <- roc(Combined_response~MHCnuggets_ClassI_eps, data=immunorx[immunorx$Disease == 'RCC',])

roc.nsclc.mhcnuggets1 <- roc(Combined_response~MHCnuggets_ClassI_eps, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.mhcnuggets2 <- roc(Combined_response~MHCnuggets_ClassI_eps, data=immunorx[immunorx$Disease == 'NSCLC',])

roc_list_all <- list(roc.all.mhcnuggets1, roc.all.mhcnuggets2)
mel_cancer_roc_list <- list(roc.mel.mhcnuggets1, roc.mel.mhcnuggets2)
rcc_cancer_roc_list <- list(roc.rcc.mhcnuggets1, roc.rcc.mhcnuggets2)
nsclc_cancer_roc_list <- list(roc.nsclc.mhcnuggets1, roc.nsclc.mhcnuggets2)
color_set <- c(brewer.pal(2, "Dark2")) 
label_set <- c("MHCnuggets Class I", "MHCnuggets Class II")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.nsclc <- ggroc(nsclc_cancer_roc_list,
                  linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("NSCLC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC', 'NSCLC'),
                             N=c(431, 302, 57, 34), 
                             MHCnuggets_1=c(sprintf("%.3f", round(roc.all.mhcnuggets1$auc, 3)),
                                            sprintf("%.3f", round(roc.mel.mhcnuggets1$auc, 3)),
                                            sprintf("%.3f", round(roc.rcc.mhcnuggets1$auc, 3)),
                                            sprintf("%.3f", round(roc.nsclc.mhcnuggets1$auc, 3))),
                             MHCnuggets_2=c(sprintf("%.3f", round(roc.all.mhcnuggets2$auc, 3)),
                                            sprintf("%.3f", round(roc.mel.mhcnuggets2$auc, 3)),
                                            sprintf("%.3f", round(roc.rcc.mhcnuggets2$auc, 3)),
                                            sprintf("%.3f", round(roc.nsclc.mhcnuggets2$auc, 3))),
                             stringsAsFactors = F)
colnames(table_contents) <- c("Cancer Type", "N", 
                              "MHC Class I Epitopes", "MHC Class II Epitopes")

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", "gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", "gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# list of cells to bold (r,c) 
cell_index <- list(c(3,3), c(4,3), c(5,3), c(6,5))

for (cell in cell_index){
  summary_tab$grobs[find_cell(summary_tab, cell[1], cell[2])][[1]][["gp"]] <- gpar(fontface="bold")
}

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"),
                         gg.nsclc+theme(legend.position = "none"),ncol= 2),
             textGrob("AUC by Neoepitope Burden and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")))


### Author-reported burden values - Supplementary Figure 18 ###

# Create models
roc.all.total <- roc(Combined_response~Original_total_mutations, data=immunorx)
roc.all.ns <- roc(Combined_response~Original_nonsynonymous_mutations, data=immunorx)
roc.all.eps <- roc(Combined_response~Original_neoantigens, data=immunorx)

roc.mel.total <- roc(Combined_response~Original_total_mutations, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.ns <- roc(Combined_response~Original_nonsynonymous_mutations, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.eps <- roc(Combined_response~Original_neoantigens, data=immunorx[immunorx$Disease == 'melanoma',])

roc.rcc.total <- roc(Combined_response~Original_total_mutations, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.ns <- roc(Combined_response~Original_nonsynonymous_mutations, data=immunorx[immunorx$Disease == 'RCC',])

roc.nsclc.total <- roc(Combined_response~Original_total_mutations, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.ns <- roc(Combined_response~Original_nonsynonymous_mutations, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.eps <- roc(Combined_response~Original_neoantigens, data=immunorx[immunorx$Disease == 'NSCLC',])

roc_list_all <- list(roc.all.total, roc.all.ns, roc.all.eps)
mel_cancer_roc_list <- list(roc.mel.total, roc.mel.ns, roc.mel.eps)
rcc_cancer_roc_list <- list(roc.rcc.total, roc.rcc.ns)
nsclc_cancer_roc_list <- list(roc.nsclc.total, roc.nsclc.ns, roc.nsclc.eps)
color_set <- c(brewer.pal(3, "Accent")) 
label_set <- c("Total mutations", "Nonsynonymous mutations", "Neoantigens")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set[1:2],
                     values = color_set[1:2])+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.nsclc <- ggroc(nsclc_cancer_roc_list,
                  linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("NSCLC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC', 'NSCLC'),
                             #N=c(409, 285, 56, 34), 
                             Total_mutations=c(sprintf("%.3f", round(roc.all.total$auc, 3)),
                                               sprintf("%.3f", round(roc.mel.total$auc, 3)),
                                               sprintf("%.3f", round(roc.rcc.total$auc, 3)),
                                               sprintf("%.3f", round(roc.nsclc.total$auc, 3))),
                             Nonsynonymous_mutations=c(sprintf("%.3f", round(roc.all.ns$auc, 3)),
                                                       sprintf("%.3f", round(roc.mel.ns$auc, 3)),
                                                       sprintf("%.3f", round(roc.rcc.ns$auc, 3)),
                                                       sprintf("%.3f", round(roc.nsclc.ns$auc, 3))),
                             Neoantigens=c(sprintf("%.3f", round(roc.all.eps$auc, 3)),
                                           sprintf("%.3f", round(roc.mel.eps$auc, 3)),
                                           "NA",
                                           sprintf("%.3f", round(roc.nsclc.eps$auc, 3))),
                             stringsAsFactors = F)
colnames(table_contents) <- c("Cancer type", "Total mutations", "Nonsynonymous mutations", "Neoantigens")

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# list of cells to bold (r,c) 
cell_index <- list(c(3,3), c(4,2), c(5,2), c(6,3))

for (cell in cell_index){
  summary_tab$grobs[find_cell(summary_tab, cell[1], cell[2])][[1]][["gp"]] <- gpar(fontface="bold")
}

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"),
                         gg.nsclc+theme(legend.position = "none"),ncol= 2),
             textGrob("AUC by Author-reported Burden and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")))


### Genomic coverage - Supplementary Figure 11 ###

# Create models
roc.all.coverage <- roc(Combined_response~Coverage, data=immunorx)
roc.mel.coverage <- roc(Combined_response~Coverage, data=immunorx[immunorx$Disease == 'melanoma',])
roc.rcc.coverage <- roc(Combined_response~Coverage, data=immunorx[immunorx$Disease == 'RCC',])
roc.nsclc.coverage <- roc(Combined_response~Coverage, data=immunorx[immunorx$Disease == 'NSCLC',])

roc_list_all <- list(roc.all.coverage)
mel_cancer_roc_list <- list(roc.mel.coverage)
rcc_cancer_roc_list <- list(roc.rcc.coverage)
nsclc_cancer_roc_list <- list(roc.nsclc.coverage)
color_set <- c("grey50")
label_set <- c("Coverage (Mbp)")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.nsclc <- ggroc(nsclc_cancer_roc_list,
                  linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("NSCLC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC', 'NSCLC'),
                             N=c(431, 302, 57, 34), 
                             Coverage=c(sprintf("%.3f", round(roc.all.coverage$auc, 3)),
                                         sprintf("%.3f", round(roc.mel.coverage$auc, 3)),
                                         sprintf("%.3f", round(roc.rcc.coverage$auc, 3)),
                                         sprintf("%.3f", round(roc.nsclc.coverage$auc, 3))),
                             stringsAsFactors = F)
colnames(table_contents)[1] <- "Cancer type" 

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", "gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", "gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"),
                         gg.nsclc+theme(legend.position = "none"),ncol= 2),
             textGrob("AUC by Coverage and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")))


### Raw mutational burden by different callers - Supplementary Figure 18 ###

# Create models
roc.all.consensus <- roc(Combined_response~Total_mutations, data=immunorx)
roc.all.muse <- roc(Combined_response~Muse_variants, data=immunorx)
roc.all.mutect <- roc(Combined_response~Mutect_variants, data=immunorx)
roc.all.pindel <- roc(Combined_response~Pindel_variants, data=immunorx)
roc.all.radia <- roc(Combined_response~Radia_variants, data=immunorx)
roc.all.somaticsniper <- roc(Combined_response~Somaticsniper_variants, data=immunorx)
roc.all.varscan <- roc(Combined_response~Varscan_variants, data=immunorx)

roc.mel.consensus <- roc(Combined_response~Total_mutations, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.muse <- roc(Combined_response~Muse_variants, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.mutect <- roc(Combined_response~Mutect_variants, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.pindel <- roc(Combined_response~Pindel_variants, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.radia <- roc(Combined_response~Radia_variants, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.somaticsniper <- roc(Combined_response~Somaticsniper_variants, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.varscan <- roc(Combined_response~Varscan_variants, data=immunorx[immunorx$Disease == 'melanoma',])

roc.rcc.consensus <- roc(Combined_response~Total_mutations, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.muse <- roc(Combined_response~Muse_variants, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.mutect <- roc(Combined_response~Mutect_variants, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.pindel <- roc(Combined_response~Pindel_variants, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.radia <- roc(Combined_response~Radia_variants, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.somaticsniper <- roc(Combined_response~Somaticsniper_variants, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.varscan <- roc(Combined_response~Varscan_variants, data=immunorx[immunorx$Disease == 'RCC',])

roc.nsclc.consensus <- roc(Combined_response~Total_mutations, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.muse <- roc(Combined_response~Muse_variants, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.mutect <- roc(Combined_response~Mutect_variants, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.pindel <- roc(Combined_response~Pindel_variants, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.radia <- roc(Combined_response~Radia_variants, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.somaticsniper <- roc(Combined_response~Somaticsniper_variants, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.varscan <- roc(Combined_response~Varscan_variants, data=immunorx[immunorx$Disease == 'NSCLC',])

roc_list_all <- list(roc.all.consensus, roc.all.muse, roc.all.mutect, roc.all.pindel,
                     roc.all.radia, roc.all.somaticsniper, roc.all.varscan)
mel_cancer_roc_list <- list(roc.mel.consensus, roc.mel.muse, roc.mel.mutect, roc.mel.pindel,
                            roc.mel.radia, roc.mel.somaticsniper, roc.mel.varscan)
rcc_cancer_roc_list <- list(roc.rcc.consensus, roc.rcc.muse, roc.rcc.mutect, roc.rcc.pindel,
                            roc.rcc.radia, roc.rcc.somaticsniper, roc.rcc.varscan)
nsclc_cancer_roc_list <- list(roc.nsclc.consensus, roc.nsclc.muse, roc.nsclc.mutect, roc.nsclc.pindel,
                              roc.nsclc.radia, roc.nsclc.somaticsniper, roc.nsclc.varscan)
color_set <- c("grey50", brewer.pal(6, "Set2")) 
label_set <- c("Consensus", "MuSE", "MuTect", "Pindel", "RADIA", "SomaticSniper", "VarScan")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.nsclc <- ggroc(nsclc_cancer_roc_list,
                  linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("NSCLC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC', 'NSCLC'),
                             N=c(431, 302, 57, 34), 
                             Consensus=c(sprintf("%.3f", round(roc.all.consensus$auc, 3)),
                                         sprintf("%.3f", round(roc.mel.consensus$auc, 3)),
                                         sprintf("%.3f", round(roc.rcc.consensus$auc, 3)),
                                         sprintf("%.3f", round(roc.nsclc.consensus$auc, 3))),
                             MuSE=c(sprintf("%.3f", round(roc.all.muse$auc, 3)),
                                    sprintf("%.3f", round(roc.mel.muse$auc, 3)),
                                    sprintf("%.3f", round(roc.rcc.muse$auc, 3)),
                                    sprintf("%.3f", round(roc.nsclc.muse$auc, 3))),
                             MuTect=c(sprintf("%.3f", round(roc.all.mutect$auc, 3)),
                                      sprintf("%.3f", round(roc.mel.mutect$auc, 3)),
                                      sprintf("%.3f", round(roc.rcc.mutect$auc, 3)),
                                      sprintf("%.3f", round(roc.nsclc.mutect$auc, 3))),
                             Pindel=c(sprintf("%.3f", round(roc.all.pindel$auc, 3)),
                                      sprintf("%.3f", round(roc.mel.pindel$auc, 3)),
                                      sprintf("%.3f", round(roc.rcc.pindel$auc, 3)),
                                      sprintf("%.3f", round(roc.nsclc.pindel$auc, 3))),
                             RADIA=c(sprintf("%.3f", round(roc.all.radia$auc, 3)),
                                     sprintf("%.3f", round(roc.mel.radia$auc, 3)),
                                     sprintf("%.3f", round(roc.rcc.radia$auc, 3)),
                                     sprintf("%.3f", round(roc.nsclc.radia$auc, 3))),
                             SomaticSniper=c(sprintf("%.3f", round(roc.all.somaticsniper$auc, 3)),
                                             sprintf("%.3f", round(roc.mel.somaticsniper$auc, 3)),
                                             sprintf("%.3f", round(roc.rcc.somaticsniper$auc, 3)),
                                             sprintf("%.3f", round(roc.nsclc.somaticsniper$auc, 3))),
                             VarScan=c(sprintf("%.3f", round(roc.all.varscan$auc, 3)),
                                       sprintf("%.3f", round(roc.mel.varscan$auc, 3)),
                                       sprintf("%.3f", round(roc.rcc.varscan$auc, 3)),
                                       sprintf("%.3f", round(roc.nsclc.varscan$auc, 3))), stringsAsFactors = F)
colnames(table_contents)[1] <- "Cancer type" 

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", "gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", "gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# list of cells to bold (r,c) 
cell_index <- list(c(3,9), c(4,5), c(5,8), c(6,4))

for (cell in cell_index){
  summary_tab$grobs[find_cell(summary_tab, cell[1], cell[2])][[1]][["gp"]] <- gpar(fontface="bold")
}

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"),
                         gg.nsclc+theme(legend.position = "none"),ncol= 2),
             textGrob("AUC by Mutation Burden and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")))


## Coverage-adjusted mutational burden by different callers - Figure 7 ##

# Create models
roc.all.consensus <- roc(Combined_response~Total_mutations_adj, data=immunorx)
roc.all.muse <- roc(Combined_response~Muse_variants_adj, data=immunorx)
roc.all.mutect <- roc(Combined_response~Mutect_variants_adj, data=immunorx)
roc.all.pindel <- roc(Combined_response~Pindel_variants_adj, data=immunorx)
roc.all.radia <- roc(Combined_response~Radia_variants_adj, data=immunorx)
roc.all.somaticsniper <- roc(Combined_response~Somaticsniper_variants_adj, data=immunorx)
roc.all.varscan <- roc(Combined_response~Varscan_variants_adj, data=immunorx)

roc.mel.consensus <- roc(Combined_response~Total_mutations_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.muse <- roc(Combined_response~Muse_variants_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.mutect <- roc(Combined_response~Mutect_variants_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.pindel <- roc(Combined_response~Pindel_variants_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.radia <- roc(Combined_response~Radia_variants_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.somaticsniper <- roc(Combined_response~Somaticsniper_variants_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.varscan <- roc(Combined_response~Varscan_variants_adj, data=immunorx[immunorx$Disease == 'melanoma',])

roc.rcc.consensus <- roc(Combined_response~Total_mutations_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.muse <- roc(Combined_response~Muse_variants_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.mutect <- roc(Combined_response~Mutect_variants_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.pindel <- roc(Combined_response~Pindel_variants_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.radia <- roc(Combined_response~Radia_variants_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.somaticsniper <- roc(Combined_response~Somaticsniper_variants_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.varscan <- roc(Combined_response~Varscan_variants_adj, data=immunorx[immunorx$Disease == 'RCC',])

roc.nsclc.consensus <- roc(Combined_response~Total_mutations_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.muse <- roc(Combined_response~Muse_variants_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.mutect <- roc(Combined_response~Mutect_variants_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.pindel <- roc(Combined_response~Pindel_variants_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.radia <- roc(Combined_response~Radia_variants_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.somaticsniper <- roc(Combined_response~Somaticsniper_variants_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.varscan <- roc(Combined_response~Varscan_variants_adj, data=immunorx[immunorx$Disease == 'NSCLC',])

roc_list_all <- list(roc.all.consensus, roc.all.muse, roc.all.mutect, roc.all.pindel,
                     roc.all.radia, roc.all.somaticsniper, roc.all.varscan)
mel_cancer_roc_list <- list(roc.mel.consensus, roc.mel.muse, roc.mel.mutect, roc.mel.pindel,
                            roc.mel.radia, roc.mel.somaticsniper, roc.mel.varscan)
rcc_cancer_roc_list <- list(roc.rcc.consensus, roc.rcc.muse, roc.rcc.mutect, roc.rcc.pindel,
                            roc.rcc.radia, roc.rcc.somaticsniper, roc.rcc.varscan)
nsclc_cancer_roc_list <- list(roc.nsclc.consensus, roc.nsclc.muse, roc.nsclc.mutect, roc.nsclc.pindel,
                              roc.nsclc.radia, roc.nsclc.somaticsniper, roc.nsclc.varscan)
color_set <- c("grey50", brewer.pal(6, "Set2")) # change num. to N-1, set2 max is 8, set3 is 12
# color_set <- brewer.pal(7, "Set2")
label_set <- c("Consensus", "MuSE", "MuTect", "Pindel", "RADIA", "SomaticSniper", "VarScan")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.nsclc <- ggroc(nsclc_cancer_roc_list,
                  linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("NSCLC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC', 'NSCLC'),
                             N=c(431, 302, 57, 34), 
                             Consensus=c(sprintf("%.3f", round(roc.all.consensus$auc, 3)),
                                         sprintf("%.3f", round(roc.mel.consensus$auc, 3)),
                                         sprintf("%.3f", round(roc.rcc.consensus$auc, 3)),
                                         sprintf("%.3f", round(roc.nsclc.consensus$auc, 3))),
                             MuSE=c(sprintf("%.3f", round(roc.all.muse$auc, 3)),
                                    sprintf("%.3f", round(roc.mel.muse$auc, 3)),
                                    sprintf("%.3f", round(roc.rcc.muse$auc, 3)),
                                    sprintf("%.3f", round(roc.nsclc.muse$auc, 3))),
                             MuTect=c(sprintf("%.3f", round(roc.all.mutect$auc, 3)),
                                      sprintf("%.3f", round(roc.mel.mutect$auc, 3)),
                                      sprintf("%.3f", round(roc.rcc.mutect$auc, 3)),
                                      sprintf("%.3f", round(roc.nsclc.mutect$auc, 3))),
                             Pindel=c(sprintf("%.3f", round(roc.all.pindel$auc, 3)),
                                      sprintf("%.3f", round(roc.mel.pindel$auc, 3)),
                                      sprintf("%.3f", round(roc.rcc.pindel$auc, 3)),
                                      sprintf("%.3f", round(roc.nsclc.pindel$auc, 3))),
                             RADIA=c(sprintf("%.3f", round(roc.all.radia$auc, 3)),
                                     sprintf("%.3f", round(roc.mel.radia$auc, 3)),
                                     sprintf("%.3f", round(roc.rcc.radia$auc, 3)),
                                     sprintf("%.3f", round(roc.nsclc.radia$auc, 3))),
                             SomaticSniper=c(sprintf("%.3f", round(roc.all.somaticsniper$auc, 3)),
                                             sprintf("%.3f", round(roc.mel.somaticsniper$auc, 3)),
                                             sprintf("%.3f", round(roc.rcc.somaticsniper$auc, 3)),
                                             sprintf("%.3f", round(roc.nsclc.somaticsniper$auc, 3))),
                             VarScan=c(sprintf("%.3f", round(roc.all.varscan$auc, 3)),
                                       sprintf("%.3f", round(roc.mel.varscan$auc, 3)),
                                       sprintf("%.3f", round(roc.rcc.varscan$auc, 3)),
                                       sprintf("%.3f", round(roc.nsclc.varscan$auc, 3))), stringsAsFactors = F)
colnames(table_contents)[1] <- "Cancer type" 

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", "gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", "gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# list of cells to bold (r,c) 
cell_index <- list(c(3,9), c(4,5), c(5,8), c(6,4))

for (cell in cell_index){
  summary_tab$grobs[find_cell(summary_tab, cell[1], cell[2])][[1]][["gp"]] <- gpar(fontface="bold")
}

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"),
                         gg.nsclc+theme(legend.position = "none"), ncol= 2),
             textGrob("AUC by Coverage-adjusted Mutation Burden and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")))


## Coverage-adjusted mutational burden by variant type - Supplementary Figure 8 ##

# Create models
roc.all.total <- roc(Combined_response~Total_mutations_adj, data=immunorx)
roc.all.snv <- roc(Combined_response~SNVs_adj, data=immunorx)
roc.all.indel <- roc(Combined_response~Total_indels_adj, data=immunorx)
roc.all.fs <- roc(Combined_response~Frameshift_indels_adj, data=immunorx)
roc.all.if <- roc(Combined_response~Inframe_indels_adj, data=immunorx)

roc.mel.total <- roc(Combined_response~Total_mutations_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.snv <- roc(Combined_response~SNVs_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.indel <- roc(Combined_response~Total_indels_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.fs <- roc(Combined_response~Frameshift_indels_adj, data=immunorx[immunorx$Disease == 'melanoma',])
roc.mel.if <- roc(Combined_response~Inframe_indels_adj, data=immunorx[immunorx$Disease == 'melanoma',])

roc.rcc.total <- roc(Combined_response~Total_mutations_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.snv <- roc(Combined_response~SNVs_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.indel <- roc(Combined_response~Total_indels_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.fs <- roc(Combined_response~Frameshift_indels_adj, data=immunorx[immunorx$Disease == 'RCC',])
roc.rcc.if <- roc(Combined_response~Inframe_indels_adj, data=immunorx[immunorx$Disease == 'RCC',])

roc.nsclc.total <- roc(Combined_response~Total_mutations_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.snv <- roc(Combined_response~SNVs_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.indel <- roc(Combined_response~Total_indels_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.fs <- roc(Combined_response~Frameshift_indels_adj, data=immunorx[immunorx$Disease == 'NSCLC',])
roc.nsclc.if <- roc(Combined_response~Inframe_indels_adj, data=immunorx[immunorx$Disease == 'NSCLC',])

roc_list_all <- list(roc.all.total, roc.all.snv, roc.all.indel, roc.all.fs, roc.all.if)
mel_cancer_roc_list <- list(roc.mel.total, roc.mel.snv, roc.mel.indel, roc.mel.fs, roc.mel.if)
rcc_cancer_roc_list <- list(roc.rcc.total, roc.rcc.snv, roc.rcc.indel, roc.rcc.fs, roc.rcc.if)
nsclc_cancer_roc_list <- list(roc.nsclc.total, roc.nsclc.snv, roc.nsclc.indel, roc.nsclc.fs, roc.nsclc.if)
color_set <- c("grey50", brewer.pal(5, "Set2")) # change num. to N-1, set2 max is 8, set3 is 12
label_set <- c("All", "SNVs", "Indels", "FS indels", "In-frame indels")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.nsclc <- ggroc(nsclc_cancer_roc_list,
                  linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Variant Caller",
                     labels = label_set,
                     values = color_set)+
  ggtitle("NSCLC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC', 'NSCLC'),
                             N=c(431, 302, 57, 34), 
                             All=c(sprintf("%.3f", round(roc.all.total$auc, 3)),
                                         sprintf("%.3f", round(roc.mel.total$auc, 3)),
                                         sprintf("%.3f", round(roc.rcc.total$auc, 3)),
                                         sprintf("%.3f", round(roc.nsclc.total$auc, 3))),
                             SNVs=c(sprintf("%.3f", round(roc.all.snv$auc, 3)),
                                    sprintf("%.3f", round(roc.mel.snv$auc, 3)),
                                    sprintf("%.3f", round(roc.rcc.snv$auc, 3)),
                                    sprintf("%.3f", round(roc.nsclc.snv$auc, 3))),
                             Indels=c(sprintf("%.3f", round(roc.all.indel$auc, 3)),
                                      sprintf("%.3f", round(roc.mel.indel$auc, 3)),
                                      sprintf("%.3f", round(roc.rcc.indel$auc, 3)),
                                      sprintf("%.3f", round(roc.nsclc.indel$auc, 3))),
                             FS_indels=c(sprintf("%.3f", round(roc.all.fs$auc, 3)),
                                      sprintf("%.3f", round(roc.mel.fs$auc, 3)),
                                      sprintf("%.3f", round(roc.rcc.fs$auc, 3)),
                                      sprintf("%.3f", round(roc.nsclc.fs$auc, 3))),
                             IF_indels=c(sprintf("%.3f", round(roc.all.if$auc, 3)),
                                     sprintf("%.3f", round(roc.mel.if$auc, 3)),
                                     sprintf("%.3f", round(roc.rcc.if$auc, 3)),
                                     sprintf("%.3f", round(roc.nsclc.if$auc, 3))), stringsAsFactors = F)
colnames(table_contents)[1] <- "Cancer type" 
colnames(table_contents)[6] <- "FS indels"
colnames(table_contents)[7] <- "In-frame indels"

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", "gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", "gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# list of cells to bold (r,c) 
cell_index <- list(c(3,4), c(4,4), c(5,7), c(6,4))

for (cell in cell_index){
  summary_tab$grobs[find_cell(summary_tab, cell[1], cell[2])][[1]][["gp"]] <- gpar(fontface="bold")
}

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"),
                         gg.nsclc+theme(legend.position = "none"), ncol= 2),
             textGrob("AUC by Mutation Burden and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")))


### Extended neoepitope burden - Figure 4 ###

merged_extended$Mutations_alleles <- merged_extended$Total_mutations*merged_extended$Tumor_HLA1_count

# Create models
roc.all.consensus <- roc(Combined_response~Total_mutations, data=merged_extended)
roc.all.consensus_adj <- roc(Combined_response~Total_mutations_adj, data=merged_extended)
roc.all.consensus_allele <- roc(Combined_response~Mutations_alleles, data=merged_extended)
roc.all.epitope <- roc(Combined_response~Total_comprehensive_neoepitopes, data=merged_extended)
roc.all.mismatch <- roc(Combined_response~Epitope_by_mismatch_burden, data=merged_extended)
roc.all.allele <- roc(Combined_response~Epitope_by_allele_burden, data=merged_extended)
roc.all.tcga <- roc(Combined_response~Epitope_by_TCGA_burden, data=merged_extended)
roc.all.mismatch_allele <- roc(Combined_response~Epitope_by_mismatch_and_allele_burden, data=merged_extended)
roc.all.mismatch_tcga <- roc(Combined_response~Epitope_by_mismatch_and_TCGA_burden, data=merged_extended)
roc.all.allele_tcga <- roc(Combined_response~Epitope_by_allele_and_TCGA_burden, data=merged_extended)
roc.all.mismatch_alleles_tcga <- roc(Combined_response~Epitope_by_mismatches_alleles_and_TCGA_burden, data=merged_extended)

roc.mel.consensus <- roc(Combined_response~Total_mutations, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.consensus_adj <- roc(Combined_response~Total_mutations_adj, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.consensus_allele <- roc(Combined_response~Mutations_alleles, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.epitope <- roc(Combined_response~Total_comprehensive_neoepitopes, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.mismatch <- roc(Combined_response~Epitope_by_mismatch_burden, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.allele <- roc(Combined_response~Epitope_by_allele_burden, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.tcga <- roc(Combined_response~Epitope_by_TCGA_burden, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.mismatch_allele <- roc(Combined_response~Epitope_by_mismatch_and_allele_burden, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.mismatch_tcga <- roc(Combined_response~Epitope_by_mismatch_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.allele_tcga <- roc(Combined_response~Epitope_by_allele_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'melanoma',])
roc.mel.mismatch_alleles_tcga <- roc(Combined_response~Epitope_by_mismatches_alleles_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'melanoma',])

roc.rcc.consensus <- roc(Combined_response~Total_mutations, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.consensus_adj <- roc(Combined_response~Total_mutations_adj, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.consensus_allele <- roc(Combined_response~Mutations_alleles, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.epitope <- roc(Combined_response~Total_comprehensive_neoepitopes, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.mismatch <- roc(Combined_response~Epitope_by_mismatch_burden, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.allele <- roc(Combined_response~Epitope_by_allele_burden, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.tcga <- roc(Combined_response~Epitope_by_TCGA_burden, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.mismatch_allele <- roc(Combined_response~Epitope_by_mismatch_and_allele_burden, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.mismatch_tcga <- roc(Combined_response~Epitope_by_mismatch_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.allele_tcga <- roc(Combined_response~Epitope_by_allele_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'RCC',])
roc.rcc.mismatch_alleles_tcga <- roc(Combined_response~Epitope_by_mismatches_alleles_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'RCC',])

roc.nsclc.consensus <- roc(Combined_response~Total_mutations, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.consensus_adj <- roc(Combined_response~Total_mutations_adj, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.consensus_allele <- roc(Combined_response~Mutations_alleles, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.epitope <- roc(Combined_response~Total_comprehensive_neoepitopes, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.mismatch <- roc(Combined_response~Epitope_by_mismatch_burden, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.allele <- roc(Combined_response~Epitope_by_allele_burden, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.tcga <- roc(Combined_response~Epitope_by_TCGA_burden, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.mismatch_allele <- roc(Combined_response~Epitope_by_mismatch_and_allele_burden, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.mismatch_tcga <- roc(Combined_response~Epitope_by_mismatch_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.allele_tcga <- roc(Combined_response~Epitope_by_allele_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'NSCLC',])
roc.nsclc.mismatch_alleles_tcga <- roc(Combined_response~Epitope_by_mismatches_alleles_and_TCGA_burden, data=merged_extended[merged_extended$Disease == 'NSCLC',])

roc_list_all <- list(roc.all.consensus, roc.all.consensus_adj, roc.all.consensus_allele, roc.all.epitope, 
                     roc.all.mismatch, roc.all.allele, roc.all.tcga, roc.all.mismatch_allele, roc.all.mismatch_tcga, 
                     roc.all.allele_tcga, roc.all.mismatch_alleles_tcga)
mel_cancer_roc_list <- list(roc.mel.consensus, roc.mel.consensus_adj, roc.mel.consensus_allele, roc.mel.epitope, 
                            roc.mel.mismatch, roc.mel.allele, roc.mel.tcga, roc.mel.mismatch_allele, roc.mel.mismatch_tcga, 
                            roc.mel.allele_tcga, roc.mel.mismatch_alleles_tcga)
rcc_cancer_roc_list <- list(roc.rcc.consensus, roc.rcc.consensus_adj, roc.rcc.consensus_allele, roc.rcc.epitope, 
                            roc.rcc.mismatch, roc.rcc.allele, roc.rcc.tcga, roc.rcc.mismatch_allele, roc.rcc.mismatch_tcga, 
                            roc.rcc.allele_tcga, roc.rcc.mismatch_alleles_tcga)
nsclc_cancer_roc_list <- list(roc.nsclc.consensus, roc.nsclc.consensus_adj, roc.nsclc.consensus_allele, roc.nsclc.epitope, 
                              roc.nsclc.mismatch, roc.nsclc.allele, roc.nsclc.tcga, roc.nsclc.mismatch_allele, roc.nsclc.mismatch_tcga, 
                              roc.nsclc.allele_tcga, roc.nsclc.mismatch_alleles_tcga)
color_set <- c(brewer.pal(11, "Set3"))
label_set <- c("Mutation burden", "Coverage-Adjusted mutation burden", "Mutations by Alleles", "Neoepitope burden", 
               "Amino acid mismatches", "Alleles presenting", "TCGA expression", "Mismatches + Alleles", 
               "Mismatches+TCGA expression","Alleles+TCGA expression", "Mismatches+Alleles+TCGA expression")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Neoepitope burden modification",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Neoepitope burden modification",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Neoepitope burden modification",
                     labels = label_set,
                     values = color_set)+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.nsclc <- ggroc(nsclc_cancer_roc_list,
                  linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "Neoepitope burden modification",
                     labels = label_set,
                     values = color_set)+
  ggtitle("NSCLC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC', 'NSCLC'),
                             N=c(431, 302, 57, 34), 
                             TMB=c(sprintf("%.3f", round(roc.all.consensus$auc, 3)),
                                   sprintf("%.3f", round(roc.mel.consensus$auc, 3)),
                                   sprintf("%.3f", round(roc.rcc.consensus$auc, 3)),
                                   sprintf("%.3f", round(roc.nsclc.consensus$auc, 3))),
                             TMB_adj=c(sprintf("%.3f", round(roc.all.consensus_adj$auc, 3)),
                                       sprintf("%.3f", round(roc.mel.consensus_adj$auc, 3)),
                                       sprintf("%.3f", round(roc.rcc.consensus_adj$auc, 3)),
                                       sprintf("%.3f", round(roc.nsclc.consensus_adj$auc, 3))),
                             TMB_HLA=c(sprintf("%.3f", round(roc.all.consensus_allele$auc, 3)),
                                       sprintf("%.3f", round(roc.mel.consensus_allele$auc, 3)),
                                       sprintf("%.3f", round(roc.rcc.consensus_allele$auc, 3)),
                                       sprintf("%.3f", round(roc.nsclc.consensus_allele$auc, 3))),
                             Epitope=c(sprintf("%.3f", round(roc.all.epitope$auc, 3)),
                                       sprintf("%.3f", round(roc.mel.epitope$auc, 3)),
                                       sprintf("%.3f", round(roc.rcc.epitope$auc, 3)),
                                       sprintf("%.3f", round(roc.nsclc.epitope$auc, 3))),
                             Mismatches=c(sprintf("%.3f", round(roc.all.mismatch$auc, 3)),
                                          sprintf("%.3f", round(roc.mel.mismatch$auc, 3)),
                                          sprintf("%.3f", round(roc.rcc.mismatch$auc, 3)),
                                          sprintf("%.3f", round(roc.nsclc.mismatch$auc, 3))),
                             Alleles=c(sprintf("%.3f", round(roc.all.allele$auc, 3)),
                                       sprintf("%.3f", round(roc.mel.allele$auc, 3)),
                                       sprintf("%.3f", round(roc.rcc.allele$auc, 3)),
                                       sprintf("%.3f", round(roc.nsclc.allele$auc, 3))),
                             TCGA=c(sprintf("%.3f", round(roc.all.tcga$auc, 3)),
                                    sprintf("%.3f", round(roc.mel.tcga$auc, 3)),
                                    sprintf("%.3f", round(roc.rcc.tcga$auc, 3)),
                                    sprintf("%.3f", round(roc.nsclc.tcga$auc, 3))),
                             Mismatch_allele=c(sprintf("%.3f", round(roc.all.mismatch_allele$auc, 3)),
                                               sprintf("%.3f", round(roc.mel.mismatch_allele$auc, 3)),
                                               sprintf("%.3f", round(roc.rcc.mismatch_allele$auc, 3)),
                                               sprintf("%.3f", round(roc.nsclc.mismatch_allele$auc, 3))),
                             Mismatch_TCGA=c(sprintf("%.3f", round(roc.all.mismatch_tcga$auc, 3)),
                                             sprintf("%.3f", round(roc.mel.mismatch_tcga$auc, 3)),
                                             sprintf("%.3f", round(roc.rcc.mismatch_tcga$auc, 3)),
                                             sprintf("%.3f", round(roc.nsclc.mismatch_tcga$auc, 3))), 
                             Allele_TCGA=c(sprintf("%.3f", round(roc.all.allele_tcga$auc, 3)),
                                           sprintf("%.3f", round(roc.mel.allele_tcga$auc, 3)),
                                           sprintf("%.3f", round(roc.rcc.allele_tcga$auc, 3)),
                                           sprintf("%.3f", round(roc.nsclc.allele_tcga$auc, 3))),
                             Mismatch_allele_TCGA=c(sprintf("%.3f", round(roc.all.mismatch_alleles_tcga$auc, 3)),
                                                    sprintf("%.3f", round(roc.mel.mismatch_alleles_tcga$auc, 3)),
                                                    sprintf("%.3f", round(roc.rcc.mismatch_alleles_tcga$auc, 3)),
                                                    sprintf("%.3f", round(roc.nsclc.mismatch_alleles_tcga$auc, 3))), 
                             stringsAsFactors = F)
colnames(table_contents) <-  c("Cancer type", "N", "TMB1", "TMB2", "TMB1*HLA",
                               "NB", "M", "A", "T", 
                               "M*A", "M*T", "A*T", "M*A*T")

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", "gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", "gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# list of cells to bold (r,c) 
cell_index <- list(c(3,3), c(4,5), c(5,7), c(5,10), c(6,6), c(6,8))

for (cell in cell_index){
  summary_tab$grobs[find_cell(summary_tab, cell[1], cell[2])][[1]][["gp"]] <- gpar(fontface="bold")
}

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"),
                         gg.nsclc+theme(legend.position = "none"), ncol= 2),
             textGrob("AUC by Burden and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")), padding=unit(0, "cm"))


##### RNA-seq features - Figure 5 #####

# Create models
roc.all.mut <- roc(Combined_response~Total_mutations, data=merged_rna)
roc.all.tvb <- roc(Combined_response~TVB, data=merged_rna)
roc.all.jx <- roc(Combined_response~Jx_burden, data=merged_rna)
roc.all.ri <- roc(Combined_response~Intron_burden, data=merged_rna)
roc.all.ri_eps <- roc(Combined_response~Binding_intron_epitope_burden, data=merged_rna)
roc.all.expression <- roc(Combined_response~Epitope_by_expression_burden, data=merged_rna)
roc.all.mismatch_expression <- roc(Combined_response~Epitope_by_mismatch_and_expression_burden, data=merged_rna)
roc.all.allele_expression <- roc(Combined_response~Epitope_by_allele_and_expression_burden, data=merged_rna)
roc.all.mismatch_alleles_expression <- roc(Combined_response~Epitope_by_mismatches_alleles_and_expression_burden, data=merged_rna)

roc.mel.mut <- roc(Combined_response~Total_mutations, data=merged_rna[merged_rna$Disease == 'melanoma',])
roc.mel.tvb <- roc(Combined_response~TVB, data=merged_rna[merged_rna$Disease == 'melanoma',])
roc.mel.jx <- roc(Combined_response~Jx_burden, data=merged_rna[merged_rna$Disease == 'melanoma',])
roc.mel.ri <- roc(Combined_response~Intron_burden, data=merged_rna[merged_rna$Disease == 'melanoma',])
roc.mel.ri_eps <- roc(Combined_response~Binding_intron_epitope_burden, data=merged_rna[merged_rna$Disease == 'melanoma',])
roc.mel.expression <- roc(Combined_response~Epitope_by_expression_burden, data=merged_rna[merged_rna$Disease == 'melanoma',])
roc.mel.mismatch_expression <- roc(Combined_response~Epitope_by_mismatch_and_expression_burden, data=merged_rna[merged_rna$Disease == 'melanoma',])
roc.mel.allele_expression <- roc(Combined_response~Epitope_by_allele_and_expression_burden, data=merged_rna[merged_rna$Disease == 'melanoma',])
roc.mel.mismatch_alleles_expression <- roc(Combined_response~Epitope_by_mismatches_alleles_and_expression_burden, data=merged_rna[merged_rna$Disease == 'melanoma',])

roc.rcc.mut <- roc(Combined_response~Total_mutations, data=merged_rna[merged_rna$Disease == 'RCC',])
roc.rcc.tvb <- roc(Combined_response~TVB, data=merged_rna[merged_rna$Disease == 'RCC',])
roc.rcc.jx <- roc(Combined_response~Jx_burden, data=merged_rna[merged_rna$Disease == 'RCC',])
roc.rcc.ri <- roc(Combined_response~Intron_burden, data=merged_rna[merged_rna$Disease == 'RCC',])
roc.rcc.ri_eps <- roc(Combined_response~Binding_intron_epitope_burden, data=merged_rna[merged_rna$Disease == 'RCC',])
roc.rcc.expression <- roc(Combined_response~Epitope_by_expression_burden, data=merged_rna[merged_rna$Disease == 'RCC',])
roc.rcc.mismatch_expression <- roc(Combined_response~Epitope_by_mismatch_and_expression_burden, data=merged_rna[merged_rna$Disease == 'RCC',])
roc.rcc.allele_expression <- roc(Combined_response~Epitope_by_allele_and_expression_burden, data=merged_rna[merged_rna$Disease == 'RCC',])
roc.rcc.mismatch_alleles_expression <- roc(Combined_response~Epitope_by_mismatches_alleles_and_expression_burden, data=merged_rna[merged_rna$Disease == 'RCC',])

roc_list_all <- list(roc.all.mut, roc.all.tvb, roc.all.jx, roc.all.ri, roc.all.ri_eps, roc.all.expression, 
                     roc.all.mismatch_expression, roc.all.allele_expression, roc.all.mismatch_alleles_expression)
mel_cancer_roc_list <- list(roc.mel.mut, roc.mel.tvb, roc.mel.jx, roc.mel.ri, roc.mel.ri_eps, roc.mel.expression, 
                            roc.mel.mismatch_expression, roc.mel.allele_expression, roc.mel.mismatch_alleles_expression)
rcc_cancer_roc_list <- list(roc.rcc.mut, roc.rcc.tvb, roc.rcc.jx, roc.rcc.ri, roc.rcc.ri_eps, roc.rcc.expression, 
                            roc.rcc.mismatch_expression, roc.rcc.allele_expression, roc.rcc.mismatch_alleles_expression)

color_set <- c(brewer.pal(9, "Set3"))
label_set <- c("Mutation burden", "TVB", "Neojx", "RI", "RI epitopes", "Patient-specific expression (E)", 
               "Amino acid mismatches (M) * E", "Alleles presenting (A) * E", "M*A*E")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "RNA feature burden",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "RNA feature burden",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "RNA feature burden",
                     labels = label_set,
                     values = color_set)+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC'),
                             N=c(106, 89, 17),
                             TMB=c(sprintf("%.3f", round(roc.all.mut$auc, 3)),
                                   sprintf("%.3f", round(roc.mel.mut$auc, 3)),
                                   sprintf("%.3f", round(roc.rcc.mut$auc, 3))),
                             TVB=c(sprintf("%.3f", round(roc.all.tvb$auc, 3)),
                                   sprintf("%.3f", round(roc.mel.tvb$auc, 3)),
                                   sprintf("%.3f", round(roc.rcc.tvb$auc, 3))),
                             Neojx=c(sprintf("%.3f", round(roc.all.jx$auc, 3)),
                                     sprintf("%.3f", round(roc.mel.jx$auc, 3)),
                                     sprintf("%.3f", round(roc.rcc.jx$auc, 3))),
                             RI=c(sprintf("%.3f", round(roc.all.ri$auc, 3)),
                                  sprintf("%.3f", round(roc.mel.ri$auc, 3)),
                                  sprintf("%.3f", round(roc.rcc.ri$auc, 3))),
                             RI_eps=c(sprintf("%.3f", round(roc.all.ri_eps$auc, 3)),
                                      sprintf("%.3f", round(roc.mel.ri_eps$auc, 3)),
                                      sprintf("%.3f", round(roc.rcc.ri_eps$auc, 3))),
                             Expression=c(sprintf("%.3f", round(roc.all.expression$auc, 3)),
                                          sprintf("%.3f", round(roc.mel.expression$auc, 3)),
                                          sprintf("%.3f", round(roc.rcc.expression$auc, 3))),
                             Mismatches=c(sprintf("%.3f", round(roc.all.mismatch_expression$auc, 3)),
                                          sprintf("%.3f", round(roc.mel.mismatch_expression$auc, 3)),
                                          sprintf("%.3f", round(roc.rcc.mismatch_expression$auc, 3))),
                             Alleles=c(sprintf("%.3f", round(roc.all.allele_expression$auc, 3)),
                                       sprintf("%.3f", round(roc.mel.allele_expression$auc, 3)),
                                       sprintf("%.3f", round(roc.rcc.allele_expression$auc, 3))),
                             All=c(sprintf("%.3f", round(roc.all.mismatch_alleles_expression$auc, 3)),
                                   sprintf("%.3f", round(roc.mel.mismatch_alleles_expression$auc, 3)),
                                   sprintf("%.3f", round(roc.rcc.mismatch_alleles_expression$auc, 3))),
                             stringsAsFactors = F)
colnames(table_contents) <-  c("Cancer type", "N", "TMB", "TVB", "Jx", "RI", "RI epitopes",
                               "E", "M*E", "A*E", "M*A*E")

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", "gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", "gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# list of cells to bold (r,c) 
cell_index <- list(c(3,3), c(4,3), c(5,5))

for (cell in cell_index){
  summary_tab$grobs[find_cell(summary_tab, cell[1], cell[2])][[1]][["gp"]] <- gpar(fontface="bold")
}

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"), ncol= 2),
             textGrob("AUC by RNA Feature Burden and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")))


##### NetCTLpan Epitopes - Supplementary Figure 10 #####

### netCTLpan epitopes ###

# Create models
roc.all.netctl <- roc(Combined_response~NetCTLpan_epitopes, data=immunorx_netctl)
roc.mel.netctl <- roc(Combined_response~NetCTLpan_epitopes, data=immunorx_netctl[immunorx_netctl$Disease == 'melanoma',])
roc.rcc.netctl <- roc(Combined_response~NetCTLpan_epitopes, data=immunorx_netctl[immunorx_netctl$Disease == 'RCC',])
roc.nsclc.netctl <- roc(Combined_response~NetCTLpan_epitopes, data=immunorx_netctl[immunorx_netctl$Disease == 'NSCLC',])

roc_list_all <- list(roc.all.netctl)
mel_cancer_roc_list <- list(roc.mel.netctl)
rcc_cancer_roc_list <- list(roc.rcc.netctl)
nsclc_cancer_roc_list <- list(roc.nsclc.netctl)
color_set <- c("grey50")
label_set <- c("NetCTLpan epitopes")

# Create individual plots by disease
gg.all <- ggroc(roc_list_all,
                linetype = 1, size=1, alpha=.9,
                legacy.axes = F)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("All Cancers")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.mel <- ggroc(mel_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("Melanoma")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.rcc <- ggroc(rcc_cancer_roc_list,
                linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("RCC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

gg.nsclc <- ggroc(nsclc_cancer_roc_list,
                  linetype = 1, size=1, alpha=.9)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")+
  theme_classic()+
  scale_color_manual(name = "",
                     labels = label_set,
                     values = color_set)+
  ggtitle("NSCLC")+theme(legend.position = "bottom")+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(face="bold", size = 16))

# Create summary table
table_contents <- data.frame(Cancer_type=c('All', 'Melanoma', 'RCC', 'NSCLC'),
                             N=c(431, 302, 57, 34), 
                             Processed_epitopes=c(sprintf("%.3f", round(roc.all.netctl$auc, 3)),
                                        sprintf("%.3f", round(roc.mel.netctl$auc, 3)),
                                        sprintf("%.3f", round(roc.rcc.netctl$auc, 3)),
                                        sprintf("%.3f", round(roc.nsclc.netctl$auc, 3))),
                             stringsAsFactors = F)
colnames(table_contents)[1] <- "Cancer type" 
colnames(table_contents)[3] <- "Processed epitopes" 

tt_custom <- ttheme_minimal(
  colhead=list(bg_params=list(fill=c("gray92", "gray92", color_set, col = NA)),
               fg_params=list(col=c("gray92", "gray92", color_set))),
  core=list(fg_params=list(fontface=c(2,rep(1,nrow(table_contents))))))

summary_tab <- tableGrob(rbind(colnames(table_contents),table_contents), rows = NULL, theme = tt_custom)
separators <- replicate(ncol(summary_tab) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
summary_tab <- gtable::gtable_add_grob(summary_tab, grobs = separators,
                                       t = 2, b = nrow(summary_tab), l = seq_len(ncol(summary_tab)-1)+1)
summary_tab <- gtable::gtable_add_grob(summary_tab, 
                                       grobs = segmentsGrob( # line across the bottom
                                         x0 = unit(0,"npc"),
                                         y0 = unit(0,"npc"),
                                         x1 = unit(1,"npc"),
                                         y1 = unit(0,"npc"),
                                         gp = gpar(lwd = 3.0)),
                                       t = 2, b = 2, l = 1, r = ncol(summary_tab))

# Arrange plots
grid.arrange(arrangeGrob(gg.all+theme(legend.position = "none"),
                         gg.mel+theme(legend.position = "none"),
                         gg.rcc+theme(legend.position = "none"),
                         gg.nsclc+theme(legend.position = "none"),ncol= 2),
             textGrob("AUC by Processed Epitopes and Cancer Type",gp=gpar(fontsize=16), just = "top"),
             summary_tab, ncol = 1,  heights=unit(c(120,10,50), c("mm")))

