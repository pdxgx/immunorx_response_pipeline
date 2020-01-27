#### Set working directory and load in the data
setwd("YOUR DIRECTORY HERE") ## UPDATE THIS DIRECTORY 

#### Load in the necessary libraries
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)

#### User defined functions
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)}


#### Load in data for response and mb coverage based on varying thresholds
dat <- fread("merged_rna_data.tsv", stringsAsFactors = F)
dat$TVB <- dat$Total_mutations + dat$Jx_burden + dat$Intron_burden

# with update use MHCnuggets epitope count not total coprehensive neoepitope count
dat$Total_comprehensive_neoepitopes <- dat$MHCnuggets_eps

# check to see how many patients are missing the info necessary for response modeling, see which ones
table(!(is.na(dat$Total_mutations) | is.na(dat$Coverage) | is.na(dat$Combined_response)))

#### define basic data to be used for response modeling
dat <- dat[which(!(is.na(dat$Total_mutations) | is.na(dat$Coverage) | is.na(dat$Combined_response))),] # basic info needed for all the models to be tested
dat$Combined_response <- as.factor(dat$Combined_response)

#### Main Question: Does total mutational burden (mutations/mb) predict immunotherapy response
## Do does the answer change based on the variant caller used?
## How does this vary across cancer types?
variant_callers <- c("Consensus", "MuSE", "Mutect", "Pindel", "RADIA", "SomaticSniper", "Varscan")


## Melanoma specifically
mel_cancer_response_list <- list()
# consensus
mel_cancer_response_list$all.mod.consensus_muts <- glm(as.factor(Combined_response)~log2(Total_mutations/Coverage), data = subset(dat, dat$Disease == "melanoma"), family = 'binomial')

# tool specific models of response
mel_cancer_response_list$all.mod.muse_muts <- glm(as.factor(Combined_response)~log2(Muse_variants/Coverage), data = subset(dat, dat$Disease == "melanoma"), family = 'binomial')
mel_cancer_response_list$all.mod.mutect_muts <- glm(as.factor(Combined_response)~log2(Mutect_variants/Coverage), data = subset(dat, dat$Disease == "melanoma"), family = 'binomial')
mel_cancer_response_list$all.mod.pindel_muts <- glm(as.factor(Combined_response)~log2(Pindel_variants/Coverage), data = subset(dat, dat$Disease == "melanoma"), family = 'binomial')
mel_cancer_response_list$all.mod.radia_muts <- glm(as.factor(Combined_response)~log2(Radia_variants/Coverage), data = subset(dat, dat$Disease == "melanoma"), family = 'binomial')
mel_cancer_response_list$all.mod.somaticsniper_muts <- glm(as.factor(Combined_response)~log2(Somaticsniper_variants/Coverage), data = subset(dat, dat$Disease == "melanoma"), family = 'binomial')
mel_cancer_response_list$all.mod.varscan_muts <- glm(as.factor(Combined_response)~log2(Varscan_variants/Coverage), data = subset(dat, dat$Disease == "melanoma"), family = 'binomial')

# show model metrics
data.frame(model = variant_callers,
           OR = round(unlist(lapply(mel_cancer_response_list, function(x){return(exp(summary(x)$coefficients[,"Estimate"][which(names(summary(x)$coefficients[,"Estimate"]) != "(Intercept)")]))})),3), 
           AIC = round(unlist(lapply(mel_cancer_response_list, function(x){return(summary(x)$aic)}))),
           pVal = round(unlist(lapply(mel_cancer_response_list, function(x){return(anova(x, test = "Chisq")$`Pr(>Chi)`[!is.na(anova(x, test = "Chisq")$`Pr(>Chi)`)])})),3),
           row.names = NULL)


## RCC
rcc_cancer_response_list <- list()

# consensus
rcc_cancer_response_list$all.mod.consensus_muts <- glm(as.factor(Combined_response)~log2(Total_mutations/Coverage), data = subset(dat, dat$Disease == "RCC"), family = 'binomial')

# tool specific models of response
rcc_cancer_response_list$all.mod.muse_muts <- glm(as.factor(Combined_response)~log2(Muse_variants/Coverage), data = subset(dat, dat$Disease == "RCC"), family = 'binomial')
rcc_cancer_response_list$all.mod.mutect_muts <- glm(as.factor(Combined_response)~log2(Mutect_variants/Coverage), data = subset(dat, dat$Disease == "RCC"), family = 'binomial')
rcc_cancer_response_list$all.mod.pindel_muts <- glm(as.factor(Combined_response)~log2(Pindel_variants/Coverage), data = subset(dat, dat$Disease == "RCC"), family = 'binomial')
rcc_cancer_response_list$all.mod.radia_muts <- glm(as.factor(Combined_response)~log2(Radia_variants/Coverage), data = subset(dat, dat$Disease == "RCC"), family = 'binomial')
rcc_cancer_response_list$all.mod.somaticsniper_muts <- glm(as.factor(Combined_response)~log2(Somaticsniper_variants/Coverage), data = subset(dat, dat$Disease == "RCC"), family = 'binomial')
rcc_cancer_response_list$all.mod.varscan_muts <- glm(as.factor(Combined_response)~log2(Varscan_variants/Coverage), data = subset(dat, dat$Disease == "RCC"), family = 'binomial')

data.frame(model = variant_callers,
           OR = round(unlist(lapply(rcc_cancer_response_list, function(x){return(exp(summary(x)$coefficients[,"Estimate"][which(names(summary(x)$coefficients[,"Estimate"]) != "(Intercept)")]))})),3), 
           AIC = round(unlist(lapply(rcc_cancer_response_list, function(x){return(summary(x)$aic)}))),
           pVal = round(unlist(lapply(rcc_cancer_response_list, function(x){return(anova(x, test = "Chisq")$`Pr(>Chi)`[!is.na(anova(x, test = "Chisq")$`Pr(>Chi)`)])})),3),
           row.names = NULL
)

## NSCLC
nsclc_cancer_response_list <- list()

# consensus
nsclc_cancer_response_list$all.mod.consensus_muts <- glm(as.factor(Combined_response)~log2(Total_mutations/Coverage), data = subset(dat, dat$Disease == "NSCLC"), family = 'binomial')

# tool specific models of response
nsclc_cancer_response_list$all.mod.muse_muts <- glm(as.factor(Combined_response)~log2(Muse_variants/Coverage), data = subset(dat, dat$Disease == "NSCLC"), family = 'binomial')
nsclc_cancer_response_list$all.mod.mutect_muts <- glm(as.factor(Combined_response)~log2(Mutect_variants/Coverage), data = subset(dat, dat$Disease == "NSCLC"), family = 'binomial')
nsclc_cancer_response_list$all.mod.pindel_muts <- glm(as.factor(Combined_response)~log2(Pindel_variants/Coverage), data = subset(dat, dat$Disease == "NSCLC"), family = 'binomial')
nsclc_cancer_response_list$all.mod.radia_muts <- glm(as.factor(Combined_response)~log2(Radia_variants/Coverage), data = subset(dat, dat$Disease == "NSCLC"), family = 'binomial')
nsclc_cancer_response_list$all.mod.somaticsniper_muts <- glm(as.factor(Combined_response)~log2(Somaticsniper_variants/Coverage), data = subset(dat, dat$Disease == "NSCLC"), family = 'binomial')
nsclc_cancer_response_list$all.mod.varscan_muts <- glm(as.factor(Combined_response)~log2(Varscan_variants/Coverage), data = subset(dat, dat$Disease == "NSCLC"), family = 'binomial')

# show model metrics
data.frame(model = variant_callers,
           OR = round(unlist(lapply(nsclc_cancer_response_list, function(x){return(exp(summary(x)$coefficients[,"Estimate"][which(names(summary(x)$coefficients[,"Estimate"]) != "(Intercept)")]))})),3), 
           AIC = round(unlist(lapply(nsclc_cancer_response_list, function(x){return(summary(x)$aic)}))),
           pVal = round(unlist(lapply(nsclc_cancer_response_list, function(x){return(anova(x, test = "Chisq")$`Pr(>Chi)`[!is.na(anova(x, test = "Chisq")$`Pr(>Chi)`)])})),3),
           row.names = NULL
)

#### predict response for each 25th and 75th percentiles from consensus models
# create consistent labeling for treatments
dat$treatment <- rep(NA, dim(dat)[1])
table(apply(dat[,c("aPD1_treatment", "aPDL1_treatment", "aCTLA4_treatment", "Other_treatment")], 1, sum)) # look at N per treatment
dat$treatment[dat$aCTLA4_treatment==1] <- "aCTLA4"
dat$treatment[dat$aPD1_treatment==1] <- "aPD1"
dat$treatment[dat$aPDL1_treatment==1] <- "aPD1"
dat$treatment[dat$Other_treatment==1] <- "other"
dat$treatment[which(apply(dat[,c("aPD1_treatment", "aPDL1_treatment", "aCTLA4_treatment", "Other_treatment")], 1, sum)>1)] <- "multi"
table(dat$treatment)

# make new model for mel that includes treatment type to separate out the two
mel_consensus_muts_tx <- glm(as.factor(Combined_response)~log2(Total_mutations/Coverage)+treatment, 
                             data = subset(dat, dat$Disease == "melanoma" & (dat$treatment == "aCTLA4" | dat$treatment == "aPD1")), 
                             family = 'binomial')

# create a data frame for holding predictions
mel_prediction_table <- data.frame(Total_mutations = rep(quantile(dat$Total_mutations[dat$Disease == "melanoma"]/dat$Coverage[dat$Disease == "melanoma"], c(.25,.75)),2), 
                                   treatment = rep(c("aCTLA4","aPD1"),each = 2))
mel_prediction_table$Coverage <- c(1,1,1,1)

# predict odds and append to the table
mel_preds <- unlist(lapply(list(mel_consensus_muts_tx), function(x){return(predict(x, mel_prediction_table))}))
mel_prediction_table$odds <- exp(mel_preds)
mel_prediction_table


# make predictions for NSCLC
nsclc_prediction_table <- data.frame(Total_mutations = quantile(dat$Total_mutations[dat$Disease == "NSCLC"]/dat$Coverage[dat$Disease == "NSCLC"], c(.25,.75)))
nsclc_prediction_table$Coverage <- c(1,1)
nsclc_preds <- unlist(lapply(nsclc_cancer_response_list[1], function(x){return(predict(x, nsclc_prediction_table))}))
nsclc_prediction_table$odds <- exp(nsclc_preds)
nsclc_prediction_table


# make predictions for RCC
rcc_prediction_table <- data.frame(Total_mutations = quantile(dat$Total_mutations[dat$Disease == "RCC"]/dat$Coverage[dat$Disease == "RCC"], c(.25,.75)))
rcc_prediction_table$Coverage <- c(1,1)
rcc_preds <- unlist(lapply(rcc_cancer_response_list[1], function(x){return(predict(x, rcc_prediction_table))}))
rcc_prediction_table$odds <- exp(rcc_preds)
rcc_prediction_table


# aggregate summary p-values for each model
model_pvals <- unlist(lapply(c(list(mel_consensus_muts_tx), 
                               nsclc_cancer_response_list[1], 
                               rcc_cancer_response_list[1]), 
                             function(x){summary(x)$coefficients["log2(Total_mutations/Coverage)","Pr(>|z|)"]}))
model_pvals


#### model response based on neoepitope burden and tumor variant burden
# create neoepitope models
neo_mods <- list()
neo_mods$melanoma <-  glm(as.factor(Combined_response)~log2(Total_comprehensive_neoepitopes/Coverage)+treatment, data = subset(dat, dat$Disease == "melanoma" & (dat$treatment == "aCTLA4" | dat$treatment == "aPD1")), family = 'binomial')
neo_mods$nsclc <-  glm(as.factor(Combined_response)~log2(Total_comprehensive_neoepitopes/Coverage), data = subset(dat, dat$Disease == "NSCLC"), family = 'binomial')
neo_mods$RCC <-  glm(as.factor(Combined_response)~log2(Total_comprehensive_neoepitopes/Coverage), data = subset(dat, dat$Disease == "RCC"), family = 'binomial')

# create TVB models
TVB_mods <- list()
TVB_mods$melanoma <- glm(as.factor(Combined_response)~log2(TVB)+treatment, data = subset(dat, dat$Disease == "melanoma" & !is.na(dat$TVB) & (dat$treatment == "aCTLA4" | dat$treatment == "aPD1")), family = 'binomial')
TVB_mods$RCC <- glm(as.factor(Combined_response)~log2(TVB), data = subset(dat, dat$Disease == "RCC" & !is.na(dat$TVB)), family = 'binomial')

# create table structure
neo_prediction_table <- data.frame(Total_comprehensive_neoepitopes = c(rep(quantile(dat$Total_comprehensive_neoepitopes[dat$Disease == "melanoma"]/dat$Coverage[dat$Disease == "melanoma"], c(.25,.75)), 2),
                                                                       quantile(dat$Total_comprehensive_neoepitopes[dat$Disease == "NSCLC"]/dat$Coverage[dat$Disease == "NSCLC"], c(.25,.75)),
                                                                       quantile(dat$Total_comprehensive_neoepitopes[dat$Disease == "RCC"]/dat$Coverage[dat$Disease == "RCC"], c(.25,.75))))

# populate table
neo_prediction_table$Coverage <- rep(1, dim(neo_prediction_table)[1])
neo_prediction_table$treatment <- c("aCTLA4","aCTLA4","aPD1","aPD1",NA,NA,NA,NA)
neo_preds <- rep(NA, dim(neo_prediction_table)[1])
neo_preds[1:4] <- unlist(lapply(neo_mods[1], function(x){return(predict(x, neo_prediction_table[1:4,]))}))
neo_preds[5:6] <- unlist(lapply(neo_mods[2], function(x){return(predict(x, neo_prediction_table[5:6,]))}))
neo_preds[7:8] <- unlist(lapply(neo_mods[3], function(x){return(predict(x, neo_prediction_table[7:8,]))}))

neo_prediction_table$odds <- exp(neo_preds)
neo_prediction_table$prob <- round(neo_prediction_table$odds/(1+neo_prediction_table$odds), 3)

neo_mod_summaries <- lapply(neo_mods, summary)
neo_model_pvals <- unlist(lapply(neo_mod_summaries, function(x){return(x$coefficients["log2(Total_comprehensive_neoepitopes/Coverage)","Pr(>|z|)"])}))

## repeat for TVB
tvb_prediction_table <- data.frame(TVB = c(rep(quantile(dat$TVB[dat$Disease == "melanoma"], c(.25,.75), na.rm = T), 2),
                                           quantile(dat$TVB[dat$Disease == "RCC"], c(.25,.75), na.rm = T)))

tvb_prediction_table$treatment <- c("aCTLA4","aCTLA4","aPD1","aPD1",NA,NA)
tvb_preds <- rep(NA, dim(tvb_prediction_table)[1])
tvb_preds[1:4] <- unlist(lapply(TVB_mods[1], function(x){return(predict(x, tvb_prediction_table[1:4,]))}))
tvb_preds[5:6] <- unlist(lapply(TVB_mods[2], function(x){return(predict(x, tvb_prediction_table[5:6,]))}))


tvb_prediction_table$odds <- exp(tvb_preds)
tvb_prediction_table$prob <- round(tvb_prediction_table$odds/(1+tvb_prediction_table$odds), 3)

tvb_mod_summaries <- lapply(TVB_mods, summary)
tvb_model_pvals <- unlist(lapply(tvb_mod_summaries, function(x){return(x$coefficients["log2(TVB)","Pr(>|z|)"])}))


#### create tables to summarize model results
# TVB
tvb_response_summary_tab <- data.frame(Cancer_type = c("Melanoma", "RCC"),
                                       aCTLA_N = c(61,NA),
                                       aCTLA_p_25 = c(tvb_prediction_table$prob[1],NA),
                                       aCTLA_p_75 = c(tvb_prediction_table$prob[2],NA),
                                       aPD1_N = c(27,17),
                                       aPD1_p_25 = c(tvb_prediction_table$prob[3],
                                                     tvb_prediction_table$prob[5]),
                                       aPD1_p_75 = c(tvb_prediction_table$prob[4],
                                                     tvb_prediction_table$prob[6]),
                                       pval = round(tvb_model_pvals, 4))

tvb_response_summary_tab

# TMB
response_summary_tab <- data.frame(Cancer_type = c("Melanoma", "NSCLC","RCC"),
                                   aCTLA_N = c(195,NA,NA),
                                   aCTLA_odds_25 = c(mel_prediction_table$odds[1],NA,NA),
                                   aCTLA_odds_75 = c(mel_prediction_table$odds[2],NA,NA),
                                   aPD1_N = c(50,33,50),
                                   aPD1_odds_25 = c(mel_prediction_table$odds[3],
                                                    nsclc_prediction_table$odds[1],
                                                    rcc_prediction_table$odds[1]),
                                   aPD1_odds_75 = c(mel_prediction_table$odds[4],
                                                    nsclc_prediction_table$odds[2],
                                                    rcc_prediction_table$odds[2]),
                                   pval = round(model_pvals, 4))
response_summary_tab


### Neoepitope
neo_response_summary_tab <- data.frame(Cancer_type = c("Melanoma", "NSCLC", "RCC"),
                                       aCTLA_N = c(195,NA,NA),
                                       aCTLA_p_25 = c(neo_prediction_table$prob[1],NA, NA),
                                       aCTLA_p_75 = c(neo_prediction_table$prob[2],NA, NA),
                                       aPD1_N = c(50,33,50),
                                       aPD1_p_25 = c(neo_prediction_table$prob[3],
                                                     neo_prediction_table$prob[5],
                                                     neo_prediction_table$prob[7]),
                                       aPD1_p_75 = c(neo_prediction_table$prob[4],
                                                     neo_prediction_table$prob[6],
                                                     neo_prediction_table$prob[8]),
                                       pval = round(neo_model_pvals, 4))

neo_response_summary_tab


#convert to response probability from odds... Neoepitope and TVB already in probability
response_summary_tab[,c(3,4,6,7)] <- round(response_summary_tab[,c(3,4,6,7)]/(1+response_summary_tab[,c(3,4,6,7)]), 3)
response_summary_tab
colnames(response_summary_tab)[c(3,4,6,7)] <- c("aCTLA_p_25", "aCTLA_p_75", "aPD1_p_25", "aPD1_p_75")
write.csv(response_summary_tab, file = "tmb-resp-summ.csv", quote = F, row.names = F)
write.csv(neo_response_summary_tab, file = "neo-resp-summ.csv", quote = F, row.names = F)
write.csv(tvb_response_summary_tab, file = "tvb-resp-summ.csv", quote = F, row.names = F)

#### add in the survival models, do iterations to show HR at diff. TMB cuts. make table of linear effects and plots w/cuts
# subset to remove those with only response info
table(is.na(dat$OS_event_status) & is.na(dat$Censoring_status))
surv_dat <- dat[which(!(is.na(dat$OS_event_status) & is.na(dat$Censoring_status)) | !(is.na(dat$OS) & is.na(dat$Censoring_days))),]
surv_dat <- subset(surv_dat, surv_dat$Study != "roh") # this study did not use a standard/comparable dosing schedule

# make new var to combine both censoring days and survival days
surv_dat$OS_time <- rep(NA, dim(surv_dat)[1])
# populate variable
surv_dat$OS_time[which(surv_dat$Censoring_status==1)] <- surv_dat$Censoring_days[which(surv_dat$Censoring_status==1)]
surv_dat$OS_time[which(surv_dat$OS_event_status==1)] <- surv_dat$OS[which(surv_dat$OS_event_status==1)]
table(is.na(surv_dat$OS_time)) # some with slightly different coding
surv_dat$OS_time[which(surv_dat$OS_event_status==0 & surv_dat$Censoring_status!=1)] <- surv_dat$OS[which(surv_dat$OS_event_status==0 & surv_dat$Censoring_status!=1)] # appears to be those that made it to end of trial
table(is.na(surv_dat$OS_time)) # a few still missing...
surv_dat[which(is.na(surv_dat$OS_time)),c("Tumor_ID", "Study","Disease", "PFS", "OS", "OS_event_status", "Censoring_status")] # looks like rixvi study is PFS only? double check


## look for melanoma stable cutoff
mel_mut_cuts <- quantile(surv_dat$Total_mutations[which(surv_dat$Disease == "melanoma")]/surv_dat$Coverage[which(surv_dat$Disease == "melanoma")], probs = seq(0,1,.02))
mel_HR_list <- rep(NA, length(mel_mut_cuts))
mel_significance <- rep(NA, length(mel_mut_cuts))

for (i in 1:length(mel_mut_cuts)){
  mod <- coxph(Surv(OS_time, OS_event_status)~(as.numeric(Total_mutations/Coverage) >= mel_mut_cuts[i]), data = subset(surv_dat, surv_dat$Disease == "melanoma"))
  mel_HR_list[i] <- summary(mod)$coefficients[,"exp(coef)"]
  mel_significance[i] <- summary(mod)$coefficients[,"Pr(>|z|)"]
}

# plot
qp_mel <- qplot(seq(0,1,.02)*100, mel_HR_list, color = (mel_significance <= .05))+
  xlim(100,0)+
  ylim(0,1.6)+
  scale_color_manual(values = c("grey50", "red"))+
  theme_classic()+guides(color = FALSE)+
  ggtitle("Melanoma")+
  xlab("TMB cutoff (percentile)")+
  ylab("Cutoff Based Hazard Ratio (above/below)")

## look for RCC cutoff
RCC_mut_cuts <- quantile(surv_dat$Total_mutations[which(surv_dat$Disease == "RCC")]/surv_dat$Coverage[which(surv_dat$Disease == "RCC")], probs = seq(0,1,.02))
RCC_HR_list <- rep(NA, length(RCC_mut_cuts))
RCC_significance <- rep(NA, length(RCC_mut_cuts))

for (i in 1:length(RCC_mut_cuts)){
  mod <- coxph(Surv(OS_time, OS_event_status)~(as.numeric(Total_mutations/Coverage) >= RCC_mut_cuts[i]), data = subset(surv_dat, surv_dat$Disease == "RCC"))
  RCC_HR_list[i] <- summary(mod)$coefficients[,"exp(coef)"]
  RCC_significance[i] <- summary(mod)$coefficients[,"Pr(>|z|)"]
}

# plot
qp_rcc <- qplot(seq(0,1,.02)*100, RCC_HR_list, color = (RCC_significance <= .05))+
  xlim(100,0)+
  ylim(0,1.6)+
  scale_color_manual(values = c("grey50", "red"))+
  theme_classic()+guides(color = FALSE)+
  ggtitle("RCC")+
  xlab("TMB cutoff (percentile)")+
  ylab("Cutoff Based Hazard Ratio (above/below)")


hr_cut_plot <- grid.arrange(qp_mel, qp_rcc)
ggsave("./figures/HR-by-TMB-cut-updated.pdf", plot=hr_cut_plot, width = 220, height = 170, units = "mm")

