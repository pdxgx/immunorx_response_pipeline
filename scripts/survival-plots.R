setwd("YOUR DIRECTORY HERE") ## UPDATE THIS DIRECTORY 

#### Load in the necessary libraries
library(ggplot2)
library(survival)
library(survminer)
library(data.table)
library(gridExtra)
library(grid)
library(RColorBrewer)

#### Load in TCGA data for later
tcga_mel_dat <- fread("skcm_tmb_survival_data.tsv", stringsAsFactors = F)
tcga_rcc_dat <- fread("kirc_tmb_survival_data.tsv", stringsAsFactors = F)

#### Load in data for response and mb coverage based on varying thresholds
dat <- fread("merged_rna_data.tsv", stringsAsFactors = F)
dat$TVB <- dat$Total_mutations + dat$Jx_burden + dat$Intron_burden

# with update use MHCnuggets epitope count not total coprehensive
dat$Total_comprehensive_neoepitopes <- dat$MHCnuggets_eps

surv_dat <- dat[which(!(is.na(dat$OS_event_status) & is.na(dat$Censoring_status)) | !(is.na(dat$OS) & is.na(dat$Censoring_days))),]
surv_dat <- subset(surv_dat, surv_dat$Study != "roh") # this study did not use a standard/comparable dosing schedule

# make new var to combine both censoring days and survival days
surv_dat$OS_time <- rep(NA, dim(surv_dat)[1])
# populate variable
surv_dat$OS_time[which(surv_dat$Censoring_status==1)] <- surv_dat$Censoring_days[which(surv_dat$Censoring_status==1)]
surv_dat$OS_time[which(surv_dat$OS_event_status==1)] <- surv_dat$OS[which(surv_dat$OS_event_status==1)]
surv_dat$OS_time[which(surv_dat$OS_event_status==0 & surv_dat$Censoring_status!=1)] <- surv_dat$OS[which(surv_dat$OS_event_status==0 & surv_dat$Censoring_status!=1)]


### make survival curves based on cutoff
# create marker for the two cutoff groups
surv_dat$mel_cut_status <- rep(NA, dim(surv_dat)[1])
surv_dat$mel_cut_status[which(surv_dat$Disease == "melanoma" & 
                                surv_dat$Total_mutations/surv_dat$Coverage > quantile(surv_dat$Total_mutations[surv_dat$Disease == "melanoma"]/surv_dat$Coverage[surv_dat$Disease == "melanoma"], .8))] <- "High TMB"
surv_dat$mel_cut_status[which(surv_dat$Disease == "melanoma" & 
                                surv_dat$Total_mutations/surv_dat$Coverage <= quantile(surv_dat$Total_mutations[surv_dat$Disease == "melanoma"]/surv_dat$Coverage[surv_dat$Disease == "melanoma"], .8))] <- "Low TMB"

#### compare to TCGA melanoma
tcga_mel_dat <- subset(tcga_mel_dat, (tcga_mel_dat$Stage %in% c("III", "IIIA", "IIIB", "IIIC", "IV")))
# censor at same follow up time as current data
tcga_mel_dat$OS_event[tcga_mel_dat$Overall_survival >= max(surv_dat$OS_time[which(surv_dat$Disease == "melanoma")])] <- 0
tcga_mel_dat$Overall_survival[tcga_mel_dat$Overall_survival >= max(surv_dat$OS_time[which(surv_dat$Disease == "melanoma")])] <- max(surv_dat$OS_time[which(surv_dat$Disease == "melanoma")])

tcga_mel_dat$mel_cut_status <- rep(NA, dim(tcga_mel_dat)[1])
tcga_mel_dat$mel_cut_status[which(tcga_mel_dat$TMB > quantile(tcga_mel_dat$TMB, .8))] <- "TCGA High TMB"
tcga_mel_dat$mel_cut_status[which(tcga_mel_dat$TMB <= quantile(tcga_mel_dat$TMB, .8))] <- "TCGA Low TMB"

mel_plot_data <- data.frame(OS_time = c(surv_dat$OS_time[which(surv_dat$Disease == "melanoma")], tcga_mel_dat$Overall_survival),
                            OS_event_status = c(surv_dat$OS_event_status[which(surv_dat$Disease == "melanoma")], tcga_mel_dat$OS_event),
                            TMB_grouping = c(surv_dat$mel_cut_status[which(surv_dat$Disease == "melanoma")], tcga_mel_dat$mel_cut_status))

mel_plot_data$TMB_grouping <- factor(mel_plot_data$TMB_grouping, levels = c("High TMB", "Low TMB", "TCGA High TMB", "TCGA Low TMB"))

# create marker for the new cutoff groups
surv_dat$rcc_cut_status <- rep(NA, dim(surv_dat)[1])
surv_dat$rcc_cut_status[which(surv_dat$Disease == "RCC" & 
                                surv_dat$Total_mutations/surv_dat$Coverage > quantile(surv_dat$Total_mutations[surv_dat$Disease == "RCC"]/surv_dat$Coverage[surv_dat$Disease == "RCC"], .8))] <- "High TMB"
surv_dat$rcc_cut_status[which(surv_dat$Disease == "RCC" & 
                                surv_dat$Total_mutations/surv_dat$Coverage <= quantile(surv_dat$Total_mutations[surv_dat$Disease == "RCC"]/surv_dat$Coverage[surv_dat$Disease == "RCC"], .8))] <- "Low TMB"

# Compare to TCGA RCC
tcga_rcc_dat <- subset(tcga_rcc_dat, tcga_rcc_dat$Stage == "IV")
# censor at same follow up time as current data
tcga_rcc_dat$OS_event[tcga_rcc_dat$Overall_survival >= max(surv_dat$OS_time[which(surv_dat$Disease == "RCC")])] <- 0
tcga_rcc_dat$Overall_survival[tcga_rcc_dat$Overall_survival >= max(surv_dat$OS_time[which(surv_dat$Disease == "RCC")])] <- max(surv_dat$OS_time[which(surv_dat$Disease == "RCC")])

tcga_rcc_dat$rcc_cut_status <- rep(NA, dim(tcga_rcc_dat)[1])
tcga_rcc_dat$rcc_cut_status[which(tcga_rcc_dat$TMB > quantile(tcga_rcc_dat$TMB, .8))] <- "TCGA High TMB"
tcga_rcc_dat$rcc_cut_status[which(tcga_rcc_dat$TMB <= quantile(tcga_rcc_dat$TMB, .8))] <- "TCGA Low TMB"


rcc_plot_data <- data.frame(OS_time = c(surv_dat$OS_time[which(surv_dat$Disease == "RCC")], tcga_rcc_dat$Overall_survival),
                            OS_event_status = c(surv_dat$OS_event_status[which(surv_dat$Disease == "RCC")], tcga_rcc_dat$OS_event),
                            TMB_grouping = c(surv_dat$rcc_cut_status[which(surv_dat$Disease == "RCC")], tcga_rcc_dat$rcc_cut_status))

rcc_plot_data$TMB_grouping <- factor(rcc_plot_data$TMB_grouping, levels = c("High TMB", "Low TMB", "TCGA High TMB", "TCGA Low TMB"))


f.mel <- survfit(Surv(OS_time, OS_event_status)~TMB_grouping, data=mel_plot_data)
survP_mel <- ggsurvplot(f.mel, 
                        conf.int = F, 
                        risk.table = T,
                        legend.title = "",
                        font.x = 16,
                        font.y = 16,
                        font.tickslab = 14,
                        font.legend = 12,
                        fontsize = 5,
                        font.title = 16,
                        font.subtitle = 16,
                        legend.labs = c("High TMB, +ICI", "Low TMB, +ICI", "High TMB, -ICI", "Low TMB, -ICI"), 
                        palette = c("red","blue","grey50","lightgrey"),
                        pval = F,
                        break.time.by = 500)

f.rcc <- survfit(Surv(OS_time, OS_event_status)~TMB_grouping, data=rcc_plot_data)
survP_rcc <- ggsurvplot(f.rcc, 
                        conf.int = F, 
                        risk.table = T,
                        legend.title = "",
                        font.x = 16,
                        font.y = 16,
                        font.tickslab = 14,
                        font.legend = 12,
                        fontsize = 5,
                        font.title = 16,
                        font.subtitle = 16,
                        legend.labs = c("High TMB, +ICI", "Low TMB, +ICI", "High TMB, -ICI", "Low TMB, -ICI"),
                        palette = c("red","blue", "grey50", "lightgrey"),
                        pval = F, 
                        break.time.by = 200)


survP_mel$table <- survP_mel$table + theme_cleantable()
survP_mel$plot <- survP_mel$plot + xlab("Time (Days)")
survP_mel

survP_rcc$table <- survP_rcc$table + theme_cleantable()
survP_rcc$plot <- survP_rcc$plot + xlab("Time (Days)")
survP_rcc

#### repeat for neoepitopes
surv_dat$mel_neo_cut_status <- rep(NA, dim(surv_dat)[1])
surv_dat$mel_neo_cut_status[which(surv_dat$Disease == "melanoma" & 
                                    surv_dat$Total_comprehensive_neoepitopes/surv_dat$Coverage > quantile(surv_dat$Total_comprehensive_neoepitopes[surv_dat$Disease == "melanoma"]/surv_dat$Coverage[surv_dat$Disease == "melanoma"], .8))] <- "High Neoepitopes"
surv_dat$mel_neo_cut_status[which(surv_dat$Disease == "melanoma" & 
                                    surv_dat$Total_comprehensive_neoepitopes/surv_dat$Coverage <= quantile(surv_dat$Total_comprehensive_neoepitopes[surv_dat$Disease == "melanoma"]/surv_dat$Coverage[surv_dat$Disease == "melanoma"], .8))] <- "Low Neoepitopes"
surv_dat$mel_neo_cut_status[which(surv_dat$Disease == "melanoma")] <- factor(surv_dat$mel_neo_cut_status[which(surv_dat$Disease == "melanoma")], levels = c("High Neoepitopes", "Low Neoepitopes"))

# create marker for the new cutoff groups
surv_dat$rcc_neo_cut_status <- rep(NA, dim(surv_dat)[1])
surv_dat$rcc_neo_cut_status[which(surv_dat$Disease == "RCC" & 
                                    surv_dat$Total_comprehensive_neoepitopes/surv_dat$Coverage > quantile(surv_dat$Total_comprehensive_neoepitopes[surv_dat$Disease == "RCC"]/surv_dat$Coverage[surv_dat$Disease == "RCC"], .8))] <- "High Neoepitopes"
surv_dat$rcc_neo_cut_status[which(surv_dat$Disease == "RCC" & 
                                    surv_dat$Total_comprehensive_neoepitopes/surv_dat$Coverage <= quantile(surv_dat$Total_comprehensive_neoepitopes[surv_dat$Disease == "RCC"]/surv_dat$Coverage[surv_dat$Disease == "RCC"], .8))] <- "Low Neoepitopes"
surv_dat$rcc_neo_cut_status[which(surv_dat$Disease == "RCC")] <- factor(surv_dat$rcc_neo_cut_status[which(surv_dat$Disease == "RCC")], levels = c("High Neoepitopes", "Low Neoepitopes"))

f.neo_mel <- survfit(Surv(OS_time, OS_event_status)~mel_neo_cut_status, data=subset(surv_dat, surv_dat$Disease == "melanoma"))
survP_neo_mel <- ggsurvplot(f.neo_mel, 
                            legend.title="",
                            conf.int = F, 
                            risk.table = T,
                            font.x = 16,
                            font.y = 16,
                            font.tickslab = 14,
                            font.legend = 12,
                            fontsize = 5,
                            font.title = 16,
                            font.subtitle = 16,
                            legend.labs = c("High Neoepitopes, +ICI","Low Neoepitopes, +ICI"), 
                            palette = c("red","blue"),
                            pval = F,
                            break.time.by = 500)

f.neo_rcc <- survfit(Surv(OS_time, OS_event_status)~rcc_neo_cut_status, data=subset(surv_dat, surv_dat$Disease == "RCC"))
survP_neo_rcc <- ggsurvplot(f.neo_rcc,
                            legend.title = "",
                            conf.int = F, 
                            risk.table = T,
                            font.x = 16,
                            font.y = 16,
                            font.tickslab = 14,
                            font.legend = 12,
                            fontsize = 5,
                            font.title = 16,
                            font.subtitle = 16,
                            legend.labs = c("High Neoepitopes, +ICI","Low Neoepitopes, +ICI"),
                            palette = c("red","blue"),
                            pval = F, 
                            break.time.by = 200)


survP_neo_mel$table <- survP_neo_mel$table + theme_cleantable()
survP_neo_mel$plot <- survP_neo_mel$plot + xlab("Time (Days)")
survP_neo_mel
survP_neo_rcc$table <- survP_neo_rcc$table + theme_cleantable()
survP_neo_rcc$plot <- survP_neo_rcc$plot + xlab("Time (Days)")
survP_neo_rcc

TVB_thresh <- quantile(surv_dat$TVB[surv_dat$Disease == "melanoma"], .8, na.rm = T)

surv_dat$mel_TVB_cut_status <- rep(NA, dim(surv_dat)[1])
surv_dat$mel_TVB_cut_status[surv_dat$Disease == "melanoma" & surv_dat$TVB > TVB_thresh] <- "High TVB"
surv_dat$mel_TVB_cut_status[surv_dat$Disease == "melanoma" & surv_dat$TVB <= TVB_thresh] <- "Low TVB"
surv_dat$mel_TVB_cut_status[surv_dat$Disease == "melanoma"] <- factor(surv_dat$mel_TVB_cut_status[surv_dat$Disease == "melanoma"], levels = c("High TVB", "Low TVB"))

f.mel_TVB <- survfit(Surv(OS_time, OS_event_status)~mel_TVB_cut_status, data=subset(surv_dat, surv_dat$Disease == "melanoma" & !is.na(surv_dat$TVB)))
survP_TVB_mel <- ggsurvplot(f.mel_TVB, 
                            legend.title = "",
                            conf.int = F, 
                            risk.table = T,
                            font.x = 16,
                            font.y = 16,
                            font.tickslab = 14,
                            font.legend = 12,
                            fontsize = 5,
                            font.title = 16,
                            font.subtitle = 16,
                            legend.labs = c("High TVB, +ICI","Low TVB, +ICI"), 
                            palette = c("red","blue"),
                            pval = F,
                            break.time.by = 500)

survP_TVB_mel$table <- survP_TVB_mel$table+theme_cleantable()
survP_TVB_mel$plot <- survP_TVB_mel$plot+xlab("Time (Days)")
survP_TVB_mel

