library(magrittr)
library(kma)

# Tumor:

base_dir <- '[YOUR DIRECTORY CONTAINING TUMOR EXPRESS SUBDIRECTORIES HERE]' ## UPDATE THIS PATH
xprs_fnames <- Sys.glob(file.path(base_dir, "*/results.xprs"))
xprs_fnames

sample_names <- sub(file.path(base_dir), "", xprs_fnames) %>% sub("results.xprs", "", .) %>% gsub("/", "", .)
sample_names

condition_names <- rep('tumor', length(sample_names2))

xprs <- read_express(xprs_fnames, sample_names, condition_names)

intron_to_trans <- data.table::fread("[YOUR DIRECTORY HERE]/intron_to_transcripts.txt", data.table = FALSE) ## UPDATE THIS PATH

ir <- newIntronRetention(xprs$tpm, intron_to_trans, xprs$condition, xprs$uniq_counts)

n_filter <- 0.25
ir_filtered_tpm <- filter_low_tpm(ir, 1, n_filter) # Only include transcripts that are expressed at a tpm >= 1 in > n_filter % of samples
ir_filtered_psi <- filter_perfect_psi(ir_filtered_tpm) # Remove transcripts in which every sample has 0 or 1 retention after rounding (artifact)
ir_filtered_frags <- filter_low_frags(ir_filtered_psi, 5, n_filter) # Only include transcripts that have at least 5 unique counts in > n_filter % of samples

write.table(ir_filtered_frags$unique_counts, 
            file='[YOUR DIRECTORY HERE]/intron_retention_read_counts.tsv', ## UPDATE THIS PATH
            sep='\t', quote=F, col.names = T, row.names = F)

# Melanocytes:

base_dir2 <- '[YOUR DIRECTORY CONTAINING MELANOCYTE EXPRESS SUBDIRECTORIES HERE]' ## UPDATE THIS PATH
xprs_fnames2 <- Sys.glob(file.path(base_dir2, "*/results.xprs"))
xprs_fnames2

sample_names2 <- sub(file.path(base_dir2), "", xprs_fnames2) %>% sub("results.xprs", "", .) %>% gsub("/", "", .)
sample_names2

condition_names2 <- rep('melanocyte', length(sample_names2))

xprs2 <- read_express(xprs_fnames2, sample_names2, condition_names2)

ir2 <- newIntronRetention(xprs2$tpm, intron_to_trans, xprs2$condition, xprs2$uniq_counts)

n_filter2 <- 0.25
ir_filtered_tpm2 <- filter_low_tpm(ir2, 1, n_filter2) # Only include transcripts that are expressed at a tpm >= 1 in > n_filter % of samples
ir_filtered_psi2 <- filter_perfect_psi(ir_filtered_tpm2) # Remove transcripts in which every sample has 0 or 1 retention after rounding (artifact)
ir_filtered_frags2 <- filter_low_frags(ir_filtered_psi2, 5, n_filter2) # Only include transcripts that have at least 5 unique counts in > n_filter % of samples

write.table(ir_filtered_frags2$unique_counts, 
            file='[YOUR DIRECTORY HERE]/melanocyte_intron_retention_read_counts.tsv', ## UPDATE THIS PATH
            sep='\t', quote=F, col.names = T, row.names = F)

