#===========================================#
#                   multi_MR                #
#===========================================#

# Set paths as variables (replace these placeholders with your actual paths)
gcta_path="<YOUR_PATH>/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"  # Path to the GCTA binary
bfile_path="<YOUR_PATH>/g1000_eur"  # Path to the reference data (e.g., 1000 Genomes EUR data)
gsmr_file_path_1="<YOUR_PATH>/risk_factor.txt"  # Path to the first GSMR file (IGF1)
gsmr_file_path_2="<YOUR_PATH>/aging_traits.txt"  # Path to the second GSMR file (Hannum)

# Output directory for results
output_dir="<YOUR_PATH>/risk_factors_aging_traits"  # Output directory where results will be saved

# Run GSMR (Genetic Summary Mendelian Randomization) analysis
${gcta_path} --bfile ${bfile_path} \
--gsmr-file ${gsmr_file_path_1} ${gsmr_file_path_2} \
--gsmr-direction 0 \
--gsmr-snp-min 10 \
--heidi-thresh 1e-700 \
--effect-plot \
--out ${output_dir} \

# Load necessary R scripts (adjust paths as needed)
source("<YOUR_PATH>/multi_MR/script/gsmr_plot_function.R")
source("<YOUR_PATH>/multi_MR/script/GSMR_effect_plot.R")
source("<YOUR_PATH>/multi_MR/script/comfun_file.R")
source("<YOUR_PATH>/multi_MR/script/MR.R")
library(survey)
library(NCmisc)
library(MASS)

setwd("<YOUR_PATH>/multi_MR/result/")
dir <- "<YOUR_PATH>/multi_MR/"
ex = c("risk_factor")
out = c("aging_traits")
res = c()

gsmr_data1 = read_gsmr_data(paste0(dir, ex, out, ".eff_plot.gz"))
resbuf = gsmr_snp_effect(gsmr_data1, ex, out)
bzx = resbuf$bzx
bzx_se = resbuf$bzx_se
bzx_pval = resbuf$bzx_pval
bzy = resbuf$bzy
bzy_se = resbuf$bzy_se
bzy_pval = resbuf$bzy_pval

# Loop through different methods for analysis
for(method in c("GSMR_heidi_v1", "GSMR_heidi_v3_stepwise", "GSMR_noheidi", "IVW", "Median", "Simple_Median", "Mode", "Egger", "Robust", "Lasso", "RAPS", "PRESSO")) {
    res_mr = all_mr(bzx = resbuf$bzx, bzx_se = resbuf$bzx_se, bzx_pval = resbuf$bzx_pval,
                    bzy = resbuf$bzy, bzy_se = resbuf$bzy_se, bzy_pval = resbuf$bzy_pval, method = method)
    res = rbind(res, data.frame(exposure = ex, outcome = out, method, p = res_mr$pval, bxy = res_mr$bxy, nr_pleio_rm = length(res_mr$pleio)))
}
print(paste(ex, "to", out, "completed", sep = " "))

# Calculate z-values and confidence intervals
res$z = p.to.Z(res$p)
res$se = res$bxy / res$z
res$z = res$z * sign(res$se)
res$se = res$se * sign(res$se)
res$association = paste(res$exposure, "_to_", res$outcome, sep = "")
res$lower = res$bxy - 1.96 * res$se
res$upper = res$bxy + 1.96 * res$se

# Use the GSMR results from CPP
gsmr_gene = read.table(paste0(dir, ex, "_Hannum.gsmr"), head = TRUE)
sigidx1 = match(ex, gsmr_gene$Exposure)
gsmr_gene_sig = gsmr_gene[sigidx1,]

res[res$method == "GSMR_heidi_v1", c("bxy", "se", "p")] = gsmr_gene_sig[, c(3:5)]
res$lower = res$bxy - 1.96 * res$se
res$upper = res$bxy + 1.96 * res$se

trait = gsmr_gene_sig$Exposure
res$association = rep(trait, each = length(table(res$method)))

# Modify method names
res[res$method == "GSMR_heidi_v1", "method"] = "GSMR"
res[res$method == "GSMR_heidi_v3_stepwise", "method"] = "GSMR_v2"
res[res$method == "Simple_Median", "method"] = "Simple\nMedian"
res[res$method == "Median", "method"] = "Weighted\nMedian"
res = res[!res$method %in% c("Lasso", "GSMR_noheidi", "GSMR_v2", "MRMix"),]

# Build output file path and write the results
output_file = file.path(getwd(), paste0(ex, "_MR_comparison_results_updated.txt"))
write.table(res, file = output_file, row.names = FALSE, col.names = TRUE, quote = TRUE, sep = "\t")

#MRlap

library(MRlap)
risk_factors <- "<your_path>/risk_factors.txt"
aging_traits <- "<your_path>/aging_traits"

# Running the MRlap analysis
A = MRlap(exposure = risk_factors,
          exposure_name = "risk_factors",
          outcome = aging_traits,
          outcome_name = "aging_traits",
          ld = "<your_path>/eur_w_ld_chr",
          hm3 = "<your_path>/w_hm3.snplist")


# 2.make plot
setwd("<YOUR_PATH>/multi_MR/")
files <- list.files(pattern = "_MR_comparison_results_updated.txt")

res <- data.frame()
for (file in files) {
  data <- read.table(file, header = TRUE) 
  res <- rbind(res, data)
}

res$method = as.factor(res$method)
res$association = as.factor(res$association)
res$association = factor(res$association, levels = levels(res$association))

# Create the plot
library(data.table)
library(ggplot2)

res$sig <- ifelse(res$p < 6.41e-4, "**", ifelse(res$p < 0.05, "*", ""))
res <- res[order(res$group, res$association), ]
res$association <- factor(res$association, levels = unique(res$association))

# Custom theme for the plot
custom_theme <- theme(
  legend.title = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "none",
  strip.text = element_text(face = "bold", size = 10),
  axis.text.x = element_text(face = "bold", size = 10),
  axis.text.y = element_text(face = "bold", size = 10),
  axis.title.y = element_text(face = "bold", size = 12),
  strip.background = element_blank()
)

# Create the plot
p1 <- ggplot(data = res, aes(y = bxy, x = method, colour = method)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = group), alpha = 0.2, color = NA) +
  geom_point(aes(fill = method), shape = 16) +  # Solid points
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.5) +  # Error bars
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(rows = vars(association), switch = "y") +
  custom_theme +
  xlab("") +
  labs(y = expression(hat(italic(b))[xy])) +
  geom_text(aes(label = sig), vjust = -1)

# Print the plot
print(p1)

# Save the plot as a PDF
ggsave("<YOUR_PATH>/multi_MR/risk_factor_MR_comparison_results.pdf", p1, width = 9, height = 10)
