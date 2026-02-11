# 1. GSMR
gcta-1.94.1 \
--bfile g1000_eur \
--gsmr-file risk_factors_path.txt aging_traits_path.txt \
--gsmr-direction 2 \
--gwas-thresh 5e-08 \
--clump-r2 0.05 \
--gsmr2-beta \
--heidi-thresh 0.01 \
--effect-plot \
--out risk_factors_aging_traits \


# 2. make scatter plot
source("gsmr_plot.r")
gsmr_data = read_gsmr_data("risk_factors_aging_traits.eff_plot.gz")
gsmr_summary(gsmr_data)      
pdf("/public/home/gaikai/risk/UKBB_data/format/gsmr/plot/risk_factors_aging_traits.pdf",height = 6, width =6)  
plot_gsmr_effect(gsmr_data, "risk_factors", "aging_traits", colors()[71]) 
dev.off()

