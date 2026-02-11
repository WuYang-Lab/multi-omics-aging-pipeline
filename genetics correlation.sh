#===============================#
#       LDSC (Linkage Disequilibrium Score Regression) Calculation
#===============================#

# Activate the environment for LDSC
conda activate <your_env>

# ============================ #
# 1. Convert the summary statistics
# ============================ #

# Path to LDSC munge_sumstats.py script
ldsc="ldsc-master/munge_sumstats.py"

# Phenotype file names
PHENO=(
  "phenotype1_format"
  "phenotype2_format"
  "phenotype3_format"
  # Add more phenotype files here as needed
)

NPHENO=$(($(echo ${#PHENO[@]}) - 1))
OUTPUT="output/ldsc_GC/"

# Loop through each pair of GWAS traits and calculate genetic correlation using LDSC
for ((o=0; o<${NGWAS}; o++))
do
    for ((p=0; p<${o}; p++))  # Limiting p to less than o
    do
       # Run LDSC directly without variable concatenation
       ldsc --rg ${gwas_name[$o]}.sumstats.gz,${trait_name[$p]}.sumstats.gz \
       --ref-ld-chr ldsc-master/LDSC/eur_w_ld_chr/ \
       --w-ld-chr ldsc-master/LDSC/eur_w_ld_chr/ \
       --out ${OUTPUT}${gwas_name[$o]}_${trait_name[$p]}
    done
done

# ============================ #
# 2. Genetic correlation calculation
# ============================ #

# Path to LDSC ldsc.py script
ldsc="ldsc-master/ldsc.py"

# GWAS summary file names for Aging related traits
gwas_aging=(
  "summary1"
  "summary2"
)

gwas_risk_factors=(
  "trait1"
  "trait2"
)

NGWAS_AGING=$((${#gwas_aging[@]}))
NGWAS_RISK_FACTORS=$((${#gwas_risk_factors[@]}))
OUTPUT="output/ldsc_GC/"

# Loop through each pair of GWAS traits and calculate genetic correlation using LDSC
for ((o=0; o<${NGWAS_AGING}; o++))
do
    for ((p=0; p<${NGWAS_RISK_FACTORS}; p++))
    do
       # Run LDSC command
       ldsc --rg ${gwas_aging[$o]}.sumstats.gz,${gwas_risk_factors[$p]}.sumstats.gz \
       --ref-ld-chr ldsc-master/LDSC/eur_w_ld_chr/ \
       --w-ld-chr ldsc-master/LDSC/eur_w_ld_chr/ \
       --out ${OUTPUT}${gwas_aging[$o]}_${gwas_risk_factors[$p]}
    done
done

# ============================ #
# 3. Merge the results
# ============================ #

# Extract the 62nd line from all log files to summarize the results
for file in *.log; do
    sed -n '62p' "$file"
done > ldsc.log
