# 1. Transform raw GWAS summary statistics into HDL input
library(data.table)
library(dplyr)
sumdir <- "<YOUR_INPUT_PATH>"
outdir <- "<YOUR_OUTPUT_PATH>"
files <- list.files(sumdir, pattern = "_format\\.txt$", full.names = TRUE)

# Process each file to compute Z-scores and prepare the data for HDL
process_file <- function(file_path) {
  data <- fread(file_path, header = TRUE, sep = "\t")
  data$Z <- with(data, b / se)  # Calculate Z-scores
  data <- data %>% rename(N = n)
  data <- data %>% select(SNP, A1, A2, N, b, se, Z)
file_name <- basename(file_path)
  output_name <- sub("_format\\.txt$", "_summary.txt", file_name)
  output_file <- file.path(outdir, output_name)
  fwrite(data, file = output_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
lapply(files, process_file)

# 2. Run HDL data wrangling for each GWAS summary file

# List of GWAS names (replace with actual file names)
gwas_names=("GWAS1" "GWAS2" "GWAS3")

# Loop over each GWAS file and run HDL data wrangling
for gwas_name in "${gwas_names[@]}"
do
    Rscript <YOUR_PATH>/HDL-master/HDL.data.wrangling.R \
    gwas.file="<YOUR_PATH>/${gwas_name}_summary.txt" \
    LD.path="<YOUR_PATH>/HDL-master/UKB_imputed_SVD_eigen99_extraction" \
    SNP=SNP A1=A1 A2=A2 N=N Z=Z \
    output.file="<YOUR_PATH>/result/${gwas_name}" \
    log.file="<YOUR_PATH>/result/${gwas_name}"
done


# 3. Estimating genetic correlation using HDL

# List of GWAS names for genetic correlation estimation
gwas_name=("GWAS1" "GWAS2" "GWAS3")
NGWAS=$((${#gwas_name[@]})) 
# Loop through each pair of GWAS traits and calculate genetic correlation using HDL
for ((o=0; o<${NGWAS}; o++))  
do
    for ((p=0; p<${o}; p++))
    do
            Rscript <YOUR_PATH>/HDL-master/HDL.run.R \
            gwas1.df=<YOUR_PATH>/result/${gwas_name[$o]}.hdl.rds \
            gwas2.df=<YOUR_PATH>/result/${trait_name[$p]}.hdl.rds \
            LD.path=<YOUR_PATH>/HDL-master/UKB_imputed_SVD_eigen99_extraction \
            output.file=<YOUR_PATH>/GC/${gwas_name[$o]}-${gwas_name[$p]}.Rout       
    done
done

# 4. Merge results from the output files

for file in <YOUR_PATH>/GC/*.Rout; do
    tail -n +24 "$file" >> hdl_result.txt
done




