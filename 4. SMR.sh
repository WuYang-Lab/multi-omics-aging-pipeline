#===========================================#
#                   smr                     #
#===========================================#

SMR="<YOUR_PATH>/smr-1.3.1-linux-x86_64/smr-1.3.1"
bfile="<YOUR_PATH>/g1000_eur"
eqtl="<YOUR_PATH>/eQTL/besd/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
GWASDIRT='<YOUR_PATH>/format/'
PHENO=("" "Phenotype1" "Phenotype2" "Phenotype3")
GWASDATA=("" "Phenotype1.txt" "Phenotype2.txt" "Phenotype3.txt")
NPHENO=$(($(echo ${#PHENO[@]}) - 1 ))
OUTPUT="<YOUR_PATH>/eqtl/"

for((p=1;p<=${NPHENO};p++))
do
${SMR} --beqtl-summary ${eqtl} --gwas-summary ${GWASDIRT}${GWASDATA[$p]} \
 --bfile ${bfile} --diff-freq-prop 0.5 --maf 0.01 --cis-wind 2000 --peqtl-smr 1e-5 \
 --out ${OUTPUT}${PHENO[$p]}_m2tsmr --thread-num 3 \

