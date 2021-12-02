#! /bin/sh
#METHOD 2 - align directly sequences to reference and get veariant calls into a database
# STEP 1: Extract sequences from large MSA, stip gaps and align to reference
#source activate seqtk 
#"/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/GISAID/29-06-2021/msa_0625/msa_0625.fasta"
samples2extract_file=$1
infasta_path=$2
reference_fasta=$3
#infasta_path="/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/MIA_mk/msa_0104/msa_0104.fasta"
#seqtk subseq $infasta_path seqtk_samples2extact.txt > extracted_fastas_subset.fasta
batch_number=$(echo $1 | grep -o "batch[0-9]\+")
out_name_suffix="_${batch_number}"

echo "grep -w -A 1 -Ff  $samples2extract_file  $infasta_path --no-group-separator > extracted_fastas_subset${out_name_suffix}.fasta"
grep -w -A 1 -Ff  $samples2extract_file  $infasta_path --no-group-separator > extracted_fastas_subset${out_name_suffix}.fasta

echo "$1 $2 $3"
echo "Working in $(pwd) directory"
out_folder="vcf_files"

mkdir -p ${out_folder}
samples_ids=($(cat ${samples2extract_file}))
echo ${#samples_ids[@]}


for sample_id in ${samples_ids[@]};do
	echo "${sample_id}" > sample2extact${out_name_suffix}.txt
	gisaid_accession=$(cat sample2extact${out_name_suffix}.txt | grep  -Eo "EPI_[[:alpha:]]+_[[:digit:]]+")
	grep -w -A 1 -Ff  sample2extact${out_name_suffix}.txt  extracted_fastas_subset${out_name_suffix}.fasta  > extracted_fasta${out_name_suffix}.fasta
	#seqtk subseq extracted_fastas_subset${out_name_suffix}.fasta sample2extact${out_name_suffix}.txt > extracted_fasta${out_name_suffix}.fasta
	
	sed -i 's/-//g' extracted_fasta${out_name_suffix}.fasta #remove alignment gap
	snippy_out_folder="snippy${out_name_suffix}"
    
    echo "snippy --ctgs extracted_fasta${out_name_suffix}.fasta --reference ${reference_fasta} --outdir $snippy_out_folder --force"
	snippy --ctgs extracted_fasta${out_name_suffix}.fasta --reference ${reference_fasta} --outdir $snippy_out_folder --force
	 
    #cat snippy/snps.vcf
    echo "COPY: cp ${snippy_out_folder}/snps.vcf ${out_folder}/${gisaid_accession}.vcf"
	cp ${snippy_out_folder}/snps.vcf ${out_folder}/${gisaid_accession}.vcf
done

#snpeff
#source activate snpeff@5.0.0
#data_directory="/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/MIA_mk/msa_0104/vcf_files/test/*"

#for vcf_file in $data_directory;do
#	echo "Processing $vcf_file sample"
#	snpEff NC_045512.2 -hgvs $vcf_file > tmp.vcf
#	mv tmp.vcf $vcf_file
	#break 
#done

#submit script batches of runs
#ls samples2extact_* | wc #number of samples to submit #number of samples to submit #372
#for i in `seq 1 153`;do 
#  echo $i; 
#  sed "s/batch1/batch$i/g" ./submit_scripts/submit_script.sh > ./submit_scripts/submit_script_${i}.sh 
#  sbatch --mem=4G -c 8 -p NMLResearch -J 0625_b${i} --time=184:00:00 --output=slurm_batch$i.out ./submit_scripts/submit_script_${i}.sh
#done;

#files completed
#grep "Using read file: /Drives/W/Projects" slurm_batch100.out | wc
#squeue -u kbessono | grep "PD" | wc
#squeue -u kbessono  -o "%.18i %.9P %.25j %.8u %.2t %.10M %.6D %R"

#sbatch --mem=4G -c 8 -p NMLResearch -J batch1 --time=120:00:00  ./submit_script.sh
# sbatch -p NMLReaserch --warp="tar -xf archive.tar.xz"
# sbatch -p NMLResearch -J tar --mem=32G -c 1 --wrap="tar -czf msa_20210825_vcfs.tar.gz gisaid_complete"
#gzip vcf files in the folder

# Then need to compress those vcf files to a single archive to prevent issues with filesystem
#diff -u /home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/29-06-2021/msa_0625/msa_0625_headers.txt /home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/msa_0822_headers.txt | grep '^\+' | sed -E 's/^\+//' > headers_diff.txt

#source activate R-3.6.1 & R
