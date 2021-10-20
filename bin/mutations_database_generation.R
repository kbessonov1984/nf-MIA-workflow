#!/usr/bin/env Rscript
require(jsonlite)
require(vcfR)
require(stringr)
args = commandArgs(trailingOnly=TRUE)




#lapply(GISAIDmetadata[1:10,'Accession ID'],function(idx){idx[,"Pango lineage"]})
#save(Accession2LineageList,file="/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/vcf_files/Accession2LineageList.R")

#query_path="/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/vcf_files/gisaid_complete/"
#vcf_files_list = list.files(query_path)
#save(vcf_files_list,file='/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/vcf_files/vcf_files_list.R')


if (length(args)==0) {
  stop("At least one argument must be supplied (range number)", call.=FALSE)
} else {
  vcf_files_list  = read.table(file=args[1])[,1] #the range number out of so many to work on
}

cat("VCF files to process",length(vcf_files_list),"\n")
GISAIDmetadatafilein=args[2]
cat("Reading GISAID metadata file ... \n")
GISAIDmetadata=read.csv2(file=GISAIDmetadatafilein, check.names = F, sep="\t", as.is = T, stringsAsFactors = F)
cat("Generating Accession2LineageList ... \n")
Accession2LineageList=as.list(GISAIDmetadata[,"Pango lineage"])
names(Accession2LineageList) = GISAIDmetadata[,"Accession ID"]


batch_n=stringr::str_extract(args[1],"batch_\\d+")
cat("Found ",length(vcf_files_list)," vcf files\n")



mutations_database=list(); num_files_processed=0; 
total_files=length(vcf_files_list)

for (filepath in vcf_files_list){
  filename=basename(filepath)  
  num_files_processed=num_files_processed+1
  accession = str_remove(filename,"\\.vcf")
  
  metadata_lineage=Accession2LineageList[[accession]]
  #metadata_lineage = GISAIDmetadata_subset[GISAIDmetadata_subset["Accession ID"] == accession,"Pango lineage"]
  if(length(metadata_lineage) == 0){
    metadata_lineage="NA"
  }
  
  #progress bar
  if(num_files_processed%%1000 == 0){ cat("Processed ",num_files_processed/total_files*100,"% on ",
                                          as.character(Sys.time()),"\n")}
  
  #cat("Reading ",filepath,"...\n")
  vcf_data = read.vcfR(filepath,verbose = FALSE)
  
  
  if(dim(vcf_data@fix)[1] == 0){
    mutations_database[[accession]] = list(lineage = unbox(metadata_lineage), mutations = NULL)
    cat("WARNING: Empty file",filename,"\n")
    next; #empty vcf file
  }
  
  mutations_vector=c()
  for(row_idx in 1:dim(vcf_data@fix)[1]) {
    mut_name = paste(vcf_data@fix[row_idx,c("REF","POS","ALT")],collapse="")
    mutations_vector = c(mutations_vector,mut_name)
  }
  
  mutations_database[[accession]] = list(lineage = unbox(metadata_lineage), 
                                         mutations = mutations_vector
  )
  
  
}

cat("Numeber of accessions added to the database ",length(names(mutations_database)),"\n")
## TODO need database on mutation effects for results mapping in future (e.g. A12162G mutation effect and ORF)

db_json_name=paste("mutataions_db_",batch_n,".json", sep="",collapse = "")
write(jsonlite::toJSON(mutations_database, pretty=T), file = db_json_name)

# submit to SLURM
# sbatch --mem=32G -c 1 -p NMLResearch -J 0625_R --time=184:00:00 --output=slurm_R.out --wrap="source activate R-3.6.1 && Rscript mutations_database_generation.R"
# for i in `seq 1 77`;do  sbatch --exclude="waffles-g-9,waffles-g-7" --mem=8G -c 1 -p NMLResearch -J R_b$i --time=24:00:00 --output=slurm_R_b$i.out --wrap="source activate R-3.6.1 && Rscript --vanilla mutations_database_generation.R $i";done


