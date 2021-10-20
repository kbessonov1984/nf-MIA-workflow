#given metadata, output the FASTA headers of samples to extract with grep for further processing

require(stringr)
require(jsonlite)

args = commandArgs(trailingOnly=TRUE) #get range number

if (length(args)==0) {
  stop("At least one argument must be supplied (e.g. path to GISAID metadata file)", call.=FALSE)
} else {
  range_num = as.numeric(args[1]) #the range number out of so many to work on
}

GISAIDmetadatafilein=args[1]
mutations_database_path = args[2]

#"/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/MIA_mk/msa_0104/GISAID_msa014_metadata.tsv"
#"/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/metadata.tsv"
#"/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/GISAID/29-06-2021/metadata.tsv"
#"/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/metadata.tsv"

#GISAIDmetadata=read.table(file=GISAIDmetadatafilein, check.names = F, sep="\t", as.is = T, stringsAsFactors = F)
cat("Openning metadata file ",GISAIDmetadatafilein,"\n")
GISAIDmetadata=read.csv2(file=GISAIDmetadatafilein, check.names = F, sep="\t", as.is = T, stringsAsFactors = F)
#idx_ingroup= which(GISAIDmetadata[,"Pango lineage"] == "B.1.617.2") 
#length(idx_outgroup)
#"/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/MIA_mk/msa_0104/msa_0104_headers.txt"
#"/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/GISAID/29-06-2021/msa_0625/msa_0625_headers.txt"

#load("/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/29-06-2021/vcf_files/session.R")


#write.table(GISAIDmetadata[GISAIDmetadata[,"Accession ID"] %in% names(mutations_database),"Virus name"],quote=F,
#            file="sample_headers.txt",
#            col.names=F,row.names = F)
#"/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/GISAID/29-06-2021/msa_0625/msa_0625_headers.txt"

#headers_fasta = read.csv2(file="/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/msa_0822_headers.txt",
#                           stringsAsFactors = F, as.is = T,colClasses = "character")[,1]
headers_fasta=paste(GISAIDmetadata[,"Virus name"],"|",GISAIDmetadata[,"Accession ID"],sep="")
#accessions=unlist(str_extract_all(headers_fasta[1:10],"EPI_\\w+_\\d+"))

#If database is provided, only missing samples will be analyzed for DB update
#Otherwise ALL samples will be used in analysis
if(is.na(mutations_database_path) == FALSE){
    cat("Identifying missing accession numbers from the current mutations database\n")
    mutations_database=jsonlite::fromJSON(mutations_database_path)
    accessions=GISAIDmetadata[,"Accession ID"]
    idx_not_analyzed = which(!accessions %in% names(mutations_database))
    headers_fasta=headers_fasta[idx_not_analyzed]
}

#accessions = unlist(lapply(str_split(headers_fasta,"\\|"),function(x){x[2]}))

#headers_fasta=headers_fasta[nchar(headers_fasta)>20]
#idx = which(unlist(str_match_all(headers_fasta,"EPI_\\w+_\\d+")) %in% GISAIDmetadata[idx_ingroup, "Accession ID"])
#length(idx)

# generation of ranges
ranges_idx = c(seq(0, length(headers_fasta),by=4999), length(headers_fasta)) #vector of range indices
ranges_list=list()
for(i in 1:length(ranges_idx)){
    if(i==length(ranges_idx)){break}
      cat(ranges_idx[i]+1,ranges_idx[i+1],"\n")
      ranges_list[[i]]=c(ranges_idx[i]+1,ranges_idx[i+1])
}

length(headers_fasta)
cat("Ranges ", length(ranges_list),"\n")

for (batch_n in 1:length(ranges_list)){
    headers_out = headers_fasta[ranges_list[[batch_n]][1]:ranges_list[[batch_n]][2]]

    #hearders_out = headers_fasta[idx] #if extracting specific lineage
    #"/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/MIA_mk/msa_0104"
    #"/Drives/W/Projects/Project_Enterics_Wastewater/Shared_Data/GISAID/29-06-2021/samples2extact_batch"
    write.table(headers_out , file=paste("samples2extact_batch",batch_n,".txt",sep=""), 
                quote=F, col.names = F, row.names = F)
}
