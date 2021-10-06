#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#GISAIDmetadatafilein="/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/metadata.tsv"
#GISAIDmetadata=read.csv2(file=GISAIDmetadatafilein, check.names = F, sep="\t", as.is = T, stringsAsFactors = F)

#Accession2LineageList=as.list(GISAIDmetadata[,"Pango lineage"]); names(Accession2LineageList) = GISAIDmetadata[,"Accession ID"]
#lapply(GISAIDmetadata[1:10,'Accession ID'],function(idx){idx[,"Pango lineage"]})
#save(Accession2LineageList,file="/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/vcf_files/Accession2LineageList.R")

#query_path="/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/vcf_files/gisaid_complete/"
#vcf_files_list = list.files(query_path)
#save(vcf_files_list,file='/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/vcf_files/vcf_files_list.R')




#if (length(args)==0) {
#  stop("At least one argument must be supplied (range number)", call.=FALSE)
#} else {
#  range_num = as.numeric(args[1]) #the range number out of so many to work on
#}

#query_path = as.character(args[1])


#if ("vcf_files_list" %in% ls() == FALSE){
#vcf_files_list = list.files(query_path)
cat(args[1],"\n")
vcf_files_list = read.table(file=args[1])[,1]

#}
ranges_idx = c(seq(0,length(vcf_files_list),by=4999),length(vcf_files_list))
#ranges_list=list()
for(i in 1:length(ranges_idx)){
  if(i==length(ranges_idx)){break}
  cat(ranges_idx[i]+1, ranges_idx[i+1],"\n")
  #print(length(as.data.frame(vcf_files_list[(ranges_idx[i]+1):(ranges_idx[i+1])])[,1]))
  write.table(as.data.frame(vcf_files_list[(ranges_idx[i]+1):(ranges_idx[i+1])])[,1],
             file=paste("vcf_file_names_batch_",i,".txt", sep="",collapse=""),
             col.names=F,row.names = F, quote=F
            )
  #ranges_list[[i]]=c(ranges_idx[i]+1,ranges_idx[i+1])
}
