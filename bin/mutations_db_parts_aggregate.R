#!/usr/bin/env Rscript
require(jsonlite) 
require(stringr)

args = commandArgs(trailingOnly=TRUE)
json_files_list=list()
for (file in args){
    if( file != "null" ){
      json_files_list=append(json_files_list,file)
    }
}

#json_files_list = list.files(query_path)[str_detect(list.files(query_path),"part\\d+\\.json")]

mutations_database=list()
for (json_part_file in json_files_list){
  cat(json_part_file,"\n")
  part_data = jsonlite::fromJSON(json_part_file)
  mutations_database = append(mutations_database,part_data)
}

length(mutations_database)
save(mutations_database, file="mutations_database_final.Rdata")


write(jsonlite::toJSON(mutations_database, pretty=T), file = "mutations_database_final.json")


#sort(table(unlist(lapply(mutations_database,function(x){x$lineage}))),decreasing = T)[1:10]