#!/usr/bin/env Rscript
require(jsonlite)

args = commandArgs(trailingOnly=TRUE)

GISAIDmetadatafilein=args[1]
json_db_file = args[2]
db_json_name="mutations_database_homogenized_final.json"

cat("Reading GISAID metadata file ... \n")
GISAIDmetadata=read.csv2(file=GISAIDmetadatafilein, check.names = F, sep="\t", as.is = T, stringsAsFactors = F)
rownames(GISAIDmetadata) = GISAIDmetadata[,"Accession ID"]
mutations_database = jsonlite::fromJSON(json_db_file)

cat("Homogenize Pango Lineages in the database and GISAID metadata as this data changes from one release to other ... \n")
GISAIDmetadata_with_sequences = GISAIDmetadata[GISAIDmetadata[,"Accession ID"] %in% names(mutations_database),] #homogenize metadata with mutations_database

mutations_database = mutations_database[GISAIDmetadata_with_sequences[,"Accession ID"]]
mutations_database=lapply(mutations_database, function(i,l){i$lineage=l[parent.frame()$i]; return(i)},l=
                            GISAIDmetadata_with_sequences[names(mutations_database),"Pango lineage"])

write(jsonlite::toJSON(mutations_database, pretty=T), file = db_json_name)
cat("Wrote final database as ",db_json_name,"\n")
save(mutations_database, file="mutations_database_homogenized_final.Rdata")