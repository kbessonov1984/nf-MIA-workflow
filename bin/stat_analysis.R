#!/usr/bin/env Rscript

library(stringr)
library(jsonlite)
if(require("islasso", quietly = TRUE, warn.conflicts = F)==FALSE){
  install.packages("islasso",repos='http://cran.us.r-project.org')
  library("islasso")
}

if(require("corrplot", quietly = TRUE, warn.conflicts = F)==FALSE){
  install.packages("corrplot",repos='http://cran.us.r-project.org')
  library("corrplot")
}

obtainLASSO_padjust_maxT_values <- function(data2fit, is_lasso_summary_df, NPERMUTATIONS){
  #mix original dataframe
  data2fit_permut=data2fit
  data2fit_permut[,"Group"]=sample(data2fit[,"Group"],length(data2fit[,"Group"]))
  
  ##PERMUTATION VALUES
  permut_Z_matrix=matrix(0,ncol=NPERMUTATIONS,nrow=length(rownames(is_lasso_summary_df)))
  rownames(permut_Z_matrix)=rownames(is_lasso_summary_df)
  for(p in 1:NPERMUTATIONS){
    cat("Permutation #",p,"\n")  
    data2fit_permut=data2fit
    data2fit_permut[,"Group"]=sample(data2fit[,"Group"],length(data2fit[,"Group"]))
    
    islasso_fit_permut = islasso(Group~0+.,data = as.data.frame(data2fit_permut), family=binomial(link = "probit"))
    islasso_fit_permut_zvals = summary(islasso_fit_permut)$coefficients[,"z value"]
    permut_Z_matrix[names(islasso_fit_permut_zvals),p]= islasso_fit_permut_zvals
    
  }
  
  
  ##STEP-DOWN MAXT p-adjusted values (strong FWER adjustment)
  is_lasso_summary_df[,"p_adj(maxT)"]=NA
  for(mutation in rownames(is_lasso_summary_df)){
    z_obs = is_lasso_summary_df[mutation,"z value"]
    b_gt_eq_tobs = length(which(sapply(1:NPERMUTATIONS,function(i){any(abs(permut_Z_matrix[,i]) >= abs(z_obs))})))  
    if(b_gt_eq_tobs == 0){b_gt_eq_tobs=1}
    p_adj = b_gt_eq_tobs/NPERMUTATIONS
    #cat(mutation," raw p-val",is_lasso_summary_df[mutation,"Pr(>|z|)"],"adjusted p-val",p_adj,"\n" )
    is_lasso_summary_df[mutation,"p_adj(maxT)"] = p_adj
  }
  
  return(list(df=is_lasso_summary_df, permut_Z_matrix=permut_Z_matrix))
}


args = commandArgs(trailingOnly=TRUE)
print(args)
GISAIDmetadatafilein=args[1] # "/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/metadata.tsv"
GISAIDdatabase_path=args[2] #"/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/vcf_files/mutations_database_final.json"
INLINEAGE = args[3] #"B.1.617.2"
NPERMUTATIONS = as.numeric(args[4]) # 100
QUERY_LOCATION = args[5] # CANADA
START_DATE = args[6]
END_DATE = args[7]


mutations_database = jsonlite::fromJSON(GISAIDdatabase_path)

GISAIDmetadata=read.csv2(file=GISAIDmetadatafilein, check.names = F, sep="\t", as.is = T, stringsAsFactors = F)
GISAIDmetadata[,"Submission date"]=as.Date(GISAIDmetadata[,"Submission date"])
rownames(GISAIDmetadata)=GISAIDmetadata[,"Accession ID"]
GISAIDmetadata_with_sequences = GISAIDmetadata[GISAIDmetadata[,"Accession ID"] %in% names(mutations_database),] #homogenize metadata with mutations_database

#filter by geographical location
cat(QUERY_LOCATION,"\n");
if( QUERY_LOCATION != "null" ){
  GISAIDmetadata_with_sequences = GISAIDmetadata_with_sequences[stringr::str_detect(
    GISAIDmetadata_with_sequences[,"Location"],QUERY_LOCATION),]
}

if( is.na(START_DATE) == FALSE & is.na(END_DATE) == FALSE ){
  idx_select=GISAIDmetadata_with_sequences[,"Submission date"] >= as.Date(START_DATE) & 
    GISAIDmetadata_with_sequences[,"Submission date"] <= as.Date(END_DATE)
  GISAIDmetadata_with_sequences=GISAIDmetadata_with_sequences[idx_select,]
}else if(is.na(START_DATE) == FALSE){
  idx_select=GISAIDmetadata_with_sequences[,"Submission date"] >= as.Date(START_DATE)
  GISAIDmetadata_with_sequences=GISAIDmetadata_with_sequences[idx_select,]
}

cat("GISAID metadata dimensions after filtering (if any):\n")
print(dim(GISAIDmetadata_with_sequences))

if(dim(GISAIDmetadata_with_sequences)[1] == 0){
  stop("The GISAID metadata dataframe is empty probably due to filtering by date or location settings\n")
}

# filter mutations database and select accessions representing IN and OUT groups
mutations_database_filtered = mutations_database[GISAIDmetadata_with_sequences[,"Accession ID"]]
idx_not_na = which(unlist(lapply(mutations_database_filtered,function(i){!is.na(i$lineage)})))
in_group_accessions = names(mutations_database_filtered[idx_not_na])[unlist(lapply(mutations_database_filtered[idx_not_na],function(i){i$lineage == INLINEAGE}))]
out_group_accessions = names(mutations_database_filtered[idx_not_na])[unlist(lapply(mutations_database_filtered[idx_not_na],function(i){i$lineage != INLINEAGE}))]

#sample 10% of the IN-GROUP if its size is > 10K of samples
n_samples=length(in_group_accessions)
if(n_samples == 0){stop("No samples were found for the following run paramenters. Relax parameters or update database?")}
if(n_samples > 10000){n_samples=round(n_samples*0.10);cat("Total chunks to run",round(n_samples/(n_samples*0.10))," with each chunk",n_samples,"\n")}

# SAMPLE STRATEGY 1: Random subsampling of the IN-GROUP and OUT-GROUP
in_group_accessions_subset  = sample(in_group_accessions,n_samples)
out_group_accessions_subset = sample(out_group_accessions,n_samples) #random sampling
cat("Subset in-group samples",length(in_group_accessions_subset),"\n")

#common mutations in the in-group found at least in 10% of the samples
ingroup_mutations=unlist(sapply(mutations_database_filtered[in_group_accessions_subset],
                                function(sample_id){sample_id[["mutations"]]}))

ingroup_mutations_freq = table(ingroup_mutations) 
ingroup_mutations_freq_select = ingroup_mutations_freq[ingroup_mutations_freq  > round(length(in_group_accessions_subset)*0.1)]
ingroup_mutations = names(ingroup_mutations_freq_select)
cat("Ingroup mutations",length(ingroup_mutations),"\n")


# populate data matrix for model fitting
data2fit=matrix(0,nrow=length(in_group_accessions_subset)+length(out_group_accessions_subset),
                ncol=length(ingroup_mutations)+1, #+1 for the response variable (group membership)
                dimnames = list(c(in_group_accessions_subset,out_group_accessions_subset),
                                c("Group",ingroup_mutations))
)

data2fit[in_group_accessions_subset,"Group"]=1 
data2fit[out_group_accessions_subset,"Group"]=0

data2fit_coln=length(colnames(data2fit))
data2fit_colnames=colnames(data2fit)

loop_counter=1
for(sample_accession in mutations_database_filtered[c(in_group_accessions_subset,out_group_accessions_subset)]){
  mutations_selected_idx = which(data2fit_colnames %in% sample_accession$mutations)
  data2fit[loop_counter,mutations_selected_idx] = 1
  loop_counter=loop_counter+1
}

# Fit Linear Probability Model using OLS
single_mutation_stats_df=data.frame()
r_adj_model_values=c()
for (mutation in colnames(data2fit)[-1]){
  cat("Mutation",mutation,"\n")
  lmfit = glm(Group~0+.,data=as.data.frame(data2fit[,c("Group",mutation)]), family=binomial(link = "logit"))
  summary_df = summary(lmfit)$coefficients
  r_adj_model_values=c(r_adj_model_values, summary(lmfit)[["adj.r.squared"]])
  single_mutation_stats_df=rbind(single_mutation_stats_df, summary_df)
}
single_mutation_stats_df$`r_adj_model_values`=r_adj_model_values
single_mutation_stats_df = single_mutation_stats_df[sort(single_mutation_stats_df[,"z value"],decreasing = T,index.return=T)$ix,]
single_mutation_stats_df[sort(single_mutation_stats_df[,"z value"],decreasing = T,index.return=T)$ix,]

#top_sig_mutations_selected_names = rownames(single_mutation_stats_df)[which(single_mutation_stats_df[,"Pr(>|t|)"]<0.001)]
#length(top_sig_mutations_selected_names)

##PERMUTATION VALUES
permut_T_matrix=matrix(0,ncol=NPERMUTATIONS,nrow=length(colnames(data2fit)[-1]))
rownames(permut_T_matrix)=colnames(data2fit)[-1]
for(p in 1:NPERMUTATIONS){
  cat("Permutation #",p,"\n")  
  data2fit_permut=data2fit
  data2fit_permut[,"Group"]=sample(data2fit[,"Group"],length(data2fit[,"Group"]))
  for (mutation in colnames(data2fit_permut)[-1]){
    lmfit = glm(Group~0+.,data=as.data.frame(data2fit_permut[,c("Group",mutation)]),family=binomial(link = "logit"))
    permut_T_matrix[mutation,p]= summary(lmfit)$coefficients[,"z value"]
  }
}

##STEP-DOWN MAXT p-adjusted values (strong FWER adjustment)
single_mutation_stats_df[,"p_adj(maxT)"]=NA
for(mutation in rownames(single_mutation_stats_df)){
  t_obs = single_mutation_stats_df[mutation,"z value"]
  b_gt_eq_tobs = length(which(sapply(1:NPERMUTATIONS,function(i){any(abs(permut_T_matrix[,i]) >= abs(t_obs))})))  
  if(b_gt_eq_tobs == 0){b_gt_eq_tobs=1}
  p_adj = b_gt_eq_tobs/NPERMUTATIONS
  cat(mutation," raw p-val",single_mutation_stats_df[mutation,"Pr(>|z|)"],"adjusted p-val",p_adj,"\n" )
  single_mutation_stats_df[mutation,"p_adj(maxT)"] = p_adj
}

pdf(file="Pearson_correlation_plot_of_all_in-group_mutations.pdf",height=10,width=10)
corrplot(cor(data2fit[,-1], method = "pearson") ,type="upper", 
         title=paste("Correlation plot of all", length(colnames(data2fit[,-1])) ,"variables used in stats models"), 
         tl.cex = 0.8,cex.main = 2, mar=c(0,0,2,0))
dev.off()


## MULTIVARIATE ANALYSIS: Penalized regression
# First round all coevicients
islasso_fit = islasso(Group ~ 0+., data = as.data.frame(data2fit), family=binomial(link = "probit"))
#islasso_fit = islasso(Group ~ 0+., data = as.data.frame(data2fit[,c("Group",names(sort(islasso_fit$coefficients, decreasing=T))[1:30])]) ) 
is_lasso_summary_df = as.data.frame(summary(islasso_fit)$coefficients)
is_lasso_summary_df = is_lasso_summary_df[sort(is_lasso_summary_df[,"Estimate"],decreasing=T,index.return=T)$ix,]
#sort_beta_idx_gt_0 = which(is_lasso_summary_df["Estimate"] > 0)
#top_betas_idx = which(is_lasso_summary_df["Estimate"] >= quantile(is_lasso_summary_df[,"Estimate"],.75))

#2nd round: SELECT only mutations that have >0 betas and re-fit a reduced model
#selected_mutations = rownames(is_lasso_summary_df[is_lasso_summary_df[,"Estimate"] > 0,])
#selected_mutations = rownames(is_lasso_summary_df[top_betas_idx,])
#islasso_fit = islasso(Group ~ 0+., data = as.data.frame(data2fit[,c("Group",selected_mutations)]), 
#                      family=binomial(link = "probit"))
#is_lasso_summary_df = as.data.frame(summary(islasso_fit)$coefficients)
#only consider beta-coefficients that are > 0 (positive correlation with the group)


is_lasso_permut_results = obtainLASSO_padjust_maxT_values(data2fit,
                                                          is_lasso_summary_df, NPERMUTATIONS)
is_lasso_summary_df = is_lasso_permut_results$df
permut_Z_matrix = is_lasso_permut_results$permut_Z_matrix

#3rd round only working with mutations that do not display "strong" signal in permutations results (below 90% quantile)
# mean_permut_z_vals_vector=apply(abs(permut_Z_matrix),1,mean)
# selected_mutations = names(mean_permut_z_vals_vector[mean_permut_z_vals_vector < quantile(mean_permut_z_vals_vector,0.75)])
# islasso_fit = islasso(Group ~ 0+., data = as.data.frame(data2fit[,c("Group",selected_mutations)]), family=binomial(link = "logit"))
# is_lasso_permut_results = obtainLASSO_padjust_maxT_values(data2fit[,c("Group",selected_mutations)],
#                                                           as.data.frame(summary(islasso_fit)$coefficients),
#                                                           NPERMUTATIONS)
# is_lasso_summary_df = is_lasso_permut_results$df




## REPORTS ##
### OVERALL REPORT ###
in_group_report_name=paste("INGROUP_",INLINEAGE,"_OVERALL_MAIN_report.txt",sep="") 
cat("IN-group accessions overall after metadata filtering (if at all) #:",length(in_group_accessions),"\n",
    file = in_group_report_name, sep="\t",append = FALSE)
cat("IN-group PANGO Lineages distribution (should be a single PANGO lineage):\n", file = in_group_report_name, 
    append = TRUE) #QC check
write.table(as.data.frame(table(unlist(lapply(mutations_database_filtered[in_group_accessions],
                                              function(item){item$lineage})))),
            file = in_group_report_name, quote=F, append = TRUE, sep="\t",row.names = F, col.names = F)
cat("Sampling strategy: Random\n",file = in_group_report_name, sep="\t",append = TRUE)
cat("Number of permutations for maxT p-value calc:",NPERMUTATIONS,"\n",file = in_group_report_name, sep="\t",append = TRUE)
cat("Metadata location filter:",QUERY_LOCATION,"\n",file = in_group_report_name, sep="\t",append = TRUE)
cat("Metadata submission date filter start date:",START_DATE,"\n",file = in_group_report_name, sep="\t",append = TRUE)
cat("Metadata submission date filter end date:",END_DATE,"\n",file = in_group_report_name, sep="\t",append = TRUE)



cat("IN-group subsampling dataset size:",length(in_group_accessions_subset)," out of ",
    length(in_group_accessions)," samples\n", sep="\t",file = in_group_report_name, append = T)
cat("OUT-group subsampling dataset size:",length(out_group_accessions_subset),"out of ",
    length(out_group_accessions),"samples\n",sep="\t", file = in_group_report_name, append = T)
cat("IN-group mutations total in subsample:",length(ingroup_mutations),"\n",file = in_group_report_name, sep="\t",append=T)
write.table(as.data.frame(sort(ingroup_mutations_freq_select, decreasing = T)),
            file = in_group_report_name, sep="\t",quote=F, append = TRUE, row.names = F, col.names = T)

cat("\nOUT-group lineages distribution in subsample:\n",file = in_group_report_name, sep="\t",append = T)
write.table(as.data.frame(sort(table(unlist(lapply(mutations_database_filtered[out_group_accessions_subset],function(item){item$lineage}))),decreasing=T)),
            file = in_group_report_name, sep="\t",quote=F, append = TRUE, row.names = F, col.names = F)


cat("Univariate analysis of ",INLINEAGE," on ",length(in_group_accessions_subset)
    +length(out_group_accessions_subset)," samples based on ",NPERMUTATIONS," permutations\n",
    file="INGROUP_univariate_single_mutations_stats.txt",sep="")
write.table(single_mutation_stats_df, file="INGROUP_univariate_single_mutations_stats.txt", 
            quote=F, sep="\t", col.names = T, append=T)

name_lasso_report_file = "INGROUP_multivariate_LASSO_mutations_constalations_stats.txt"
cat("Multivariate LASSO analysis of ",INLINEAGE," on ",length(in_group_accessions_subset)
    +length(out_group_accessions_subset)," samples\n",file=name_lasso_report_file,sep="")
cat(paste("Final model: Group ~ ",paste(rownames(is_lasso_summary_df),collapse="+"),"\n",sep=""),
    file=name_lasso_report_file,sep="", append = T)
cat(paste("Number of variables in this model: ",length(row.names(is_lasso_summary_df)),"\n",sep=""),
    file=name_lasso_report_file,sep="", append = T)
cat(paste("Number of permutations: ",NPERMUTATIONS,"\n",sep=""),
    file=name_lasso_report_file,sep="", append = T)
write.table(is_lasso_summary_df, file=name_lasso_report_file, 
            quote=F, sep="\t", col.names = T, append=T)



## PLOTS ##
GISAID_lineages_freq=table(GISAIDmetadata[names(mutations_database),"Pango lineage"],useNA="no")
GISAID_lineages_freq_perc = GISAID_lineages_freq/sum(GISAID_lineages_freq) #percentages of total GISAID dataset
GISAID_lineages_freq_perc[GISAID_lineages_freq_perc>1]=1 #nomralization in case the counts are about totals

pdf(file="Top100_lineages_distribution_abs_freq.pdf",height=4,width=12)
plot(sort(table(GISAIDmetadata_with_sequences[,"Pango lineage"]),decreasing = T)[1:100], 
    main=paste("Top 100 PANGO lineages distribution overall based on ",dim(GISAIDmetadata_with_sequences)[1]," samples",
               sep=""), las=2,cex.axis=0.5,ylab="Frequency")
dev.off()


pdf(file="Top100_PANGO_LINEAGES_distribution_GISAID_percent_cumulative.pdf",height=4, width=12)
plot(sort(GISAID_lineages_freq_perc, decreasing = T)[1:100],
     main=paste("Top 100 PANGO Lineages distribution in the GISAID (",dim(GISAIDmetadata_with_sequences)[1]," samples)",sep=""),las=2,cex.axis=0.5,ylab="Percent total (%)",xlab="")
dev.off()

pdf(file="Top100_OUT_group_lineages_distribution.pdf",width=12,height = 5)
plot(sort(table(GISAIDmetadata[out_group_accessions_subset,"Pango lineage"]),decreasing = T)[1:100],
     main=paste("Top 100 PANGO lineages distribution in the OUT-group subset based on ",
                length(out_group_accessions_subset)," samples\n",sep=""),las=2,cex.axis=0.5,ylab="Frequency" )
dev.off()

pdf(file="INGROUP_mutations_freq_distribution.pdf",height=4, width=10)
plot(sort(ingroup_mutations_freq_select, decreasing = T),
     main=paste("IN-group ",length(ingroup_mutations_freq_select)," mutations frequency distribution in the IN-group (",
                length(in_group_accessions_subset),") samples",sep=""),
     las=2,cex.axis=0.6,ylab="Frequency",xlab="")
dev.off()
getwd()


#FIND in how many lineages mutations are found in
mutation_freq=list()
dir.create("INGROUP-MUTATIONS-INFO-PLOTS")
for(mutation in ingroup_mutations){
  mutation_freq[[mutation]]=c()
  mutation_freq[[mutation]]=table(unlist(lapply(mutations_database, function(item,m){
    #cat(m,parent.frame()$i,"\n")
    if(m %in% item$mutations){return(item$lineage)}
  },m=mutation)))
  
  cat(names(mutation_freq[mutation]),"\n")
  
  pdf(file=paste("INGROUP-MUTATIONS-INFO-PLOTS/INGROUP_mutation_",mutation,"_freq_lineages.pdf",sep=""),
      height=4, width=20)
   shared_lineages_counts = sort(mutation_freq[[mutation]]/GISAID_lineages_freq[names(mutation_freq[[mutation]])],decreasing = T)
   plot(sort(shared_lineages_counts, decreasing = T)[1:100],
       main=paste("IN-group mutation ",mutation," frequency found in the top 100 PANGO lineages across the entire GISAID (",
                  length(mutations_database),")", sep=""),
       las=2,cex.axis=0.5, cex.main=1, ylab="% of lineage samples",xlab="")
  dev.off()
  
}

getwd()