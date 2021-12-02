#!/usr/bin/env nextflow
/*
========================================================================================
    nf-MIA-GISAID-pipeline
    =================================================================
    Github : https://github.com/kbessonov1984/nf-MIA-pipeline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

def helpMessagePipeline() {
  log.info"""
  ==================================================================
  ${workflow.manifest.name}   ~  version ${workflow.manifest.version}
  ==================================================================
  Usage:
  Given GISAID database JSON file, the typical command for running the analytical pipeline is as follows:
  
    nextflow run main.nf \\
     --GISAIDMetaDataFile "<path to the GISAID metadata file>" \\
     --GISAIDdatabase  "<path to the JSON database file>" \\
     --in_group_lineage "<PANGO Lineage defining the IN-Group>" \\
     --outdir <name of the output folder>
     --n_permutations 1000
     --analysis
     -profile conda
     
   To create or update GISAID previous database JSON file (define --GISAIDdatabase), the typical command for running the analytical pipeline is as follows:
  
    nextflow run main.nf \\
     --GISAIDMetaDataFile "<path to the GISAID metadata file>" \\
     --outdir <name of the output folder>
     --dbupdate_create
     -profile conda  -resume
     
     
  Options:
    --outdir              Output directory name where results will be stored
    --GISAIDMetaDataFile  Path to the GISAID meatadate file (metadata.tsv)
    --GISAIDdatabase      Path to the accessions to mutations and PANGO Lineages database in JSON format (see --dbupdate_create)
    --GISAID_fasta_file   Path to the GISAID MSA file in FASTA format (msa_0822.fasta)
    --reference_fasta     ${params.reference_fasta}
    --analysis            Perform stastistical analysis on IN- and OUT-Groups (requires --GISAIDdatabase definition)
    --dbupdate_create     Create GISAID accession to mutations and PANGO Lineages database from scratch (requires --GISAIDMetaDataFile  and --GISAID_fasta_file)
    --in_group_lineage    PANGO Lineage defining the IN-Group
    --location_string     Optional location filter for GISAID medata data (e.g. "Canada"). If not defined no metadata filtering (i.e. global context)
    --n_permutations      Number of permutations for step-down MAXT p-value caluclation (default: ${params.n_permutations})
    --start_date          Optional start date filter for GISAID metadata in Year-Month-Day format (e.g. 2021-01-30)
    --end_date            Optional end date filter for GISAID metadata in Year-Month-Day format (e.g. 2021-02-15). If not defined, the most recent GISAID metadata record is used as upper time limit
  """
}

def checkFileExists(file_path) {
  if (file_path){
      f = file(file_path)
      if ( !f.isFile() || !f.exists() ) {
        exit 1, "File '$file_path' does not exist!"
      }
  }
}

if (params.help){
  helpMessagePipeline()
  exit 0
}


process CreateSampleBatches{
    publishDir "${params.outdir}/.sample_batches",
    pattern: "*.txt",
    mode: 'symlink'
    
    
    input:
    val GISAIDMetaDataFile
    val GISAIDdatabase
    
    output:
    path 'samples2extact_batch*.txt', emit: outfiles
    stdout emit: messages
    
    script:
    //log.info "Creating batches of sample IDs ..."
    """
    Rscript ${baseDir}/bin/select_samples.R $GISAIDMetaDataFile $GISAIDdatabase
    """
 
}

process GetVCFSample{
    echo false
    
    publishDir "${params.outdir}",
    mode: 'symlink'
    
    input:
    val file_batch_ids
    path GISAID_fasta_file
    path reference_fasta
    
    output:
    path "vcf_files/*.vcf", emit: vcf_files
    //stdout emit: messages
    
    script:
    //log.info "Get mutations for batch   ${baseDir}/bin/get_mutations_vcf.sh ${file_batch_ids} ${GISAID_fasta_file} ${reference_fasta} ..."
    """
    # echo -n $file_batch_ids
     ${baseDir}/bin/get_mutations_vcf.sh ${file_batch_ids} ${GISAID_fasta_file} ${reference_fasta}
    """
    
}


process GenerateSample2MutationsDBparts{
    echo false
    maxRetries 2
    
    publishDir "${params.outdir}/mutations_db_parts",
    pattern: "*.json",
    mode: 'symlink'
    
    input:
    val vcf_path_batch_file
    val GISAIDMetaDataFile
    
    output:
    path "*.json", emit: json_database_parts

    script:
    //log.info "$vcf_path_batch_file"
    """
    ${baseDir}/bin/mutations_database_generation.R ${vcf_path_batch_file} ${GISAIDMetaDataFile}
    """
}

process GenerateSample2MutationsDBaggregate{
    echo false
    
    input:
    val vcf_path_batch_files
    val previous_mutations_database
    
    output:
    path "*.json", emit: json_database_raw_file

    
    script:
    vcf_path_batch_files = vcf_path_batch_files.join(' ')
    """
    ${baseDir}/bin/mutations_db_parts_aggregate.R ${vcf_path_batch_files} ${previous_mutations_database}
    """
}

process GenerateSample2MutationsDBlineagesUpdate{
    echo false
    
    publishDir "${params.outdir}",
    mode: 'copy'
    
    input:
    val GISAIDMetaDataFile
    val mutations_database
    
    output:
    path "*.json", emit: json_database_file
    path "*.Rdata", emit: Rdata_database_file

    
    script:
    """
    ${baseDir}/bin/mutations_db_homogenize.R ${GISAIDMetaDataFile} ${mutations_database}
    """
}

process GenerateVCFBatchLists{
    echo false
    
    publishDir "${params.outdir}/.vcf_file_batches",
    pattern: "*.txt",
    mode: 'symlink'
    
    input:
    val vcf_files_path_list
    
    output:
    path "*.txt", emit: vcf_file_path_batches
    
    
    script:
    //log.info " ${vcf_files_path_list} ..."
    log.info "${params.outdir}" 
    //FileOut = file()
    File FileOut = new File("${baseDir}/${params.outdir}/vcf_list_files.txt")
    FileOut.withWriter{out -> vcf_files_path_list.each {out.println it}}
    //println vcf_files_path_list.getClass() 

    """
    ${baseDir}/bin/create_vcf_batches.R  ${baseDir}/${params.outdir}/vcf_list_files.txt 
    """  
}

process PerformAnalysisGenerateReport{
    
    publishDir "${params.outdir}/analysis",
   // pattern: "*.{pdf,txt}",
    mode: 'copy'
    
    input:
    val GISAIDMetaDataFile
    val GISAIDdatabase
    val ingroupPangoLineage
    val n_permutations
    val location_string
    val start_date
    val end_date
 
    output:
    path '*', emit: outfiles
    stdout emit: messages
    
    script:
    """
    Rscript ${baseDir}/bin/stat_analysis.R $GISAIDMetaDataFile $GISAIDdatabase $ingroupPangoLineage $n_permutations $location_string $start_date $end_date
    """

}
workflow {
   checkFileExists(params.GISAIDMetaDataFile)
   checkFileExists(params.GISAIDdatabase)
   
   if (params.dbupdate_create){
   checkFileExists(params.reference_fasta)
   
   CreateSampleBatches("${params.GISAIDMetaDataFile}", "${params.GISAIDdatabase}")
   //CreateSampleBatches.out.messages.view()
   

   //log.info "Obtaining mutations for each sample batch ..."
   //channel.fromPath("results/.sample_batches/*.txt").view()
   GetVCFSample(CreateSampleBatches.out.outfiles.flatMap(), "${params.GISAID_fasta_file}", "${params.reference_fasta}")
  
   GenerateVCFBatchLists(GetVCFSample.out.vcf_files.collect())
   //GenerateVCFBatchLists(channel.fromPath("results/vcf_files/*.vcf").take(1000).toList()) // TEST CODE
   GenerateSample2MutationsDBparts(GenerateVCFBatchLists.out.vcf_file_path_batches.flatMap(),"${params.GISAIDMetaDataFile}")
   GenerateSample2MutationsDBaggregate(GenerateSample2MutationsDBparts.out.json_database_parts.collect(), "${params.GISAIDdatabase}")
   GenerateSample2MutationsDBlineagesUpdate("${params.GISAIDMetaDataFile}", GenerateSample2MutationsDBaggregate.out.json_database_raw_file)
   }
   
   if (params.analysis && params.GISAIDdatabase){
     if(params.in_group_lineage == null){exit 1, "Parameter 'in_group_lineage' not defined!"}
     PerformAnalysisGenerateReport("${params.GISAIDMetaDataFile}", "${params.GISAIDdatabase}", 
       "$params.in_group_lineage", "$params.n_permutations", "$params.location_string", "$params.start_date", "$params.end_date")
       
   }
   
   if ( params.dbupdate_create == null &&  params.analysis == null ){
     log.error """MIssing parameter(s). Specify either --dbupdate_create, --analysis or both"""  
   }
   
}

//EXECUTION PIPELINE LAUNCH PARAMS
log.info """
=======================================================
${workflow.manifest.name} v${workflow.manifest.version}
Results Dir: ${file(params.outdir)}
GISAIDMetaDataFile: ${params.GISAIDMetaDataFile}
GISAIDdatabase: ${params.GISAIDdatabase}
======================================================="""

workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir)}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

/* 

nextflow run  main.nf --GISAIDMetaDataFile "/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/27-09-2021/metadata.tsv" --GISAIDdatabase "/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/25-08-2021/vcf_files/mutations_database_final.json" --GISAID_fasta_file "/home/CSCScience.ca/kbessono/SARS-COV-2/Shared_Data/GISAID/27-09-2021/msa_0927/msa_0927.fasta"  --dbupdate_create -profile conda -resume 

*/
