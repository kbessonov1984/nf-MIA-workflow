# MIA (Mutation Importance Assigner)

Mutation Importance Assigner SARS-COV-2 analytical pipeline estimates mutation "importance" by estimating its association strength to the selected IN-group against the OUT-group GISAID database samples (https://www.gisaid.org/). The aim is to facilitate more effective and objective discrimination of SARS-COV-2 lineages using mutation sets under different data selection criteria (time, geography, PANGO lineage) in light of the highly dynamic everchanging SARS-COV-2 data landscape.


# Algorithm

1. Define the IN-Group that is most frequently represented by the PANGO Lineage (e.g. B.1.1.7 or B.1.617.2). 
2. Define OUT-Group represented by the remainder GISAID data randomly sampled with sample size equal to IN-Group. 
3. Encode each sample group membership as either 1 or 0 depending on whether the sample is part of the IN or OUT Group (1 or 0, respectively). The group membership defines the response variable Y 
3. The predictors (X) are defined as the IN-Group mutations/SNVs that are found at least in the 10% of the IN-Group samples. The mutation presence/absence encoded as 0 and 1 for each sample, respectively.
4. Fit univariate logistic regression model for each mutation and extract the t-test statistic values corresponding to each beta coefficient (Tobs)
5. Permute the Y variable (but not X) and re-fit the same logistic regression model on permuted data. We collect the t-test statistic values (Tpermut). 
6. After B number of permutations, we obtain an empirical distribution of the Tpermut values and calculate FWER adjusted p-values using the step-down maxT algorithm for each mutation. The adjusted p-values for the maxT procedure is discussed on p. 50 and 114 of Westfall & Young (1993) [1].

# Install

Install NextFlow pipeline manager natively (not as conda environment) by running `curl -s https://get.nextflow.io | bash` in the destination folder. A binary file `nextflow` will be created. In addition, install conda if the pipeline will run under conda environment (recommended).

# Requirements

* Nextflow >= 21.04.0
* Conda environment
* R language >= 4.0
* R libraries: stringr, islasso, jsonlite, corrplot

# Usage
Note: we recommend running pipeline under conda profile (`-profile conda`) and in cache mode in case of failure (`-resume`) especially during the accessions to mutations and lineages database update or creation stage.
```
N E X T F L O W  ~  version 21.04.3
Launching `main.nf` [focused_hugle] - revision: f7b43ae852

  ==================================================================
  nf-MIA-GISAID-pipeline   ~  version 1.0.0
  ==================================================================
  Usage:
  Given GISAID database JSON file, the typical command for running the analytical pipeline is as follows:
  
    nextflow run main.nf \
     --GISAIDMetaDataFile "<path to the GISAID metadata file>" \
     --GISAIDdatabase  "<path to the JSON database file>" \
     --in_group_lineage "<PANGO Lineage defining the IN-Group>" \
     --outdir <name of the output folder>
     --n_permutations 1000
     --analysis
     -profile conda
     
   To create or update GISAID previous database JSON file (define --GISAIDdatabase), the typical command for running the analytical pipeline is as follows:
  
    nextflow run main.nf \
     --GISAIDMetaDataFile "<path to the GISAID metadata file>" \
     --outdir <name of the output folder>
     --dbupdate_create
     -profile conda  -resume
     
     
  Options:
    --outdir              Output directory name where results will be stored
    --GISAIDMetaDataFile  Path to the GISAID meatadate file (metadata.tsv)
    --GISAIDdatabase      Path to the accessions to mutations and PANGO Lineages database in JSON format (see --dbupdate_create)
    --GISAID_fasta_file   Path to the GISAID MSA file in FASTA format (msa_0822.fasta)
    --reference_fasta     references/NC_045512.2.fasta
    --analysis            Perform stastistical analysis on IN- and OUT-Groups (requires --GISAIDdatabase definition)
    --dbupdate_create     Create GISAID accession to mutations and PANGO Lineages database from scratch (requires --GISAIDMetaDataFile  and --GISAID_fasta_file)
    --in_group_lineage    PANGO Lineage defining the IN-Group
    --location_string     Optional location filter for GISAID medata data (e.g. "Canada"). If not defined no metadata filtering (i.e. global context)
    --n_permutations      Number of permutations for step-down MAXT p-value caluclation (default: 1000)
    --start_date          Optional start date filter for GISAID metadata in Year-Month-Day format (e.g. 2021-01-30)
    --end_date            Optional end date filter for GISAID metadata in Year-Month-Day format (e.g. 2021-02-15). If not defined, the most recent GISAID metadata record is used as upper time limit
```
# Output
The pipeline produced a number of outputs in the `results` folder defined by the `--outdir ` including:
1. `INGROUP_<LINEAGE>_OVERALL_MAIN_report.txt ` - Overall report defining running parameters and composition of the IN- and OUT- Groups in terms of the sample sizes and PANGO Lineages composition
2.  `Pearson correlation_plot_of_all_in-group_mutations.pdf` - Pearson correlation plot showing pairwise correlation values between all IN-Group mutations in the entire dataset used for univariate and multivariate model fitting
3.  `INGROUP_univariate_single_mutations_stats.txt` - univariate analysis results on set of logistic regression models (single mutation per model)
4.  `INGROUP_multivariate_LASSO_mutations_constalations_stats.txt` - multivariate penalized regression analysis LASSO results on the selected set of mutations after univariate analysis (with betas > 0) together with adjusted p-value based on observed and permuted regressors Z-scores
5.  `INGROUP-MUTATION-REPORTS-INFO` Folder - for each IN-GROUP mutation find corresponding frequency plot in entire GISAID metadata (no filtering applied) allowing to see the "mutation exclusivity", if any, across all PANGO Lineages in the database. Only mutation abundances in the top 100 PANGO Lineages are shown.
6. `INGROUP_mutations_freq_distribution.pdf` - Mutation abundances in the selected IN-Group plot
7. `Top100_PANGO_LINEAGES_distribution_GISAID_percent_cummulative.pdf` - The top 100 most frequent PANGO Lineages found in the entire GISAID data expressed as a cumulative percent of the total samples
8.  `Top100_lineages_distribution_abs_freq.pdf` - The top 100 most frequent PANGO Lineages found in the entire GISAID data expressed as absolute sample counts
9.  `Top100_OUT_group_lineages_distribution.pdf` - The top 100 PANGO Lineages found in the OUT-Group subsample

# Examples
1. There is a need to identify mutations that are significant in the Canadian context for the Delta variant (B.1.617.2) in the time window from 2021-01-01 to the most recent record available in the supplied database. Given updated most recent GISAID data release, run the following command to assess mutations importance.

    ````
    nextflow run main.nf --GISAIDMetaDataFile "GISAID/25-08-2021/metadata.tsv"  --GISAIDdatabase "25-08-2021/mutations_database_final.json" --in_group_lineage "B.1.617.2" --location_string "Canada"  --start_date 2021-01-01  --outdir results --analysis -profile conda --n_permutations 1000
    ````
2. To update an existing database with new GISAID data downloaded from the GISAID portal run the following command. The update could take 5 or more days depending on the resources and amount of data.
    ```
    nextflow run  main.nf --GISAIDMetaDataFile "27-09-2021/metadata.tsv" --GISAIDdatabase "25-08-2021/mutations_database_final.json" --GISAID_fasta_file "27-09-2021/msa_0927/msa_0927.fasta"  --dbupdate_create -profile conda  --outdir updated_database -resume
    ```

# License
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# References
[1] Westfall, Peter H., and S. Stanley Young. 1993. Resampling-Based Multiple Testing: Examples and
Methods for p-Value Adjustment. New York: John Wiley & Sons.


