ARCHIVE (if present)
> Scripts and files that are either outdated or unused in final analyses; stored here as back-up

DATA
> Files that are necessary for running MATLAB script to calculate granularity and mean affect for each data set
> Files that are necessary for running R scripts to perform meta-analysis, integrative data analysis (IDA), and sensitivity (power) analyses

LIWC (license/download and general instructions available at http://liwc.wpengine.com/)
> Custom LIWC2015 dictionaries created for Studies 2 and 3 based on (translations) of Study 1 themes
> Outputs from scoring Study 2 and 3 data using custom LIWC dictionaries (please note that some columns/variables have been removed and the file naming modified from the LIWC defaults)

MEH (free download and general instructions available at https://www.ryanboyd.io/software/meh/)
> Settings used to analyze Study 1 texts and generate frequency list to assist with Study 3 LIWC dictionary translation
> OUTPUT_STUDY 1 subfolder includes files for all 27 MEH/PCA analyses (>20, >25, >30 words; 100, 150, 200 unigrams; 5, 10, 15 components)

SCRIPTS
> MATLAB and R scripts used in data analysis
> FUNCTIONS subfolder includes MATLAB functions that are necessary for running scripts
- granularity_affect measures.m = uses rating data to calculate granularity and mean affect for negative and positive emotions; run separately for each dataSet
- Study#_data.m = cleans texts for each dataSet, removing dropped participants and event descriptions below minimum word count; outputs texts for use in MEH
- Exploratory PCA_parallel analysis.R = uses MEH outputs to perform initial PCA to determine how many components to extract per binary matrix
- PCA_score texts_changepoints.R = uses MEH outputs to a) perform PCA on binary matrix, b) score texts (via verbose matrix) for extracted PCA components, and c) run changepoints analysis to threshold texts considered to evidence themes
- Study1_analyses.m = uses outputs from 'granularity_affect measures.m' and 'PCA_score texts_changepoints.R' to a) calculate diversity of themes for each participant in dataSet1, b) run regression models for the relationship between granularity and thematic diversity, and c) generate figures for study results as well as example participants
- Study2_analyses.m = uses outputs from 'granularity_affect measures.m' and LIWC2015 to a) calculate diversity of themes for each participant in dataSet2, b) run regression models for the relationship between granularity and thematic diversity, and c) generate figures for study results
- Study3_analyses.m = uses outputs from 'granularity_affect measures.m' and LIWC2015 to a) calculate diversity of themes for each participant in dataSet3, b) run regression models for the relationship between granularity and thematic diversity, and c) generate figures for study results
- Study1-3_IDA.R = combines standardized data from all three studies in an integrative data analysis
- Study1-3_meta-analysis.R = performs fixed effects meta-analysis of effect sizes from Studies 1-3, including meta-analytic sensitivity/power analysis
- Study1-3_sensitivity analyses.R = performs sensitivity (power) analyses for regression analyses (both single-study and IDA)