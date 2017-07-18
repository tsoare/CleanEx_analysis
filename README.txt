README for the CleanEx data and analysis:

The goal of CleanEx: Tissue-specific expression predictor is to predict whether a novel promoter sequence is likely to be expressed specifically in an organ or tissue of interest.
The input is a FASTA file (e.g. *.fasta, *.fa) of sequences that are 1000 base pairs long and the output is the probability that each sequence is likely to be "liver-specific" or
expressed primarily in the liver as opposed to the rest of the human body.

Under the hood, features are derived from the sequences (counts of N-grams), normalized, and fed into a logistic regression model to return the probability of being "liver-specific".
The accuracy in the test data was 60%, so testing novel sequences in this way is a useful filter before carrying forward promising sequences to be assessed experimentally in cell lines.

Working in both Python and R, I cleaned up the input data (promoter sequences and gene expression across 37 tissue types), derived the features, and tested a range of models.

Contents:
hg19_lnENr.fa - raw promoter sequences (obtained from EPD: http://epd.vital-it.ch/human/human_database.php)
rna_tissue.csv - raw gene expression (too big to be stored here, can be obtained from HPA: http://www.proteinatlas.org/)
input_data_top-bottom_classify.csv - cleaned input data: top 10% of genes with liver-specific expression and bottom 10% of genes without liver-specific expression
CleanEx_cleaning_promoter_sequences.ipynb - jupyter notebook for checking alphabet of promoter sequences
CleanEx_model_selection.ipynb - jupyter notebook for statistical analysis and model selection
bash_script_to_calc_PWM_match_in_MOODS.txt - bash script to calculate match scores for transcription factor binding site (TFBS) position weight matrices (PWMs) in MOODS (https://www.cs.helsinki.fi/group/pssmfind/)
parse_expression_data.R - R script for parsing expression data
MOODS_results.R - R script for compiling results of PWM-matching in MOODS
data_cleaning for statistical analysess.R - R script for cleaning data and taking top 10% of genes with liver-specific expression and bottom 10% of genes without liver-specific expression
graphs_for_presentation.R - R script for making figures for slide presentation

For more background information, please see cleanex.press.  Feel free to reach out with any questions!