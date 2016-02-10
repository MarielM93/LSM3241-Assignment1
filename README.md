# LSM3241 Assignment 1

**microarray data used:** GEO4498

**study title:** *"Modification of gene expression of the small airway epithelium in response to cigarette smoking"* by Harvey et al

**submitted by:** Maniclang Mariel Anne F [A0099537L]

## Assignment1_codev5.R
- R Script to initialize microarray analysis. Also performs gene annotation via biomaRt.
- Exports data into csv files

## sorterv3.py
- python script which merges and compares the results from RMA and MAS5 normalization
- ensure probe mapping csv files have been properly cut first via bash/shell (exact bash commands are included inside the python script for convenience)
- to export to a txt file, simply input "python sorterv3.py > filename.txt" into command line