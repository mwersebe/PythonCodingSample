#!/bin/python3
########################################
## Matthew Wersebe
## University of Oklahoma, Norman
## Coding Sample; Proficient in Python3
## July 19, 2022
## matthew.wersebe@ou.edu
##########################################

# Purpose: Use information in a gff file to calculate some descriptive statistics for a genome annotation. 

# Usage: run in conda environment, with required packages. 

# Import required libs for running the script 

import pandas

import os

import matplotlib as plt

import seaborn as sns

# Shell saves the day by parsing unimportant info from the gff file:

os.system("zcat GCF_021234035.1_SC_F0-13Bv2_genomic.gff.gz|grep -v '#'|cut -f 1,2,3,4,5 > fixed.gff")


# import the parsed gff into pandas:

feature_table = pandas.read_csv('fixed.gff', sep = '\t', names = ['Chrom', 'Source', 'Feature', 'Start', 'Stop'])

# lets look at average gene size for this species:

# Ssubset down to the genes listed in the table:

genes = feature_table[feature_table["Feature"] == "gene"]

# Calculate the length of the genes:

Length = genes["Stop"] - genes["Start"]

# Insert this back into the dateframe:
genes.insert(5, "Length", Length, True)

# Calc the mean
mean_length = genes["Length"].mean()

# Print the important info to the screen:
print('Mean gene length is', mean_length, 'base pairs.')

# Work iteratively over different feature types and calculate Mean lengths:


feature_types = ['gene', 'exon', 'lnc_RNA', 'tRNA']

for feature in feature_types:
	print('Doing work on', feature)
	table = feature_table[feature_table["Feature"] == feature]
	Length = table["Stop"] - table["Start"]
	table.insert(5, "Length", Length, True)
	mean_length = table["Length"].mean()
	print('Mean', feature, 'length is', mean_length, 'base pairs')


# Create a boxplot for each chromosome displaying the gene lengths:

# This just gets rid of the pesky unplaced contigs:

os.system("zcat GCF_021234035.1_SC_F0-13Bv2_genomic.gff.gz|grep -v '#'|cut -f 1,2,3,4,5|grep -v 'h1tg' > placed.gff")

#Repeating steps above with new gff style file:
placed_features = pandas.read_csv('placed.gff', sep = '\t', names = ['Chrom', 'Source', 'Feature', 'Start', 'Stop'])

genes = placed_features[placed_features["Feature"] == "gene"]

Length = genes["Stop"] - genes["Start"]

genes.insert(5, "Length", Length, True)

#Plot the data and save as a png file for later use. 

sns.boxplot(x = 'Length', y = 'Chrom', data = genes)

plt.pyplot.savefig('GeneLength.png')


