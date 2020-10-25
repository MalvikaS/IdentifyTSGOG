cTaG
============================
cTaG is a tool used to identify tumour suppressor genes(TSGs) and oncogenes (OGs) using somatic mutation data.

## Table of content

- [Description](##description)
- [Overview of cTaG](##overview-of-ctag)
- [Folder structure](##folder-structure)
- [Links](#links)

## Description

The model is built using mutation data from COSMIC (v79). The pre-processed data is used to generate features. Cancer Gene Census (CGC) is used to label genes as TSG or OG. The pan-cancer model is trained on these genes.



## Folder structure
The top-level directories contain code, data and output folders. 

### The directory layout

    .
    ├── code                   		# All the code for analysis
    ├── data
    │   ├── FeatureMat			# pre-processed data, feature matrix
    │   │   └── tissues		 	# Tissue-specific feature matrices
    │   └── Functional_Analysis		# List of genes used for function analysis
    ├── output
    │   ├── evalRandIter        	# Results for verying number of random iterations
    │   ├── FunctionalAnalysis  	# Results for functional analysis
    │   ├── MutSigCVComparisson 	# Results of MutSigCV and comparisson of mutation rates 
    │   ├── RandomForest		# Model for each CV and consensus of results
    │   └── TissueSpecificAnalyisis	# Tissue specific predictions and consensus of results
    └── README.md

The code folder containes all the files used for building the feature matrix, building the models and and the tissue-specific analysis.
    .
    ├── ...
    ├── code                   			# All the code for analysis
    │   ├── consolidateModelPredictions.py	# Consolidate results for pan-cancer analysis
    │   ├── consolidateTissueAnalysis.py	# Consolidate results for tissue-specific analysis
    │   ├── create_feature_mat_keepv2.py	# Creates feature matrix from COSMIC data
    │   ├── evaluate_randiter.py		# Evaluates performance different number of random iterations
    │   ├── identifyTSGOG.py			# Module for building feature matrix
    │   ├── plotMutRate.py			# Plots mutation rate for genes redicted by cTaG and MutSigCV
    │   ├── TSG_OG_Feat_final.py		# Loads feature matrices and build models
    │   └── TSG_OG_tissueAnalysis.py		# Creates tissue-specific feature matrices and makes predictions
    └── ...

## Links
