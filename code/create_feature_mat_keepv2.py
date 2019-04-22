#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 19:39:13 2018

@author: malvika
Generate stats for simulating gene using COSMIC data
"""

import os
import pandas as pd
import numpy as np
import identifyTSGOG as ito
import glob
import pickle

# TODO: Set path
PATH = "/home/malvika/Documents/code/IdentificationOfTSG-OG"
DATAPATH = "/home/malvika/Documents/code/data/IdentificationOfTSG-OG"
#PATH = "/home/symec-02-01/Documents/IdentificationOfTSG-OG"
#DATAPATH = "/home/symec-02-01/Documents/data/IdentificationOfTSG-OG"
os.chdir(PATH)

# Read COSMIC data and store in dataframe
filename = DATAPATH + "/COSMIC/v79/" + "CosmicWGS_MutantExport.tsv"
data = ito.parseCDmutCosmic(filename, data=None)

# Filter samples with large number of mutations
data = ito.filterSamples(data, numMutations=2000)
# Filter data to keep rows without missing data and other criteria
data = ito.keep_v2(data)

# Change GRCh 38 to 37
# TODO: Change path and filename if required
os.makedirs(PATH + "/data/liftover/keepv2", exist_ok=True)
os.chdir(PATH + "/data/liftover/keepv2")
liftOverIP = "LiftOverIP_v2.txt"
# Get LiftOver input and store it in file(liftOverIP) and also return as list
[indices, pos_list] = ito.getLiftOverInput(ito.getMissenseData(data), 38,
                                           liftOverIP)

# TODO: Give file created as input to LiftOver tool.
# Save 2 files, .bed and .err for GRcH 37 locations and error file with
# unmapped locations respectively.

# Read LiftOver output from files(converted and error file) and add to data
# TODO: change folder and file names of LiftOver output
os.chdir(PATH + "/data/liftover/Output/keepv2")
converFile = "LO_op_keepv2_20181203.bed"
errFile = "LO_op_keepv2_20181203.err"
[converList, converListIndex, errList] = ito.readLiftOverOutput(
                                    converFile, errFile, [indices, pos_list])
GRCh37 = [None] * len(data)
for idx, pos in zip(converListIndex, converList):
    GRCh37[idx] = pos
data["GRCh37"] = GRCh37

# Create input files for PolyPhen2
# TODO: change path
os.makedirs(PATH + "/data/PolyPhen2/keepv2", exist_ok=True)
os.chdir(PATH + "/data/PolyPhen2/keepv2")
[PP2list, PP2files_ip] = ito.getPolyPhen2input(data)

# TODO: Give the files generated as input to PolyPhen2.
# Save short file in path where files need to read from. The full, SNPs,
# and log file can be stored for fututre if needed.

# Reads PolyPhen2 output
# TODO: change path
os.chdir(PATH + "/data/PolyPhen2/keepv2/Output")
PP2_opfiles = glob.glob('*.txt')
PP2scores = ito.readPolyPhen2output(PP2_opfiles, PP2list)
PP2prob1 = [np.nan] * len(data)
PP2prob2 = [np.nan] * len(data)
for idx in PP2scores:
    if len(idx[1]) == 2:
        PP2prob1[idx[0]] = idx[1][0]
        PP2prob2[idx[0]] = idx[1][1]
    else:
        PP2prob1[idx[0]] = idx[1][0]
        PP2prob2[idx[0]] = np.nan
data["PP2_prob1"] = PP2prob1
data["PP2_prob2"] = PP2prob2

# Find splicing mutations
os.chdir(DATAPATH)
SSitefile = "SpliceSites.txt"
spliceSite = ito.getSplicingMut(data, SSitefile)
data["SpliceSite"] = spliceSite

# Save filtered and complete data
os.makedirs(DATAPATH + "/FeatureMat/keepv2", exist_ok=True)
os.chdir(DATAPATH + "/FeatureMat/keepv2")
fname = "filteredData.pkl"
with open(fname, 'wb') as f:
    pickle.dump(data, f)

# read TSG and OG gene lists
os.chdir(PATH + "/data")
TSGfile = "onlyTSG.txt"
OGfile = "onlyOG.txt"
with open(TSGfile) as f:
    TSGlist = [tsg.strip() for tsg in f.readlines()]
with open(OGfile) as f:
    OGlist = [og.strip() for og in f.readlines()]

# Create multiple feature matrices for different MAX_MULTIPLIER
os.makedirs(DATAPATH + "/FeatureMat/keepv2", exist_ok=True)
os.chdir(DATAPATH + "/FeatureMat/keepv2")
max_mulultipliers = [2]
for max_mul in max_mulultipliers:
    # Create feature matrix
    features_cd = ito.getCdMutFeatures_v2(data, max_mul)

    # Labels for all genes
    Label = [None] * len(features_cd)
    for idx, gene in enumerate(features_cd.index):
        if gene in TSGlist:
            Label[idx] = "TSG"
        elif gene in OGlist:
            Label[idx] = "OG"
        else:
            Label[idx] = "Neutral"
    # for all feature matrix versions
    features_cd["Label"] = Label

    # save as pickle and also write as a csv file
    # TODO: change the path to where you want the files to be saved
    featFName = "allGenes_keepv2_MM{:03d}.txt".format(max_mul)
    features_cd.to_csv(featFName, sep='\t')
    featFName = "feat_keepv2_MM{:03d}.pkl".format(max_mul)
    with open(featFName, 'wb') as f:
        pickle.dump(features_cd, f)


