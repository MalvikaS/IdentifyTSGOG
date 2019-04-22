#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 20:37:25 2018

@author: malvika
Loads data required and creates feature matrix specific for tissues.
Loads the Pan cancer model and standard scaler to classifies genes for a
given each tissue
"""
import os
import pickle
import glob
import pandas as pd
import identifyTSGOG as ito

# TODO: Set path
PATH = "/home/malvika/Documents/code/IdentificationOfTSG-OG"
DATAPATH = "/home/malvika/Documents/code/data/IdentificationOfTSG-OG"
os.chdir(PATH)
#PATH = "/home/symec-02-01/Documents/IdentificationOfTSG-OG"
#DATAPATH = "/home/symec-02-01/Documents/data/IdentificationOfTSG-OG"
#os.chdir(PATH)
folderPath = "/TSG_OG_classifier/tissues"
os.chdir(DATAPATH + "/FeatureMat/keepv2")
fname = "filteredData.pkl"
with open(fname, 'rb') as f:
    data = pickle.load(f)

os.makedirs(DATAPATH + "/FeatureMat/keepv2/tissues", exist_ok=True)
os.chdir(DATAPATH + "/FeatureMat/keepv2/tissues")
samp_thresh = 1000
for tissue in set(data["Primary site"]):
    num_samp = len(set(data[data["Primary site"] == tissue]["ID_sample"]))
    print("{} : {}".format(num_samp, tissue))
    if num_samp >= samp_thresh:
        data_tiss = data[data["Primary site"] == tissue]
        features_cd = ito.getCdMutFeatures_v2(data_tiss, 2)
        features_cd["Label"] = ["Unlabeled"] * len(features_cd)
        featFName = "feat_{}.pkl".format(tissue)
        with open(featFName, 'wb') as f:
            pickle.dump(features_cd, f)

Kfolds = 5
sc = [None] * Kfolds
rfc = [None] * Kfolds
for fold in range(Kfolds):
    os.chdir(PATH + "/TSG_OG_classifier/RandomForest/CV/{}".format(fold))
    # Load scaler fit
    fname = "cosmicStdScale.pkl"
    with open(fname, 'rb') as f:
        sc[fold] = pickle.load(f)
    # Load the model
    fname = "randIterRFv5Model_allFeat.pkl"
    with open(fname, 'rb') as f:
        rfc[fold] = pickle.load(f)

os.chdir(DATAPATH + "/FeatureMat/keepv2/tissues")
tiss_files = glob.glob('*.pkl')
for file in tiss_files:
    os.chdir(DATAPATH + "/FeatureMat/keepv2/tissues")
    with open(file, 'rb') as f:
        features_cd = pickle.load(f)
    # Drop rows where all entries are Nan
    tsg_og_feat = features_cd[:-1].dropna(subset=list(features_cd.columns[0:-1]))

    for fold in range(Kfolds):
        # Standard scaling
        all_std = pd.DataFrame(sc[fold].transform(tsg_og_feat.iloc[:, 0:-1]),
                               index=tsg_og_feat.index,
                               columns=tsg_og_feat.iloc[:, 0:-1].columns)

        # Prediction of novel driver genes
        novel_pred = rfc[fold].predict(all_std)
        novel_prob = pd.DataFrame(rfc[fold].predict_proba(all_std),
                                  index=all_std.index,
                                  columns=rfc[fold].classes_).sort_values(
                                          by=["TSG"], ascending=False)
        novel_logp = pd.DataFrame(rfc[fold].predict_log_proba(all_std),
                                  index=all_std.index,
                                  columns=rfc[fold].classes_).sort_values(
                                          by=["TSG"], ascending=False)
        # Print to file
        os.makedirs("{}{}/{}".format(PATH, folderPath, file[5:-4]),
                    exist_ok=True)
        os.chdir("{}{}/{}".format(PATH, folderPath, file[5:-4]))
        filename = "prob_{}{}.txt".format(file[5:-4], fold)
        novel_prob.to_csv(filename, sep="\t", header=True, index=True)
