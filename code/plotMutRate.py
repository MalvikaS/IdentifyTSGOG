#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:23:35 2019

@author: malvika
"""

import pandas as pd
import numpy as np
import os
from matplotlib import pyplot
from scipy import stats

PATH = "/home/malvika/Documents/code/IdentifyTSGOG"
os.chdir(PATH)

# ----------- Load data and add or rename required columns ----------- #
# MutSigCV predictions load
os.chdir(PATH + "/mutsigCV/all_set")
fname = "all_output.sig_genes.xlsx"
data_mutsig = pd.read_excel(fname, sheet_name="all_output.sig_genes")
# rename columns
data_mutsig.rename(columns={"gene": "Gene", 'u x 10^6': "u"}, inplace=True)
# add columns
data_mutsig["Consensus"] = [5] * len(data_mutsig)
data_mutsig["Label"] = ["Driver"] * len(data_mutsig)
data_mutsig["Source"] = ["MutSigCV"] * len(data_mutsig)

# Our predictions load
os.chdir(PATH + "/RandomForest/CV")
fname = "CVpredictions_oldfeat.xlsx"
data_old = pd.read_excel(fname, sheet_name="CVpredictions_oldfeat")
# rename columns
data_old.rename(columns={"Gene name": "Gene"}, inplace=True)
# add columns
data_old["Source"] = ["Our predictions"] * len(data_old)

# Old feature predictions load
os.chdir(PATH + "/RandomForest/CV")
fname = "CVpredictions.xlsx"
data_novel = pd.read_excel(fname, sheet_name="CVpredictions")
# rename columns
data_novel.rename(columns={"Gene name": "Gene"}, inplace=True)
# add columns
data_novel["Source"] = ["Our predictions"] * len(data_novel)

# Training data
fname = "Training_predictions.xlsx"
data_train = pd.read_excel(fname, sheet_name="ConsolidatedCV_results")
# rename columns
data_train.rename(columns={"Gene name": "Gene"}, inplace=True)
# add columns
data_train["Source"] = ["Training"] * len(data_train)

# ----------------- Data for plotting ---------------- ##
# Data requires following columns
# "Source" can be {MutSigCV, Training, Our predictions"}
cols = ["Gene", "u", "Consensus", "Label", "Source"]
data_all = pd.DataFrame(columns=cols)

# Filter data
# Our predictions
data_temp = data_novel[data_novel["Consensus"] >= 3]
data_temp = pd.merge(data_temp, data_mutsig[['Gene', 'u']], how="left",
                     on=["Gene"])
data_temp = data_temp[cols]
data_temp = data_temp.dropna()
data_all = pd.concat([data_all, data_temp])

# Training predictions
data_temp = data_train[data_train["Consensus"] >= 0]
data_temp = pd.merge(data_temp, data_mutsig[['Gene', 'u']], how="left",
                     on=["Gene"])
data_temp = data_temp[cols]
data_temp = data_temp.dropna()
data_all = pd.concat([data_all, data_temp])

# MutSigCV predictions
data_temp = data_mutsig[(data_mutsig["p"] <= 0.005) &
                        (data_mutsig["q"] <= 0.01)]
data_temp = data_temp[cols]
data_temp = data_temp.dropna()
data_all = pd.concat([data_all, data_temp])


# ---------------- Plot Histogram -------------------##
# plot histogram for MutSigCv (602 genes), our model (429) and training (208)
# MutSigCV: pvalue >=0.005 & q >= 0.01
# Our model: Consensus >=3
# Training: Consensus >=0
data_plot = data_all
minbin = min(data_plot["u"])
maxbin = max(data_plot["u"])
bins = np.linspace(minbin, maxbin, 50)
# plot
pyplot.hist(data_plot[data_plot["Source"] == "MutSigCV"]["u"],
            bins, alpha=0.5, label='MutSigCv')
pyplot.hist(data_plot[data_plot["Source"] == "Our predictions"]['u'],
            bins, alpha=0.5, label="Our model")
pyplot.hist(data_plot[data_plot["Source"] == "Training"]['u'],
            bins, alpha=0.5, label="Training")
pyplot.legend(loc='upper right')
pyplot.xlabel("Mutation rate (1e-6)", fontsize =14)
pyplot.ylabel("Frequency", fontsize =14)
pyplot.title("Histogram of mutation rates of genes predicted")
os.chdir(PATH + "/mutsigCV/all_set")
pyplot.savefig("mut_rate_hist_.jpg")
pyplot.close()

# --------------- Kolmogorov-Smirnov statistic -------------- ##
# Calculate statistic and pvalue between MutSigCV and (Our method + Training)
KS_stat, pval = stats.ks_2samp(data_plot[data_plot["Source"] ==
                                         "MutSigCV"]["u"],
                               data_plot[data_plot["Source"] !=
                                         "MutSigCV"]['u'])
print("KS statistic = {:0.3f}\np-value = {:0.3f}".format(KS_stat, pval))

# Compare training dirstribution to MutSigCv and our prediciton distributions
for source in ["MutSigCV", "Our predictions"]:
    KS_stat, pval = stats.ks_2samp(data_plot[data_plot["Source"] ==
                                             source]["u"],
                                   data_plot[data_plot["Source"] ==
                                             "Training"]['u'])
    print("Kolmogorov statistic for {} = {:0.3f}".format(source, KS_stat))
    print("P-value for {} = {:0.3f}".format(source, pval))


# ------------- Plot fraction of genes for varying mutation rates ------#
num_genes = 60
# Consensus defines how our model predictions are filtered
consensus = 5
# MutSigCV predictions ranked and filtered
cols = ["Gene", "u", "Consensus", "Label", "Source"]
data_plot = data_mutsig.sort_values(by=["q", "p"]).iloc[:num_genes, :]
data_plot = data_plot[cols]
# Data concatenated
data_plot = pd.concat([data_plot,
                       data_all[((data_all["Source"] == "Training") |
                                (data_all["Consensus"] >= consensus)) &
                                (data_all["Source"] != "MutSigCV")]])
data_plot["log_u"] = round(np.log(data_plot["u"]), 2)
# Populate fractions to be plotted
temp = []
for source, u in zip(data_plot["Source"], data_plot["log_u"]):
    tot = len(data_plot[(data_plot["Source"] == source)]["log_u"])
    fraction = len(data_plot[(data_plot["Source"] == source) &
                             (data_plot["log_u"] <= u)]["log_u"]) / tot
    temp.append(fraction)
data_plot["Fraction"] = temp

# plot
fig = pyplot.figure(figsize=(11, 8))
ax1 = fig.add_subplot(111)
for source, c, m in zip(list(set(data_plot["Source"])), ['c', 'b', 'r'],
                        ['x', 'o', 's']):
    temp = data_plot[data_plot["Source"] == source]
    ax1.scatter(temp["log_u"], temp["Fraction"], label=source,
         color=c, marker=m)
pyplot.xlabel('Log mutation rate', fontname='Calibri', fontsize =16)
pyplot.ylabel('Fraction of genes predicted below mutation rate',
              fontname='Calibri', fontsize =16)
handles, labels = ax1.get_legend_handles_labels()
pyplot.xticks(fontsize=12)
pyplot.xticks(fontsize=12)
lgd = ax1.legend(handles, labels, loc='upper left',
                 bbox_to_anchor=(0.01, 1), fontsize =16)
ax1.grid('on')
os.chdir(PATH + "/mutsigCV/all_set")
pyplot.savefig('log_fraction_scatter_cv60.png')
pyplot.close()



# Compare MutSigCV to Training set but vary threshold for q-value
for thres in [0.010, 0.001, 0.0001, 0.00001]:
    # Filter data
    fil_mutsig = data_mutsig[(data_mutsig["p"] <= .005) &
                             (data_mutsig["q"] <= thres)]["u"]
    # Calculate KS stat and pval and print
    KS_stat, pval = stats.ks_2samp(fil_mutsig,
                                   data_plot[data_plot["Source"] ==
                                             "Training"]['u'])
    print("Number of genes = {}".format(len(fil_mutsig)))
    print("Kolmogorov statistic for q-value threshold {} = {:0.3f}".format(
            thres, KS_stat))
    print("p-value for q-value threshold {} = {:0.3f}".format(thres, pval))
    # Plot
    # For histogram uncomment next 3 lines
#    pyplot.hist(fil_mutsig, bins, alpha=0.5, label='MutSigCv')
#    pyplot.hist(data_plot[data_plot["Source"] == "Our predictions"]['u'],
#                bins, alpha=0.5, label="Our model")
#    pyplot.hist(data_plot[data_plot["Source"] == "Training"]['u'], bins,
#                alpha=0.5, label="Training")
#    pyplot.legend(loc='upper right')
#    pyplot.xlabel("Mutation rate (1e-6)")
#    pyplot.ylabel("Frequency")
#    pyplot.title("Histogram of mutation rates of genes predicted")
#    os.chdir(PATH + "/mutsigCV/all_set")
#    pyplot.savefig("hist_filtered_qval{}.jpg".format(thres))
#    pyplot.close()
    # For sorted scatter uncomment next 3 lines
#    pyplot.scatter(x=range(len(fil_mutsig)), y=sorted(fil_mutsig),
#                   alpha=0.5, label='MutSigCv')
#    pyplot.scatter(x=range(len(data_plot[data_plot["Source"] ==
#                                         "Our predictions"]['u'])),
#                   y=sorted(data_plot[data_plot["Source"] ==
#                                      "Our predictions"]['u']),
#                   alpha=0.5, label="Our model")
#    pyplot.scatter(x=range(len(data_plot[data_plot["Source"] ==
#                                         "Training"]['u'])),
#                   y=sorted(data_plot[data_plot["Source"] ==
#                                      "Training"]['u']),
#                   alpha=0.5, label="Training")
#    pyplot.legend(loc='upper right')
#    pyplot.xlabel("Genes")
#    pyplot.ylabel("Mutation rate (1e-6)")
#    pyplot.title("Sorted scatter plot of mutation rates of genes predicted")
#    os.chdir(PATH + "/mutsigCV/all_set")
#    pyplot.savefig("scatter_filtered_qval{}.jpg".format(thres))
#    pyplot.close()


# -------------------- Plot of fractions ------------------------ #
# Plots the fraction of genes predicted upon the total genes predicted
# X-axis: Varying threshold of mutation rates
# Y-axis: Fraction of genes predicted
# MutSigCV:
#  For c_cutoff == 5,4,3 -> all 602 genes predicted considered (p<=0.005,
#  q<=0.01)
#  For c_cutoff == 126, 60, 161, 429 -> genes predicted were ranked and top
#  c_cutoff value number of genes considered
# Our predictions:
#  For c_cutoff == 5,126, 60 -> genes with consensus >= 5 considered
#  For c_cutoff == 4, 161 -> genes with consensus >= 4 considered
#  For c_cutoff == 3, 429 -> genes with consensus >= 3 considered
# Training:
#  All 208 genes were considered

# Consensus list defined the way the data is filtered and ploted
consensus_list = [5, 4, 3, 126, 60, 161, 429]
for c_cutoff in consensus_list:
    if c_cutoff in [126, 60, 161, 429]:
        # MutSigCV predictions ranked and filtered
        cols = ["Gene", "u", "Consensus", "Label", "Source"]
        data_plot = data_mutsig.sort_values(by=["q", "p"]).iloc[:c_cutoff, :]
        data_plot = data_plot[cols]
        # Consensus defines how our model predictions are filtered
        if c_cutoff == 126 or c_cutoff == 60:
            consensus = 5
        elif c_cutoff == 161:
            consensus = 4
        elif c_cutoff == 429:
            consensus = 3
        # Data concatenated
        data_plot = pd.concat([data_plot,
                               data_all[((data_all["Source"] == "Training") |
                                        (data_all["Consensus"] >= consensus)) &
                                        (data_all["Source"] != "MutSigCV")]])
    else:
        # Data filtered to retain what is considered for plotting
        data_plot = data_all
        data_plot = data_plot[(data_plot["Source"] == "Training") |
                              (data_plot["Consensus"] >= c_cutoff)]
    data_plot["log_u"] = round(np.log(data_plot["u"]), 2)

    # Compare training dirstribution to MutSigCv and our prediciton
    # distributions
    for source in ["MutSigCV", "Our predictions"]:
        KS_stat, pval = stats.ks_2samp(data_plot[data_plot["Source"] ==
                                                 source]["u"],
                                       data_plot[data_plot["Source"] ==
                                                 "Training"]['u'])
        print("Kolmogorov statistic for {} ({} genes)= {:0.3f}".format(source,
              len(data_plot[data_plot["Source"] == source]["u"]), KS_stat))
        print("P-value for {} ({} genes)= {:0.3f}".format(source,
              len(data_plot[data_plot["Source"] == source]["u"]), pval))

    # Define variables to bin data
    tot_size = len(data_plot)
    num_bins = 10
    size = int(round(tot_size / num_bins, 0))
    # Calculate threshold for bins by sorting mutation rate and find
    # equidistant u values
    thres_list = [round(u, 2) for idx, u in enumerate(sorted(data_plot["log_u"]))
                  if (idx + 1) % size == 0]
    # Last bin threshold is max mutation rate
    if len(thres_list) < num_bins:
        thres_list.append(round(max(data_plot["log_u"]), 2))
    else:
        thres_list[num_bins-1] = round(max(data_plot["log_u"]), 2)
#    print(c_cutoff, tot_size, size, len(thres_list))
    # Initialise variable and populate table of fractions to be plotted
    cols = list(sorted(set(data_plot["Source"])))
    data_frac = pd.DataFrame(index=thres_list,
                             columns=cols)
    for thres in thres_list:
        # Fractions calculated for all 3 (MutSigCv, Training, Our model)
        for source in cols:
            tot = len(data_plot[(data_plot["Source"] == source)]["log_u"])
            fraction = len(data_plot[(data_plot["Source"] == source) &
                                     (data_plot["log_u"] <= thres)]["log_u"]) / tot
            data_frac.loc[thres, source] = fraction
    # plot
    fig = pyplot.figure(figsize=(11, 8))
    ax1 = fig.add_subplot(111)
#    ax1.plot(range(num_bins), data_frac["Training"], label='Training',
#             color='c', marker='x')
#    ax1.plot(range(num_bins), data_frac["MutSigCV"], label='MutSigCV',
#             color='b', marker='o')
#    ax1.plot(range(num_bins), data_frac["Our predictions"],
#             label='Our model', color='r', marker='s')
    ax1.plot(thres_list, data_frac["Training"], label='Training',
             color='c', marker='x')
    ax1.plot(thres_list, data_frac["MutSigCV"], label='MutSigCV',
             color='b', marker='o')
    ax1.plot(thres_list, data_frac["Our predictions"],
             label='Our model', color='r', marker='s')
#    pyplot.xticks(range(num_bins), tuple(thres_list))
    pyplot.xlabel('Mutation rate')
    pyplot.ylabel('Fraction of genes predicted below mutation rate threshold')
    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='upper left',
                     bbox_to_anchor=(0.01, 1))
    ax1.grid('on')
    os.chdir(PATH + "/mutsigCV/all_set")
    pyplot.savefig('log_fraction_plot_cv{}.png'.format(c_cutoff))
    pyplot.close()

    # Data 
    

# Plot Boxplot
# read TSG and OG gene lists
os.chdir("/home/malvika/Documents/code/IdentificationOfTSG-OG/data")
TSGfile = "onlyTSG.txt"
OGfile = "onlyOG.txt"
with open(TSGfile) as f:
    TSGlist = [tsg.strip() for tsg in f.readlines()]
with open(OGfile) as f:
    OGlist = [og.strip() for og in f.readlines()]

data = data_plot
dropidx = [idx for val, idx in zip(data["Gene"], data.index) if val not in TSGlist and val not in OGlist]
data.drop(dropidx)
pyplot.boxplot(data=data, x="u")
pyplot.boxplot(data_plot[data_plot["Source"] == "Our predictions"]['u'],)
pyplot.boxplot(data_plot[data_plot["Source"] == "Training"]['u'],)

boxplot = data.boxplot(column="u", by="Source")

KS_stat_, pval_ = stats.ks_2samp(data[data["Label"] == "Driver"]["u"],
                                 data[(data["Label"] == "TSG") |
                                         (data["Label"] == "OG")]['u'])

