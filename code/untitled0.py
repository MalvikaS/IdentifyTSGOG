#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 17:33:04 2019

@author: malvika
"""

from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
from inspect import signature
from sklearn.metrics import average_precision_score
import pandas as pd
import pickle as pkl
import numpy as np
import os
from matplotlib import pyplot
from scipy import stats

PATH = "/home/malvika/Documents/code/IdentifyTSGOG/"
os.chdir(PATH)

fig, ax = plt.subplots(3, 2, figsize=(6, 6), dpi=500)
st = fig.suptitle('Precision-Recall curves', fontsize="x-large")
for cv in range(0, 5):
    fpath = PATH + "RandomForest/CV_/{}/".format(cv)
    os.chdir(fpath)
    df = pd.read_csv("TrainingTest_predictions", sep="\t")
    df = df[df.Data == "test"]
    df[df=="TSG"] = 1
    df[df=="OG"] = 0
    precision, recall, _ = precision_recall_curve(list(df.Label),
                                                  list(df.Predictions))
    average_precision = average_precision_score(list(df.Label),
                                                list(df.Predictions))
    # In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
    step_kwargs = ({'step': 'post'}
                   if 'step' in signature(plt.fill_between).parameters
                   else {})
    c = cv % 2
    r = cv // 2
    ax[r, c].step(recall, precision, color='b', alpha=0.1,
             where='post')
    ax[r, c].fill_between(recall, precision, alpha=0.1, color='b', **step_kwargs)
    ax[r, c].set(xlabel='Recall', ylabel='Precision')
    ax[r, c].set_title("CV{}: AP={:0.3f}".format(cv + 1,
                          average_precision))
plt.tight_layout()
st.set_y(0.95)
fig.subplots_adjust(top=0.85)
fig.delaxes(ax[2, 1])

fname = PATH + "RandomForest/CV/PR_curve_test.png"
fig.savefig(fname)
