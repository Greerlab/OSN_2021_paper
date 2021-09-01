#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pandas as pd
import numpy as np
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.feature_selection import f_regression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import sys
import pickle 
from collections import Counter

def lr(files):
    model = LinearRegression()
    random_state = 12
    min_cells = 7
    df = pd.read_csv("data/transformed_coords_1to2.txt",sep = "\t", header=None, names=["OR", "X", "Y", "EXP"])
    df = df.groupby("OR").agg({"X" : "mean", "Y" : "mean", "EXP" : "sum"})
    df = df[df["Y"] < 0]
    ORs = sorted(list(df.index))
    coords = df[["X", "Y"]].to_dict("index")
    scaler_y = MinMaxScaler().fit(df[["X", "Y"]].values)
    df = files
    n_col = len(df.columns)
    df.rename(columns={"observed": "OR"}, inplace=True)
    df = df[df["OR"].isin(ORs)]
    ORs_to_exclude = [x for x in ORs if len(df[df["OR"] == x]) < min_cells]
    df = df[~df["OR"].isin(ORs_to_exclude)]
    ORs = list(df.OR.unique())
    scaler_x = StandardScaler().fit(df.values[:,0:(n_col - 1)])
    X_train = df.values[:,0:(n_col - 1)]
    X_train = scaler_x.transform(X_train)
    train_ors = np.array(df.iloc[:,(n_col - 1)])
    counter = Counter(train_ors)
    oversample = SMOTE(random_state=random_state)
    X_train, train_ors = oversample.fit_resample(X_train, train_ors)
    counter = Counter(train_ors)
    max_cells = np.max(list(counter.values()))
    y_train = []
    for i in range(train_ors.shape[0]):
        y_train.append([coords[train_ors[i]]["X"], coords[train_ors[i]]["Y"]])
    y_train = np.array(y_train)
    y_train = scaler_y.transform(y_train)
    model.fit(X_train, y_train)
    # Make predictions on all other ORs
    df = files
    df.rename(columns={"observed": "OR"}, inplace=True)
    n_col = len(df.columns)
    all_ORs = list(df.OR.unique())
    ORs_of_interest = [x for x in all_ORs if len(df[df["OR"] == x]) >= min_cells]
    with open("data/lr_out.txt", "w") as f:
        for OR in ORs_of_interest:
            X = df[df["OR"] == OR].values[:,0:(n_col - 1)]
            n = X.shape[0]
            X = scaler_x.transform(X)
            pred = scaler_y.inverse_transform(model.predict(X))
            pred = np.mean(pred, axis=0)
            x_pred, y_pred = pred
            print(OR,n,"predicted",x_pred,y_pred, file=f, sep="\t")
            if OR in ORs:
                print(OR, n, "slide-seq", coords[OR]["X"], coords[OR]["Y"], file=f, sep="\t")
    out = pd.read_csv("data/lr_out.txt",sep = "\t", header=None)
    return out

