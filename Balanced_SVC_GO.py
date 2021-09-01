#!/usr/bin/env python
# coding: utf-8

# In[82]:


#!/usr/bin/python
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import accuracy_score, f1_score, make_scorer, balanced_accuracy_score
import numpy as np
from sklearn.svm import SVC
import sys
F = int(sys.argv[1])
R = int(sys.argv[2])


# In[83]:


data = pd.read_csv("SPCA_{}.csv".format(F))
data = data.reset_index(drop = True)
PCs = data.shape[1]-1
X = data.values[:,:PCs]
y = data["observed"].values
le = LabelEncoder()
y = le.fit_transform(y)


# In[146]:


balanced_score = []


# In[147]:


for i in range(R*100,R*100+5,1):
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=i)
    sss.get_n_splits(X, y)
    train_index, test_index = next(sss.split(X, y))
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    model = SVC(C = 0.01,kernel = "linear",class_weight =  "balanced", decision_function_shape = "ovo")
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    balanced_score.append(balanced_accuracy_score(y_test, y_pred))


# In[157]:




# In[158]:


results = pd.DataFrame({"balanced_score" : balanced_score})
results.to_csv('Set_{}_clf_{}_PCs_balanced_accuracy.csv'.format(F,PCs))


# In[ ]:




