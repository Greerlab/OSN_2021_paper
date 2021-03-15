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
R = int(sys.argv[1])


# In[83]:


data = pd.read_csv('svm_654_rm_all_50SPCs_3000G.csv')
data = data.reset_index(drop = True)
X = data.values[:,:42]
y = data["observed"].values
le = LabelEncoder()
y = le.fit_transform(y)
z = data["barcodes"].values


# In[146]:


test_barcode = []
test_OR = []
pred_OR = []
balanced_score = []
rep = []


# In[147]:


for i in range(R*100,R*100+10,1):
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=i)
    sss.get_n_splits(X, y)
    train_index, test_index = next(sss.split(X, y))
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    z_test = z[test_index]
    model = SVC(C = 0.01,kernel = "linear",class_weight =  "balanced", decision_function_shape = "ovo")
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    balanced_score.append(balanced_accuracy_score(y_test, y_pred))
    test_barcode.append(z_test)
    test_OR.append(le.inverse_transform(y_test))
    pred_OR.append(le.inverse_transform(y_pred))
    rep.append(np.repeat(i, len(y_test), axis=0))


# In[157]:


results = pd.DataFrame({"barcodes" : np.concatenate(test_barcode),
                       "test" : np.concatenate(test_OR),
                       "pred" : np.concatenate(pred_OR),
                       "rep" : np.concatenate(rep)})
results.to_csv('clf_{}_results.csv'.format(R))


# In[158]:


results = pd.DataFrame({"balanced_score" : balanced_score})
results.to_csv('clf_{}_balanced_accuracy.csv'.format(R))


# In[ ]:




