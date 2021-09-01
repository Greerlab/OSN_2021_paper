#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
from sklearn.metrics import balanced_accuracy_score
def balanced_score(test, pred):
    return balanced_accuracy_score(test,pred)

