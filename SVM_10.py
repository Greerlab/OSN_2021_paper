import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from imblearn.under_sampling import RandomUnderSampler 
import numpy as np
from sklearn.metrics import balanced_accuracy_score, accuracy_score
import argparse

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('-i','--input', help='input name', required=True, type=str)
args = parser.parse_args()
ipt = args.input

train_file_name = "train_"+ipt+".csv"
test_file_name = "test_"+ipt+".csv"

train = pd.read_csv(train_file_name, index_col = 0)
train = train.reset_index(drop = True)
X_train = train.values[:,:8]
y_train = train["OR"].values

test = pd.read_csv(test_file_name, index_col = 0)
X_test = test.values[:,:8]
y_test = test["OR"].values

le = LabelEncoder()
y_train = le.fit_transform(y_train)

undersample = RandomUnderSampler(sampling_strategy='majority')
X_train, y_train = undersample.fit_resample(X_train, y_train)

model = SVC(C = 0.01, class_weight = 'balanced', decision_function_shape = 'ovo', kernel = 'linear')
model.fit(X_train, y_train)

y_pred = model.predict(X_test)
y_pred = le.inverse_transform(y_pred)
print(balanced_accuracy_score(y_test, y_pred))

dicts = {'barcode': test.index, 'observed': y_test, 'predicted': y_pred} 
df = pd.DataFrame(dicts)
df.to_csv("Result_"+ipt+".csv")
