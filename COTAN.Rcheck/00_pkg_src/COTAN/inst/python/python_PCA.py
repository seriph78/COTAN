import pandas as pd
import sklearn
from sklearn.decomposition import PCA

def python_PCA(file):
  f = PCA(n_components=10)
  f_principalComponents = f.fit_transform(file)
  return f_principalComponents

