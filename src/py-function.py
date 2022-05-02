# Load module
import numpy as np
import pandas as pd
import os
import datetime
from requests import get
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import RFECV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split, GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
from sklearn.linear_model import Lasso
from sklearn.impute import KNNImputer
from sklearn.svm import SVC


def feature_selection_svm_rfecv(X,Y):
  # model conf.
  k_fold = StratifiedKFold(n_splits=5, shuffle=True, random_state=331)
  svc = SVC(kernel = 'linear', probability = True)
  rfecv = RFECV(estimator=svc, cv=k_fold, scoring='roc_auc')
  
  ## search parameter
  param_grid = {
      # SVC - linear kernel의 경우 'C'만 parameter로 확인하면 됨
      'estimator__C': np.arange(0.01,100,0.01)
  }
  
  # CV conf.
  CV_rfc = RandomizedSearchCV(estimator=rfecv, 
                      param_distributions=param_grid, 
                      n_iter=50,
                      cv= k_fold, 
                      scoring = 'roc_auc', 
                      verbose=1,
                      random_state=331,
                      n_jobs=15)
  CV_rfc.fit(X, Y.values.ravel())
  
  return CV_rfc.best_estimator_.support_
  
  
def feature_selection_LASSO(X, Y):
  # ML pipeline
  pipeline = Pipeline([
                     ('scaler',StandardScaler()),
                     ('model',Lasso(max_iter=9999))])
                    
  # grid search
  search = GridSearchCV(pipeline,
                      {'model__alpha':np.arange(0.001,5,0.001)},
                      cv = 10, scoring="neg_mean_squared_error",verbose=1
                      )
                      
  search.fit(X,Y.values.ravel())

  print(search.best_params_)
  coefficients = search.best_estimator_.named_steps['model'].coef_
  importance = np.abs(coefficients)
  
  return importance



def load_tcga_dataset(pkl_path, raw_path, cancer_type):   
  # subfunction
  def non_zero_column(DF):
    sample_cnt = int(len(DF.columns) * 0.2)
    zero_row = dict(DF.isin([0]).sum(axis=1))
    non_remove_feature = list()

    for key, value in zero_row.items():
        if value < sample_cnt:
            non_remove_feature.append(key)
    
    return non_remove_feature
  
  def cancer_select(cols, cancer_type, raw_path):
      # phenotype
      phe1 = pd.read_csv(raw_path + "GDC-PANCAN.basic_phenotype.tsv", sep="\t")
      phe1 = phe1.loc[phe1.program == "TCGA", :].loc[:, ['sample', 'sample_type', 'project_id']].drop_duplicates(['sample'])
      phe1['sample'] =  phe1.apply(lambda x : x['sample'][:-1], axis=1)
      phe2 = pd.read_csv(raw_path + "TCGA_phenotype_denseDataOnlyDownload.tsv", sep="\t")
      ph_join = pd.merge(left = phe2 , right = phe1, how = "left", on = "sample").dropna(subset=['project_id'])
      
      if cancer_type == "PAN" or cancer_type == "PANCAN":
          filterd = ph_join.loc[ph_join['sample_type_y'] == "Primary Tumor", :]
          sample_barcode = filterd["sample"].tolist()
      else:
          filterd = ph_join.loc[((ph_join['sample_type_y'] == "Primary Tumor")|(ph_join['sample_type_y'] == "Solid Tissue Normal")) & (ph_join['project_id'] == "TCGA-" + cancer_type) , :]
          sample_barcode = filterd["sample"].tolist()
          
      intersect_ = list(set(cols).intersection(sample_barcode))
      
      return intersect_

# main
  if os.path.isfile(pkl_path + cancer_type + "_rna.pkl"):
    # sep
    rna = pd.read_pickle(pkl_path  + cancer_type + "_rna.pkl")
      
  else :
      # create dir
      Path(pkl_path).mkdir(parents=True, exist_ok=True)
      Path(raw_path).mkdir(parents=True, exist_ok=True)
      
      # file name
      mrna_f = 'tcga_RSEM_gene_fpkm.gz'
      
      # file url
      mrna_url = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_fpkm.gz"
       # to make GET request
      def download(url, file_name):
          with open(file_name, "wb") as file:   # open in binary mode
              response = get(url)               # get request
              file.write(response.content)      # write to file
  
      # mrna
      if os.path.isfile(raw_path + mrna_f) == False:
          download(mrna_url, raw_path + mrna_f)
      
      # RNA gene expression
      if os.path.isfile(pkl_path + cancer_type + "_rna.pkl") == False:
          col = pd.read_csv(raw_path + mrna_f, sep = "\t", index_col=0, nrows=0).columns.to_list()
          use_col = ['sample'] + cancer_select(cols=col, cancer_type=cancer_type, raw_path=raw_path)
          df_chunk = pd.read_csv(raw_path + mrna_f,
                       sep = "\t", index_col=0, iterator=True, chunksize=50000, usecols=use_col)
          rna = pd.concat([chunk for chunk in df_chunk])
          rna = rna[rna.index.isin(non_zero_column(rna))]
  
          rna.to_pickle(pkl_path + cancer_type + "_rna.pkl")
      else : 
          rna = pd.read_pickle(pkl_path  + cancer_type + "_rna.pkl")
      
      # set same column for merge
  
      # pickle save
      rna.to_pickle(pkl_path + "/" + cancer_type + "_rna.pkl")
  
  # set index
  rna_index = rna.index.to_list()
  
  # missing impute
  imputer = KNNImputer(n_neighbors=10)
  rna_impute = imputer.fit_transform(rna)
  
  omics = pd.DataFrame(rna_impute, columns=rna.columns)
  rna.index = rna_index
  
  return rna
