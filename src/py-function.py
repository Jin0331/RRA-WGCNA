# Load module
# module
import numpy as np
from numpy import interp
import pandas as pd
import os
import datetime
from itertools import cycle
from requests import get
from pathlib import Path
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.feature_selection import RFECV
from sklearn.model_selection import train_test_split, GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold, RandomizedSearchCV, KFold, cross_val_score
from sklearn.linear_model import Lasso
from sklearn.impute import KNNImputer
from sklearn.svm import SVC,SVR
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, auc
import matplotlib.pyplot as plt
from scipy import interp
from imblearn.over_sampling import SMOTE

# roc curve
def roc_auc_function(y_test, y_score, num_class):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(len(num_class)): 
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i]) 
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area 
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel()) 
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(len(num_class))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(len(num_class)):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    # Finally average it and compute AUC
    mean_tpr /= len(num_class)
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    return fpr, tpr, roc_auc
def roc_acu_calculator(DF, feature_name, log_save, over_sampling):
  X = DF.iloc[:, 1:]
  y = DF.iloc[:, 0]
  
  # SMOTE oversampling for minority
  if over_sampling:
    sm = SMOTE("minority", random_state=331)
    X, y = sm.fit_resample(X,y)
    
  # multi-class detection
  num_class = set(y)
  if len(num_class) > 2:
    y = label_binarize(y, classes=list(num_class))
  
  # Learn to predict each class against the other 
  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = .3, random_state=331)
  classifier = OneVsRestClassifier(SVC(kernel='linear', probability=True, random_state=331)) 
  y_score = classifier.fit(X_train, y_train).decision_function(X_test)
  
  if len(num_class) > 2:    
      
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(len(num_class)): 
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i]) 
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area 
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel()) 
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(len(num_class))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(len(num_class)):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    # Finally average it and compute AUC
    mean_tpr /= len(num_class)
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    lw = 2
    plt.figure(figsize=(10,10))
    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.4f})'
                   ''.format(roc_auc["micro"]),
             color='deeppink', linestyle=':', linewidth=4)

    plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.4f})'
                   ''.format(roc_auc["macro"]),
             color='navy', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'darkgreen', 'violet', 'peru', 'gold'])
    for i, color in zip(range(len(num_class)), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                 label='ROC curve of class {0} (AUC={1:0.4f})'
                 ''.format(i, roc_auc[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(feature_name + " (Multi-class)")
    plt.legend(loc="lower right")
    plt.savefig(log_save + "/" + feature_name + '_ROC_AUC.png')
    plt.clf()
    
    # calculation
    y_prob = classifier.predict_proba(X_test)
    macro_roc_auc_ovr = roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro")
    micro_roc_auc_ovr = roc_auc_score(
        y_test, y_prob, multi_class="ovr", average="micro")
    return {'macro' : macro_roc_auc_ovr, 'micro': micro_roc_auc_ovr}
    return macro_roc_auc_ovr
  else :
    fpr, tpr, thres = roc_curve(y_test, y_score)
    roc_auc = roc_auc_score(y_test, y_score)

    # roc curve
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.4f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(feature_name + " (Binary-class")
    plt.legend(loc="best")
    plt.savefig(log_save + "/" + feature_name + '_ROC_AUC.png')
    plt.clf()
    
    return roc_auc
def roc_acu_calculator_cv(DF, feature_name, log_save, over_sampling):
    X = DF.iloc[:, 1:]
    y = DF.iloc[:, 0]

    # SMOTE oversampling for minority
    if over_sampling:
        sm = SMOTE(sampling_strategy="minority", random_state=331)
        X, y = sm.fit_resample(X,y)

    # multi-class detection
    num_class = set(y)
    if len(num_class) > 2:
        y = label_binarize(y, classes=list(num_class))

    fold_cnt = 5
    kfold = KFold(n_splits=fold_cnt)
    scores = []
    pipe_line = make_pipeline(StandardScaler(),
                             OneVsRestClassifier(SVC(kernel='linear', probability=True, random_state=331)))

    roc_auc_fold = dict()
    fpr_fold = dict()
    tpr_fold = dict()

    for k, (train, test) in enumerate(kfold.split(DF)):
        y_score = pipe_line.fit(X.iloc[train, :], y[train]).decision_function(X.iloc[test, :])

        if len(num_class) > 2:
            fpr, tpr, roc_auc = roc_auc_function(y_test=y[test], y_score=y_score, num_class=num_class)
            roc_auc_fold[k] = roc_auc['micro']
            fpr_fold[k] = fpr['micro']
            tpr_fold[k] = tpr['micro']

        else :
            fpr, tpr, thres = roc_curve(y[test], y_score, drop_intermediate=False)
            roc_auc_fold[k] = roc_auc_score(y[test], y_score)
            fpr_fold[k] = fpr
            tpr_fold[k] = tpr


    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr_fold[i] for i in range(len(num_class))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(len(num_class)):
        mean_tpr += interp(all_fpr, fpr_fold[i], tpr_fold[i])

    # Finally average it and compute AUC
    mean_tpr /= len(num_class)
    fpr_fold["macro"] = all_fpr
    tpr_fold["macro"] = mean_tpr
    roc_auc_fold["macro"] = auc(fpr_fold["macro"], tpr_fold["macro"])

    # Plot all ROC curves
    lw = 2
    plt.figure(figsize=(10,10))
    plt.plot(fpr_fold["macro"], tpr_fold["macro"],
             label='macro-average ROC curve (AUC = {0:0.3f})'
                   ''.format(roc_auc_fold["macro"]),
             color='red', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'darkgreen', 'violet', 'peru', 'gold'])
    for i, color in zip(range(len(num_class)), colors):
        plt.plot(fpr_fold[i], tpr_fold[i], color=color, lw=lw, alpha=0.2,
                 label='ROC curve fold-{0} (AUC={1:0.3f})'
                 ''.format(i + 1, roc_auc_fold[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(feature_name + "_Multi-class-CV")
    plt.legend(loc="lower right")
    plt.savefig(log_save + "/" + feature_name + '_ROC_AUC_CV.png')
    plt.close()

    return roc_auc_fold["macro"]
  
# gene selection
def feature_selection_svm_rfecv(X, Y, over_sampling):
  # SMOTE oversampling for minority
  if over_sampling:
    sm = SMOTE("minority",random_state=331)
    X, y = sm.fit_resample(X,y)  

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
                      n_jobs=10)
  CV_rfc.fit(X, Y.values.ravel())
  
  return CV_rfc.best_estimator_.support_
def feature_selection_LASSO(X, y, over_sampling=None):
    num_class = set(y)
    
    # SMOTE oversampling for minority
    if over_sampling:
      sm = SMOTE("minority",random_state=331)
      X, y = sm.fit_resample(X,y)
    
    if len(num_class) > 2:
        y = label_binarize(y, classes=list(num_class))

        # ML pipeline
        pipeline = Pipeline([
                     ('scaler',StandardScaler()),
                     ('model', OneVsRestClassifier(Lasso(max_iter=99999)))])

        # grid search
        search = GridSearchCV(pipeline,
                      {'model__estimator__alpha':np.arange(0.01,3,0.001)},
                      cv = 10, scoring="neg_mean_squared_error",verbose=1
                      )

        search.fit(X,y)

        print(search.best_params_)
        coefficients = pd.DataFrame(search.best_estimator_.named_steps.model.coef_)
        return coefficients.abs().sum().to_list()

    else :
        pipeline = Pipeline([
                     ('scaler',StandardScaler()),
                     ('model',Lasso(max_iter=99999))])

        # grid search
        search = GridSearchCV(pipeline,
                          {'model__alpha':np.arange(0.01,3,0.001)},
                          cv = 10, scoring="neg_mean_squared_error",verbose=1
                          )

        search.fit(X,y.values.ravel())

        print(search.best_params_)
        coefficients = search.best_estimator_.named_steps['model'].coef_
        return coefficients
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
def load_tcga_dataset(pkl_path, raw_path, cancer_type):   
  # subfunction
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
