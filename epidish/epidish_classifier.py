import os
import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal
from matplotlib import pyplot as plt

from classifier_utils import *

# NOTE(milo): Make sure to change this once we figure out where sensitized belongs!
MARTINO2015_LABEL_MAP = {"nonallergic": 0, "allergic": 1, "sensitized": 2}
MARTINO2018_LABEL_MAP = {"control": 0, "allergic": 1, "resolved": 2}


def predict_bulk_dnam(M_control, M_control_var, M_disease, M_disease_var, cell_fracs,
                      cpg_subset=None):
  """
  Predict the bulk DNAm beta values that we would expect for a control and disease person. Also
  compute the expected variance of those predictions.

  M_control (pd.DataFrame) : Each row is a CpG and each col is a cell type.
  M_control_var (pd.DataFrame) : Each row is a CpG and each col is a cell type.
  M_disease (pd.DataFrame) : Each row is a CpG and each col is a cell type.
  M_disease_var (pd.DataFrame) : Each row is a CpG and each col is a cell type.
  cell_fracs (pd.DataFrame) : Each row is a patient and each col is a cell type.
  cpg_subset (list) : A list of CpG names that we care about predicting bulk values for (i.e ones
                      with significant p-value).
  """
  cell_fracs_np = cell_fracs.transpose().to_numpy()

  k, N = cell_fracs_np.shape

  # The control matrix has variables for the Intercept and each cell type.
  # Multiply (C x k+1) DNAm coeffs by augmented cell fractions (k+1 x N) to get a "bulk" DNA vector.
  M_control_np = M_control.to_numpy()
  M_control_var_np = M_control_var.to_numpy()

  # Add a leading vector of ones which will be multiplied against the intercept (effective adding it).
  intercept_and_cell_fracs_np = np.concatenate((np.ones((1, N)), cell_fracs_np), axis=0)

  # B_control = intercept + coeff * cell_frac
  B_control = np.dot(M_control_np, intercept_and_cell_fracs_np)
  B_control_var = np.dot(M_control_var_np, intercept_and_cell_fracs_np ** 2)

  # Multiply (C x k) DNAm delta coeffs by cell fractions (k x N) to get a "bulk" DNA delta vector.
  # The disease bulk vector = B_control + B_disease_delta.
  M_disease_np = M_disease.to_numpy()
  M_disease_var_np = M_disease_var.to_numpy()

  B_disease_delta = np.dot(M_disease_np, cell_fracs_np)
  B_disease_delta_var = np.dot(M_disease_var_np, cell_fracs_np ** 2)

  B_disease = B_control + B_disease_delta
  B_disease_var = B_control_var + B_disease_delta_var

  # Convert the numpy arrays back to pd.Dataframe. The rows are CpG locations and the columns are
  # the names of individuals.
  B_control = B_control.clip(min=0, max=1)
  B_disease = B_disease.clip(min=0, max=1)

  B_control = pd.DataFrame(B_control, index=M_control.index, columns=cell_fracs.index)
  B_control_var = pd.DataFrame(B_control_var, index=M_control.index, columns=cell_fracs.index)
  B_disease = pd.DataFrame(B_disease, index=M_control.index, columns=cell_fracs.index)
  B_disease_var = pd.DataFrame(B_disease_var, index=M_control.index, columns=cell_fracs.index)

  # Optionally take a subset of the CpG locations (rows) if we only care about some of them for
  # classification.
  if cpg_subset is not None:
    B_control = B_control.loc[cpg_subset,]
    B_control_var = B_control_var.loc[cpg_subset,]
    B_disease = B_disease.loc[cpg_subset,]
    B_disease_var = B_disease_var.loc[cpg_subset,]

  return B_control, B_control_var, B_disease, B_disease_var


def classify_patients(B_control, B_control_var, B_disease, B_disease_var,
                      observed_beta_or_mvalues, bulk=False):
  """
  Classify a person as control or disease using their observed bulk beta or M-values.

  B_control (pd.DataFrame) : Row for each CpG, col for each patient.
  B_control_var (pd.DataFrame) : Row for each CpG, col for each patient.
  B_disease (pd.DataFrame) : Row for each CpG, col for each patient.
  B_disease_var (pd.DataFrame) : Row for each CpG, col for each patient.
  observed_beta_or_mvalues (pd.DataFrame) : Row for each CpG, column for each patient.

  Returns a (dict) where keys are patient names, and values are the ratio of disease likelihood to
  control likelihood. Ratios > 1 indicate that someone is more likely to be disease than control.
  """
  patient_names = observed_beta_or_mvalues.columns
  signif_cpg_names = B_control.index
  likelihood_ratios = {}

  for i, patient in enumerate(patient_names):
    # For cell-specific data, we have a predicted bulk vector for each patient, since they all have
    # different cell fractions. For bulk data, the same predicted control and disease vector is used
    # for all patients.
    if bulk:
      pred_control_beta = B_control
      pred_control_var = B_control_var

      pred_disease_beta = B_disease
      pred_disease_var = B_disease_var

    else:
      pred_control_beta = B_control.loc[:,patient]
      pred_control_var = B_control_var.loc[:,patient]

      pred_disease_beta = B_disease.loc[:,patient]
      pred_disease_var = B_disease_var.loc[:,patient]

    # NOTE(milo): Way too slow.
    # control_ll = 0
    # disease_ll = 0

    # for cpg_name, j in enumerate(signif_cpg_names):
    #   obs_beta_or_m = observed_beta_or_mvalues.loc[signif_cpg_names,patient]

    #   control_ll += multivariate_normal.logpdf(
    #     obs_beta_or_m, mean=pred_control_beta[cpg_name], cov=np.sqrt(pred_control_var[cpg_name]))

    #   disease_ll += multivariate_normal.logpdf(
    #     obs_beta_or_m, mean=pred_disease_beta[cpg_name], cov=np.sqrt(pred_disease_var[cpg_name]))

    # NOTE(milo): Much faster, but doesn't work when there are thousands of CpG locations.
    # Doesn't make sense to allocate a massive matrix with only the diagonal filled.
    control_likelihood = multivariate_normal.pdf(
      observed_beta_or_mvalues.loc[signif_cpg_names,patient],
      mean=pred_control_beta,
      cov=np.diag(np.sqrt(pred_control_var)))

    disease_likelihood = multivariate_normal.pdf(
      observed_beta_or_mvalues.loc[signif_cpg_names,patient],
      mean=pred_disease_beta,
      cov=np.diag(np.sqrt(pred_disease_var)))

    likelihood_ratios[patient] = disease_likelihood / control_likelihood

    # likelihood_ratios[patient] = np.exp(disease_ll) / np.exp(control_ll)
    # print("Finished patient {}".format(i))

  return likelihood_ratios


def run_cell_specific_classifier(analysis_folder, cell_types, use_mvalues=True, p_value_thresh=0.05):
  print("=============== RUNNING *CELL SPECIFIC* CLASSIFIER ================")
  print(">> Parameters:")
  print("  >> analysis_folder:", analysis_folder)
  print("  >> use_mvalues:", use_mvalues)
  print("  >> p_value_thresh:", p_value_thresh)

  # STEP 1: Load in data from R output.
  coe_control, coe_change, cell_fracs, phenotypes, beta = \
    load_epidish_results(analysis_folder, use_mvalues)

  # STEP 2: Extract significant CpG locations (adjP < 0.05).
  signif_cpg = report_significant_cpgs(cell_types, coe_change, p_value_thresh=p_value_thresh)
  print("==> All significant CpG locations:")
  print(signif_cpg)

  # STEP 3: Make matrices with cell-specific methylation beta values.
  Mc, Mc_var, Md, Md_var = cell_methylation_matrices(coe_control, coe_change, cell_types)

  # STEP 4: Multiply methylation matrices by cell fractions to predict the class conditional bulk
  # beta values for each patient.
  B_control, B_control_var, B_disease, B_disease_var = \
      predict_bulk_dnam(Mc, Mc_var, Md, Md_var, cell_fracs, cpg_subset=signif_cpg)

  # STEP 5: Classify each person based on how well their measured bulk beta values matches either
  # the control or disease class-conditioned ones.
  likely_ratios = classify_patients(B_control, B_control_var, B_disease, B_disease_var, beta)

  return likely_ratios


def run_bulk_classifier(analysis_folder, use_mvalues=True, p_value_thresh=0.05):
  print("=============== RUNNING *BULK* CLASSIFIER ================")
  print(">> Parameters:")
  print("  >> analysis_folder:", analysis_folder)
  print("  >> use_mvalues:", use_mvalues)
  print("  >> p_value_thresh:", p_value_thresh)

  # STEP 1: Load in data from R output.
  coe_control, coe_change, cell_fracs, phenotypes, beta_or_mvalues = \
    load_epidish_results(analysis_folder, use_mvalues, has_cellfrac=False)

  assert(cell_fracs is None) # Expect these to be None.

  # STEP 2: Extract significant CpG locations (adjP < 0.05).
  signif_cpg = report_significant_cpgs_bulk(coe_change, p_value_thresh=p_value_thresh)
  print("==> All significant CpG locations:")
  print(signif_cpg)

  # STEP 3: From DMC linear regression, compute the mean and variance of the control and disease
  # bulk M-values.
  B_control, B_control_var, B_disease, B_disease_var = \
      bulk_control_and_disease_mean(coe_control, coe_change)

  B_control = B_control[signif_cpg]
  B_control_var = B_control_var[signif_cpg]
  B_disease = B_disease[signif_cpg]
  B_disease_var = B_disease_var[signif_cpg]

  # STEP 4: Classify each person based on how well their measured bulk beta values matches either
  # the control or disease class-conditioned ones.
  likely_ratios = classify_patients(B_control, B_control_var, B_disease, B_disease_var,
                                    beta_or_mvalues, bulk=True)

  return likely_ratios


def cs_2015_classifier(p_value_thresh, likelihood_ratio):
  """
  Extract significant CpG locations and classify patients from Martino 2015 using them.
  """
  analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_all/" 
  # analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_neutro/"
  # analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_eosino/"
  # analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_only_pbmc/"
  # NOTE(milo): Only use PBMC cell types here.
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono"]
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro"]
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Eosino"]
  cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]

  likely_ratios = run_cell_specific_classifier(
      analysis_folder, cell_types_m2015, use_mvalues=True, p_value_thresh=p_value_thresh)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "../phenotypes.csv"), index_col=0)

  with open(os.path.join(analysis_folder, "likelihood_ratios.txt"), "w") as f:
    for patient in likely_ratios:
      f.write("{} {}\n".format(patient, likely_ratios[patient]))

  labels = {}
  for patient in likely_ratios:
    predicted_label = 1 if likely_ratios[patient] > likelihood_ratio else 0
    pheno_str = str(phenotypes.loc[patient,"challenge outcome:ch1"])
    true_label = MARTINO2015_LABEL_MAP[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("\n====== CLASSIFICATION RESULTS =====")
  print("==> Precision:", precision)
  print("==> Recall:", recall)

  return precision, recall


def bulk_2015_classifier(p_value_thresh, likelihood_ratio):
  test_patients = list(pd.read_csv("../analysis/martino2015/test_set.txt", header=None)[0])
  analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_bulk/"

  likely_ratios = run_bulk_classifier(analysis_folder, use_mvalues=True, p_value_thresh=p_value_thresh)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "../phenotypes.csv"), index_col=0)

  with open(os.path.join(analysis_folder, "likelihood_ratios.txt"), "w") as f:
    for patient in likely_ratios:
      f.write("{} {}\n".format(patient, likely_ratios[patient]))

  labels = {}
  for patient in test_patients:
    predicted_label = 1 if likely_ratios[patient] > likelihood_ratio else 0
    pheno_str = str(phenotypes.loc[patient,"challenge outcome:ch1"])
    true_label = MARTINO2015_LABEL_MAP[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("\n====== CLASSIFICATION RESULTS =====")
  print("==> Precision:", precision)
  print("==> Recall:", recall)

  # print(labels)
  return precision, recall


def bulk_2018_classifier(p_value_thresh, likelihood_ratio):
  test_patients = list(pd.read_csv("../analysis/martino2018/test_set.txt", header=None)[0])
  analysis_folder = "../analysis/martino2018/Mvalues_control_vs_allergic_bulk/"

  likely_ratios = run_bulk_classifier(analysis_folder, use_mvalues=True, p_value_thresh=p_value_thresh)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "../phenotypes.csv"), index_col=0)

  with open(os.path.join(analysis_folder, "likelihood_ratios.txt"), "w") as f:
    for patient in likely_ratios:
      f.write("{} {}\n".format(patient, likely_ratios[patient]))

  labels = {}
  for patient in test_patients:
    predicted_label = 1 if likely_ratios[patient] > likelihood_ratio else 0
    pheno_str = str(phenotypes.loc[patient,"allergy status:ch1"])
    true_label = MARTINO2018_LABEL_MAP[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("\n====== CLASSIFICATION RESULTS =====")
  print("==> Precision:", precision)
  print("==> Recall:", recall)

  # print(labels)
  return precision, recall


def cs_2018_classifier(p_value_thresh, likelihood_ratio):
  test_patients = list(pd.read_csv("../analysis/martino2018/test_set.txt", header=None)[0])
  analysis_folder = "../analysis/martino2018/Mvalues_control_vs_allergic/"
  cell_types_m2018 = ["CD4T", "CD8T"]

  likely_ratios = run_cell_specific_classifier(
      analysis_folder, cell_types_m2018, use_mvalues=True, p_value_thresh=p_value_thresh)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "../phenotypes.csv"), index_col=0)

  with open(os.path.join(analysis_folder, "likelihood_ratios.txt"), "w") as f:
    for patient in likely_ratios:
      f.write("{} {}\n".format(patient, likely_ratios[patient]))

  labels = {}
  for patient in test_patients:
    predicted_label = 1 if likely_ratios[patient] > likelihood_ratio else 0


    pheno_str = str(phenotypes.loc[patient, "allergy status:ch1"])

    true_label = MARTINO2018_LABEL_MAP[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("\n====== CLASSIFICATION RESULTS =====")
  print("==> Precision:", precision)
  print("==> Recall:", recall)

  return precision, recall


def precision_recall_vs_pvalue_thresh(function_to_run):
  """
  Look at the classifier performance as we vary the p_value_threshold for choosing CpG features.
      
  Outputs a text file where each line has:
  p_value_thresh    precision    recall
  """
  likelihood_ratio = 1

  # for p_value_thresh in np.arange(0.05, 0.55, 0.05):
  for p_value_thresh in np.arange(0.05, 0.15, 0.01):
    pr, re = function_to_run(p_value_thresh, likelihood_ratio)
    
    with open("{}_pr_pvalue_thresh.txt".format(function_to_run.__name__), "a") as f:
      f.write("{} {} {}\n".format(p_value_thresh, pr, re))


def precision_recall_vs_cutoff(lr_files):
  """
  Look at classifier performance as we vary the likelihood ratio for classifying disease.

  Outputs a text file with format:
  likelihood_ratio    precision   recall
  """
  fig, axs = plt.subplots()

  for name in lr_files:
    precisions = []
    recalls = []

    fname, dataset = lr_files[name]
    patient_and_ratio = pd.read_csv(fname, header=None, index_col=0, delim_whitespace=True, names=["patient", "ratio"])

    test_set_fname = "../analysis/martino{}/test_set.txt".format(dataset)
    test_patients = pd.read_csv(test_set_fname, header=None)[0]

    phenotypes = pd.read_csv("../analysis/martino{}/phenotypes.csv".format(dataset), index_col=0)

    for lr in np.arange(0.5, 2.0, 0.05):
      labels = {}

      for patient in test_patients:
        predicted_label = 1 if patient_and_ratio.loc[patient, "ratio"] >= lr else 0

        pheno_str = str(phenotypes.loc[patient,"allergy status:ch1" if dataset == 2018 else "challenge outcome:ch1"])
        true_label = MARTINO2018_LABEL_MAP[pheno_str] if dataset == 2018 else MARTINO2015_LABEL_MAP[pheno_str]
        # print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
        labels[patient] = (predicted_label, true_label)

      pr, re = compute_precision_recall(labels)

      precisions.append(pr)
      recalls.append(re)

    axs.plot(recalls, precisions, label=name)

  axs.legend()
  axs.set_title("Classifier Precision and Recall")
  axs.set_xlabel("recall")
  axs.set_ylabel("precision")
  plt.show()


def make_param_precision_recall_plot(pr_files):
  """
  Each file has some precision-recall values for a classifier as we vary an input parameter.
  """
  fig, axs = plt.subplots()
  axs.set_xlabel("recall")
  axs.set_ylabel("precision")
  axs.set_title("Classifier Precision vs. Recall")

  for name in pr_files:
    fname = pr_files[name]
    results = np.loadtxt(fname)
    pr = results[:,1]
    re = results[:,2]
    axs.plot(re, pr, label=name, marker="o", linewidth=0, markersize=10)

  axs.legend()
  plt.show()


if __name__ == "__main__":  
  # cs_2015_classifier(0.05, 1) # Do for PBMC.
  # cs_2015_classifier(0.06, 1) # Do for ALL.

  # bulk_2018_classifier(0.02, 1)
  # cs_2018_classifier(0.25, 1)

  # precision_recall_vs_pvalue_thresh(bulk_2018_classifier)
  # precision_recall_vs_pvalue_thresh(cs_2018_classifier)

  # precision_recall_vs_pvalue_thresh(bulk_2015_classifier)
  # precision_recall_vs_pvalue_thresh(cs_2015_classifier)

  # pr_files = {
  #   "ALL_2015" : "cs_2015_classifier_pr_pvalue_thresh_all.txt",
  #   "PBMC_2015" : "cs_2015_classifier_pr_pvalue_thresh_pbmc.txt",
  #   "CD4_2018" : "bulk_2018_classifier_pr_pvalue_thresh.txt",
  #   "CD4_CD8_2018": "cs_2018_classifier_pr_pvalue_thresh.txt"
  # }
  # make_param_precision_recall_plot(pr_files)

  # Best p-value cutoffs:
  # 2015 ALL: 0.06
  # 2015 PBMC: 0.05
  # 2018 CD4: 0.02
  # 2018 CD4 and CD8: 0.25

  lr_files = {
    "2015_all": ("../analysis/likelihood_ratios/2015_all.txt", 2015),
    "2015_pbmc": ("../analysis/likelihood_ratios/2015_pbmc.txt", 2015),
    "2018_cd4": ("../analysis/likelihood_ratios/2018_cd4.txt", 2018),
    "2018_cd4_cd8": ("../analysis/likelihood_ratios/2018_cd4_cd8.txt", 2018)
  }
  precision_recall_vs_cutoff(lr_files)
