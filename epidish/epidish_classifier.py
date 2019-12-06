import os
import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal

from classifier_utils import *

# NOTE(milo): Make sure to change this once we figure out where sensitized belongs!
MARTINO2015_LABEL_MAP = {"nonallergic": 0, "allergic": 1,"sensitized": 2}
MARTINO2018_LABEL_MAP = {"control": 0, "allergic": 1,"resolved": 2}


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


# def compute_bulk_dnam_variance(B_control, B_disease, observed_beta_values, phenotypes):
#   signif_cpg = B_control.index

#   control_patients = phenotypes[phenotypes["challenge outcome:ch1"] == "nonallergic"].index
#   disease_patients = phenotypes[phenotypes["challenge outcome:ch1"] == "allergic"].index

#   print("Fitting variance using {} control and {} disease patients".format(
#     len(control_patients), len(disease_patients)))

#   sq_control_resid = (B_control.loc[:,control_patients].to_numpy() - \
#                       observed_beta_values.loc[signif_cpg,control_patients].to_numpy()) ** 2
#   control_var = np.sum(sq_control_resid, axis=1) / (control_patients.shape[0] - 1)

#   sq_disease_red = (B_disease.loc[:,disease_patients].to_numpy() - \
#                     observed_beta_values.loc[signif_cpg,disease_patients].to_numpy()) ** 2
#   disease_var = np.sum(sq_control_resid, axis=1) / (disease_patients.shape[0] - 1)

#   control_var = pd.DataFrame(control_var, index=B_control.index, columns=["control_variance"])
#   disease_var = pd.DataFrame(disease_var, index=B_disease.index, columns=["disease_variance"])

#   return control_var, disease_var


def classify_using_beta_values(B_control, B_control_var, B_disease, B_disease_var,
                               observed_beta_values):
  """
  Classify a person as control or disease using their observed bulk beta values.

  B_control (pd.DataFrame) : Row for each CpG, col for each patient.
  B_control_var (pd.DataFrame) : Row for each CpG, col for each patient.
  B_disease (pd.DataFrame) : Row for each CpG, col for each patient.
  B_disease_var (pd.DataFrame) : Row for each CpG, col for each patient.
  observed_beta_values (pd.DataFrame) : Row for each CpG, column for each patient.

  Returns a (dict) where keys are patient names, and values are the ratio of disease likelihood to
  control likelihood. Ratios > 1 indicate that someone is more likely to be disease than control.
  """
  patient_names = observed_beta_values.columns
  signif_cpg_names = B_control.index
  likelihood_ratios = {}

  for i, patient in enumerate(patient_names):
    pred_control_beta = B_control.loc[:,patient]
    pred_control_var = B_control_var.loc[:,patient]
    # pred_control_var = np.squeeze(B_control_var)

    pred_disease_beta = B_disease.loc[:,patient]
    pred_disease_var = B_disease_var.loc[:,patient]
    # pred_disease_var = np.squeeze(B_disease_var)

    control_likelihood = multivariate_normal.pdf(
      observed_beta_values.loc[signif_cpg_names,patient],
      mean=pred_control_beta,
      cov=np.diag(np.sqrt(pred_control_var)))

    disease_likelihood = multivariate_normal.pdf(
      observed_beta_values.loc[signif_cpg_names,patient],
      mean=pred_disease_beta,
      cov=np.diag(np.sqrt(pred_disease_var)))

    likelihood_ratios[patient] = disease_likelihood / control_likelihood

  return likelihood_ratios


def run_classifier(analysis_folder, cell_types, use_mvalues=False, p_value_thresh=0.05):
  print("=============== RUNNING CLASSIFIER ================")
  print(">> Parameters:")
  print("  >> use_mvalues:", use_mvalues)
  print("  >> p_value_thresh:", p_value_thresh)

  # STEP 1: Load in data from R output.
  coe_control, coe_change, cell_fracs, phenotypes, beta = \
    load_epidish_results(analysis_folder, use_mvalues)

  # coe_control = rename_control_cols(coe_control, cell_types_m2015)

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
  likely_ratios = classify_using_beta_values(B_control, B_control_var, B_disease, B_disease_var, beta)

  return likely_ratios


def classify_martino2015():
  """
  Extract significant CpG locations and classify patients from Martino 2015 using them.
  """
  analysis_folder = "../analysis/martino2015/nonallergic_vs_allergic_with_eosino/" 
  # NOTE(milo): Only use PBMC cell types here.
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono"]
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro"]
  cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Eosino"]
  likely_ratios = run_classifier(analysis_folder, cell_types_m2015)

  print(likely_ratios)

  labels = {}
  for patient in likely_ratios:
    predicted_label = 1 if likely_ratios[patient] > 1 else 0
    pheno_str = str(phenotypes.loc[patient,"challenge outcome:ch1"])
    true_label = MARTINO2015_LABEL_MAP[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("\n====== CLASSIFICATION RESULTS =====")
  print("==> Precision:", precision)
  print("==> Recall:", recall)


def classify_martino2018():
  """
  Extract significant CpG locations and classify patients from Martino 2015 using them.
  """
  analysis_folder = "../analysis/martino2018/control_vs_allergic/"
  cell_types_m2018 = ["CD4T", "CD8T"]
  likely_ratios = run_classifier(analysis_folder, cell_types_m2018)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "phenotypes.csv"), index_col=0)

  labels = {}
  for patient in likely_ratios:
    predicted_label = 1 if likely_ratios[patient] > 1 else 0
    pheno_str = str(phenotypes.loc[patient,"allergy status:ch1"])
    true_label = MARTINO2018_LABEL_MAP[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("\n====== CLASSIFICATION RESULTS =====")
  print("==> Precision:", precision)
  print("==> Recall:", recall)


def classify_martino2015_Mvalues():
  """
  Extract significant CpG locations and classify patients from Martino 2015 using them.
  """
  # analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_all/" 
  # analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_neutro/"
  analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_eosino/"
  # NOTE(milo): Only use PBMC cell types here.
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono"]
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro"]
  cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Eosino"]
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]
  likely_ratios = run_classifier(analysis_folder, cell_types_m2015, use_mvalues=True, p_value_thresh=0.05)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "phenotypes.csv"), index_col=0)

  labels = {}
  for patient in likely_ratios:
    predicted_label = 1 if likely_ratios[patient] > 1 else 0
    pheno_str = str(phenotypes.loc[patient,"challenge outcome:ch1"])
    true_label = MARTINO2015_LABEL_MAP[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("\n====== CLASSIFICATION RESULTS =====")
  print("==> Precision:", precision)
  print("==> Recall:", recall)


def classify_martino2018_Mvalues():
  analysis_folder = "../analysis/martino2018/Mvalues_control_vs_allergic/"
  cell_types_m2018 = ["CD4T", "CD8T"]
  likely_ratios = run_classifier(analysis_folder, cell_types_m2018, use_mvalues=True, p_value_thresh=0.1)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "phenotypes.csv"), index_col=0)

  labels = {}
  for patient in likely_ratios:
    predicted_label = 1 if likely_ratios[patient] > 1 else 0
    pheno_str = str(phenotypes.loc[patient,"allergy status:ch1"])
    true_label = MARTINO2018_LABEL_MAP[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("\n====== CLASSIFICATION RESULTS =====")
  print("==> Precision:", precision)
  print("==> Recall:", recall)


if __name__ == "__main__":
  # classify_martino2015()
  # classify_martino2018()
  classify_martino2015_Mvalues()
  # classify_martino2018_Mvalues()
