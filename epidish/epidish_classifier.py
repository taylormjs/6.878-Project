import os
import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal

from classifier_utils import *


def load_epidish_results(study="martino2015"):
  """
  Load R results from the analysis folder.
  """
  folder = "../analysis/{}/nonallergic_vs_allergic/".format(study)
  # NOTE(milo): index_col = 0 sets the first column as the row names.
  coe_control = pd.read_csv(os.path.join(folder, "{}_coe_control.csv".format(study)), index_col=0)
  coe_change = pd.read_csv(os.path.join(folder, "{}_coe_change.csv".format(study)), index_col=0)
  cell_frac = pd.read_csv(os.path.join(folder, "{}_cellfrac.csv".format(study)), index_col=0)
  pheno = pd.read_csv(os.path.join(folder, "{}_phenotypes.csv".format(study)), index_col=0)
  beta = pd.read_csv(os.path.join(folder, "{}_beta.csv".format(study)), index_col=0)

  return coe_control, coe_change, cell_frac, pheno, beta


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
  B_control_var = np.dot(M_control_var_np, intercept_and_cell_fracs_np)

  # Multiply (C x k) DNAm delta coeffs by cell fractions (k x N) to get a "bulk" DNA delta vector.
  # The disease bulk vector = B_control + B_disease_delta.
  M_disease_np = M_disease.to_numpy()
  M_disease_var_np = M_disease_var.to_numpy()

  B_disease_delta = np.dot(M_disease_np, cell_fracs_np)
  B_disease_delta_var = np.dot(M_disease_var_np, cell_fracs_np)

  B_disease = B_control + B_disease_delta
  B_disease_var = B_control_var + B_disease_delta_var

  # Convert the numpy arrays back to pd.Dataframe. The rows are CpG locations and the columns are
  # the names of individuals.
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


def classify_using_beta_values(B_control, B_control_var, B_disease, B_disease_var, observed_beta_values):
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

    pred_disease_beta = B_disease.loc[:,patient]
    pred_disease_var = B_disease_var.loc[:,patient]

    control_likelihood = multivariate_normal.pdf(
      observed_beta_values.loc[signif_cpg_names,patient],
      mean=pred_control_beta,
      cov=np.sqrt(pred_control_var))

    disease_likelihood = multivariate_normal.pdf(
      observed_beta_values.loc[signif_cpg_names,patient],
      mean=pred_disease_beta,
      cov=np.sqrt(pred_disease_var))

    likelihood_ratios[patient] = disease_likelihood / control_likelihood

  return likelihood_ratios


def run_martino2015():
  """
  Extract significant CpG locations and classify patients from Martino 2015 using them.
  """
  # STEP 1: Load in data from R output.
  coe_control, coe_change, cell_fracs, phenotypes, beta = load_epidish_results("martino2015")

  # NOTE(milo): Had to throw away Neutro in R because it was causing errors.
  cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]
  # coe_control = rename_control_cols(coe_control, cell_types_m2015)

  # STEP 2: Extract significant CpG locations (adjP < 0.05).
  signif_cpg = report_significant_cpgs(cell_types_m2015, coe_change)
  print("==> All significant CpG locations:")
  print(signif_cpg)

  # STEP 3: Make matrices with cell-specific methylation beta values.
  Mc, Mc_var, Md, Md_var = cell_methylation_matrices(coe_control, coe_change, cell_types_m2015)

  # STEP 4: Multiply methylation matrices by cell fractions to predict the class conditional bulk
  # beta values for each patient.
  B_control, B_control_var, B_disease, B_disease_var = \
      predict_bulk_dnam(Mc, Mc_var, Md, Md_var, cell_fracs, cpg_subset=signif_cpg)

  # # STEP 5: Classify each person based on how well their measured bulk beta values matches either
  # # the control or disease class-conditioned ones.
  # likely_ratios = classify_using_beta_values(B_control, B_control_var, B_disease, B_disease_var, beta)

  # # NOTE(milo): Make sure to change this once we figure out where sensitized belongs!
  # MARTINO2015_LABEL_MAP = {"nonallergic": 0, "allergic": 1,"sensitized": 2}

  # labels = {}
  # for patient in likely_ratios:
  #   predicted_label = 1 if likely_ratios[patient] > 1 else 0
  #   pheno_str = str(phenotypes.loc[patient,"challenge outcome:ch1"])
  #   true_label = MARTINO2015_LABEL_MAP[pheno_str]
  #   print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
  #   labels[patient] = (predicted_label, true_label)

  # precision, recall = compute_precision_recall(labels)
  # print("\n====== CLASSIFICATION RESULTS =====")
  # print("==> Precision:", precision)
  # print("==> Recall:", recall)


if __name__ == "__main__":
  run_martino2015()
