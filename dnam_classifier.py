import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal


def load_epidish_results_martino2015():
  # NOTE(milo): index_col = 0 sets the first column as the row names.
  coe_control = pd.read_csv("./analysis/martino2015_coe_control.csv", index_col=0)
  coe_change = pd.read_csv("./analysis/martino2015_coe_change.csv", index_col=0)
  cell_frac = pd.read_csv("./analysis/martino2015_cellfrac.csv", index_col=0)
  pheno = pd.read_csv("./analysis/martino2015_phenotypes.csv", index_col=0)
  beta = pd.read_csv("./analysis/martino2015_beta.csv", index_col=0)
  return coe_control, coe_change, cell_frac, pheno, beta


def cell_methylation_matrices(coe_control, coe_change, cell_types):
  """
  Build a two matrices with a row for each CpG location, and a column for each cell type. Each entry
  M_ij in the matrices represents an average beta value for CpG location i and cell type j.
  """
  assert(coe_control.shape[0] == coe_change.shape[0])

  num_cpg_locations = coe_control.shape[0]
  num_cell_types = len(cell_types)

  estimate_colnames = ["{}.Estimate".format(cellname) for cellname in cell_types]
  se_colnames = ["{}.SE".format(cellname) for cellname in cell_types]

  # Make the control matrix.
  M_control = coe_control[estimate_colnames]
  M_control_var = coe_control[se_colnames] ** 2
  assert(M_control.shape == (num_cpg_locations, num_cell_types))
  assert(M_control_var.shape == (num_cpg_locations, num_cell_types))

  # Make the disease matrix (control + change).
  M_disease = M_control + coe_change[estimate_colnames]
  M_disease_var = M_control_var + coe_change[se_colnames] ** 2
  assert(M_disease.shape == (num_cpg_locations, num_cell_types))
  assert(M_disease_var.shape == (num_cpg_locations, num_cell_types))

  return M_control, M_control_var, M_disease, M_disease_var


def report_significant_cpgs(cell_types, coe_change, p_value_thresh=0.05):
  """
  Print out the number of significant CpG locations that were found for each cell type.
  """
  signif_set = set()

  for cell_type in cell_types:
    adjP_colname = "{}.adjP".format(cell_type)
    p_values = coe_change[coe_change[adjP_colname] <= p_value_thresh]
    num_signif = p_values.shape[0]
    print("Found {} significant CpGs for {}".format(cell_type, num_signif))

    for cpg_name in p_values.index:
      signif_set.add(cpg_name)
  
  return sorted(list(signif_set))


def rename_control_cols(coe_control, cell_types):
  """
  Index(['Unnamed: 0', 'Estimate', 'SE', 't', 'p', 'adjP', 'Estimate.1', 'SE.1',
       't.1', 'p.1', 'adjP.1', 'Estimate.2', 'SE.2', 't.2', 'p.2', 'adjP.2',
       'Estimate.3', 'SE.3', 't.3', 'p.3', 'adjP.3', 'Estimate.4', 'SE.4',
       't.4', 'p.4', 'adjP.4', 'Estimate.5', 'SE.5', 't.5', 'p.5', 'adjP.5'],
      dtype='object')
  """
  colname_map = {}
  for i in range(len(cell_types)):
    if i == 0:
      colname_map["Estimate"] = "{}.Estimate".format(cell_types[i])
      colname_map["SE"] = "{}.SE".format(cell_types[i])
      colname_map["t"] = "{}.t".format(cell_types[i])
      colname_map["p"] = "{}.p".format(cell_types[i])
      colname_map["adjP"] = "{}.adjP".format(cell_types[i])
    else:
      colname_map["Estimate.{}".format(i)] = "{}.Estimate".format(cell_types[i])
      colname_map["SE.{}".format(i)] = "{}.SE".format(cell_types[i])
      colname_map["t.{}".format(i)] = "{}.t".format(cell_types[i])
      colname_map["p.{}".format(i)] = "{}.p".format(cell_types[i])
      colname_map["adjP.{}".format(i)] = "{}.adjP".format(cell_types[i])

  coe_control = coe_control.rename(columns=colname_map)
  return coe_control


def predict_bulk_dnam(M_control, M_control_var, M_disease, M_disease_var, cell_fracs, cpg_subset=None):
  """
  Predict the bulk DNAm beta values that we would expect for a control and disease person. Also
  compute the expected variance of those predictions.
  """
  cell_fracs_np = cell_fracs.transpose().to_numpy()

  # Multiply DNAm coeff (C x k) by cell fractions (k x N) to get a "bulk" DNAm vector for each
  # person (C x N). Each column represents the bulk DNAm beta vector we would expect for a control
  # or disease person.
  M_control_np = M_control.to_numpy()
  M_control_var_np = M_control_var.to_numpy()
  B_control = np.dot(M_control_np, cell_fracs_np)
  B_control_var = np.dot(M_control_var_np, cell_fracs_np)

  M_disease_np = M_disease.to_numpy()
  M_disease_var_np = M_disease_var.to_numpy()
  B_disease = np.dot(M_disease_np, cell_fracs_np)
  B_disease_var = np.dot(M_disease_var_np, cell_fracs_np)

  # Convert the numpy arrays back to pd.Dataframe. The rows are CpG locations and the columns are
  # the names of individuals.
  B_control = pd.DataFrame(B_control, index=M_control.index, columns=cell_fracs.index)
  B_control_var = pd.DataFrame(B_control_var, index=M_control.index, columns=cell_fracs.index)
  B_disease = pd.DataFrame(B_disease, index=M_control.index, columns=cell_fracs.index)
  B_disease_var = pd.DataFrame(B_disease_var, index=M_control.index, columns=cell_fracs.index)

  if cpg_subset is not None:
    B_control = B_control.loc[cpg_subset,]
    B_control_var = B_control_var.loc[cpg_subset,]
    B_disease = B_disease.loc[cpg_subset,]
    B_disease_var = B_disease_var.loc[cpg_subset,]

  return B_control, B_control_var, B_disease, B_disease_var


def classify_using_beta_values(B_control, B_control_var, B_disease, B_disease_var, observed_beta_values):
  """
  Classify a person as control or disease using their observed bulk beta values.

  observed_beta_values (pd.Dataframe) : Row for each CpG, column for each patient.
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


def compute_precision_recall(labels):
  """
  labels (dict) : Keys are patients, values are (pred_label, true_label).
  """
  tp = 0
  fp = 0
  tn = 0
  fn = 0

  for l in labels:
    pred, true = labels[l]

    if pred == true:
      if pred == 0:
        tn += 1
      elif pred == 1:
        tp += 1

    if pred != true:
      if pred == 0:
        fn += 1
      elif pred == 1:
        fp += 1

  precision = tp / (tp + fp)
  recall = tp / (tp + fn)

  return precision, recall


def run_martino2015():
  coe_control, coe_change, cell_fracs, phenotypes, beta = load_epidish_results_martino2015()

  # NOTE(milo): Had to throw away Neutro in R because it was causing errors.
  cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Eosino"]
  coe_control = rename_control_cols(coe_control, cell_types_m2015)

  signif_cpg = report_significant_cpgs(cell_types_m2015, coe_change)
  print("==> All significant CpG locations:")
  print(signif_cpg)

  Mc, Mc_var, Md, Md_var = cell_methylation_matrices(coe_control, coe_change, cell_types_m2015)

  # print("============= Control Methylation by Cell Type ==============")
  # print(Mc.head())
  # print(Mc_stdev.head())

  # print("============= Disease Methylation by Cell Type ==============")
  # print(Md.head())
  # print(Md_stdev.head())

  # print(cell_frac.head())
  # print(phenotypes.head())

  B_control, B_control_var, B_disease, B_disease_var = \
      predict_bulk_dnam(Mc, Mc_var, Md, Md_var, cell_fracs, cpg_subset=signif_cpg)
  # print("B_control:", B_control.shape)
  # print("B_control_var:", B_control_var.shape)

  # print("==> Bulk predictions (CONTROL):")
  # print(B_control)
  # print(B_control_var)

  # print("==> Bulk predictions (DISEASE):")
  # print(B_disease)
  # print(B_disease_var)

  likely_ratios = classify_using_beta_values(B_control, B_control_var, B_disease, B_disease_var, beta)

  label_map = {
    "nonallergic": 0,
    "allergic": 1,
    "sensitized": 2
  }

  labels = {}
  for patient in likely_ratios:
    predicted_label = 1 if likely_ratios[patient] > 1 else 0
    pheno_str = str(phenotypes.loc[patient,"challenge outcome:ch1"])
    true_label = label_map[pheno_str]
    print("Patient {}: predicted={} true={}".format(patient, predicted_label, true_label))
    labels[patient] = (predicted_label, true_label)

  precision, recall = compute_precision_recall(labels)
  print("==> Precision=", precision)
  print("==> Recall=", recall)


if __name__ == "__main__":
  run_martino2015()
