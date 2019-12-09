import os
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal


def load_epidish_results(folder, mvalues=True, has_cellfrac=True):
  """
  Load R results from the analysis folder.
  """
  # NOTE(milo): index_col = 0 sets the first column as the row names.
  coe_control = pd.read_csv(os.path.join(folder, "coe_control.csv"), index_col=0)
  coe_change = pd.read_csv(os.path.join(folder, "coe_change.csv"), index_col=0)

  if has_cellfrac:
    cell_frac = pd.read_csv(os.path.join(folder, "cellfrac.csv"), index_col=0)
  else:
    cell_frac = None
  
  pheno = pd.read_csv(os.path.join(folder, "../phenotypes.csv"), index_col=0)

  if not mvalues:
    beta = pd.read_csv(os.path.join(folder, "../beta.csv"), index_col=0)
    return coe_control, coe_change, cell_frac, pheno, beta
  else:
    M = pd.read_csv(os.path.join(folder, "../Mvalues.csv"), index_col=0)
    return coe_control, coe_change, cell_frac, pheno, M


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
      else:
        raise NotImplementedError()

    if pred != true:
      if pred == 0:
        fn += 1
      elif pred == 1:
        fp += 1
      else:
        raise NotImplementedError()

  total = tp + fp + tn + fn
  assert(total == len(labels))

  precision = tp / (tp + fp) if (tp + fp) > 0 else 1.0
  recall = tp / (tp + fn) if (tp + fn) > 0 else 1.0

  return precision, recall


def report_significant_cpgs(cell_types, coe_change, p_value_thresh=0.05):
  """
  Print out the number of significant CpG locations that were found for each cell type.
  """
  signif_set = set()

  for cell_type in cell_types:
    adjP_colname = "{}.adjP".format(cell_type)
    p_values = coe_change[coe_change[adjP_colname] <= p_value_thresh]
    num_signif = p_values.shape[0]
    print("Found {} significant CpGs for {}".format(num_signif, cell_type))

    for cpg_name in p_values.index:
      signif_set.add(cpg_name)
  
  return sorted(list(signif_set))


def report_significant_cpgs_bulk(coe_change, p_value_thresh=0.05):
  """
  Print out the number of significant CpG locations found from bulk data.

  Returns: A list of CpG locations, sorted by their cg code.
  """
  signif_set = set()

  adjP_colname = "adjP"
  p_values = coe_change[coe_change[adjP_colname] <= p_value_thresh]
  num_signif = p_values.shape[0]
  print("Found {} significant CpGs".format(num_signif))

  for cpg_name in p_values.index:
    signif_set.add(cpg_name)
  
  return sorted(list(signif_set))


def report_cell_specific_DMCs(cell_types, coe_change, p_value_thresh=0.05):
  cpgs = defaultdict(lambda: [])

  dmcs = pd.DataFrame(columns=cell_types)

  for cell_type in cell_types:
    adjP_colname = "{}.adjP".format(cell_type)
    p_values = coe_change[coe_change[adjP_colname] <= p_value_thresh]

    num_signif = p_values.shape[0]
    print("Found {} significant CpGs for {}".format(num_signif, cell_type))

    for cpg_name in p_values.index:
      coeff = coe_change[coe]
      cpgs[cpg_name].append(cell_type)

  return cpgs


def cell_methylation_matrices(coe_control, coe_change, cell_types):
  """
  Build two matrices with a row for each CpG location, and a column for each cell type. Each entry
  M_ij in the matrices represents an average beta value for CpG location i and cell type j.

  coe_control (pd.DataFrame) : Has a row for each CpG and a columns for the linear regression
                               results for the Intercept and each cell type.
  coe_change (pd.DataFrame) : Has a row for each CpG and a columns for the linear regression
                              results for each cell type.
  cell_types (list of str) : A list of cell type names that we care about.
  """
  assert(coe_control.shape[0] == coe_change.shape[0])

  num_cpg_locations = coe_control.shape[0]
  num_cell_types = len(cell_types)

  # Make the control matrix.
  c_est_colnames = ["frac.m{}.Estimate".format(cellname) for cellname in cell_types]
  c_se_colnames = ["frac.m{}.SE".format(cellname) for cellname in cell_types]
  c_est_colnames.insert(0, "(Intercept).Estimate")
  c_se_colnames.insert(0, "(Intercept).SE")

  M_control = coe_control[c_est_colnames]
  M_control_var = coe_control[c_se_colnames] ** 2

  # NOTE(milo): We include an (Intercept) column, so there is a +1.
  assert(M_control.shape == (num_cpg_locations, 1 + num_cell_types))
  assert(M_control_var.shape == (num_cpg_locations, 1 + num_cell_types))

  # Make the disease matrix (change).
  estimate_colnames = ["{}.Estimate".format(cellname) for cellname in cell_types]
  se_colnames = ["{}.SE".format(cellname) for cellname in cell_types]
  M_disease = coe_change[estimate_colnames]
  M_disease_var = coe_change[se_colnames] ** 2
  assert(M_disease.shape == (num_cpg_locations, num_cell_types))
  assert(M_disease_var.shape == (num_cpg_locations, num_cell_types))

  return M_control, M_control_var, M_disease, M_disease_var


def bulk_control_and_disease_mean(coe_control, coe_change):
  """
  Use bulk DMC linear regression results to compute the mean M-values we would expect for a control
  and disease person.

  coe_control (pd.DataFrame) : Has a row for each CpG and a columns for the linear regression
                               results for the Intercept and each cell type.
  coe_change (pd.DataFrame) : Has a row for each CpG and a columns for the linear regression
                              results for each cell type.
  """
  B_control = coe_control["(Intercept).Estimate"]
  B_control_var = coe_control["(Intercept).SE"] ** 2

  B_disease = B_control + coe_change["Estimate"]
  B_disease_var = B_control_var + coe_change["SE"] ** 2

  return B_control, B_control_var, B_disease, B_disease_var


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
