import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal


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


def cell_methylation_matrices(coe_control, coe_change, cell_types):
  """
  Build a two matrices with a row for each CpG location, and a column for each cell type. Each entry
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
