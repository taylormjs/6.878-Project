import numpy as np
import pandas as pd


def load_epidish_results_martino2015():
  coe_control = pd.read_csv("./analysis/martino2015_coe_control.csv")
  coe_change = pd.read_csv("./analysis/martino2015_coe_change.csv")
  return coe_control, coe_change


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
  M_control_stdev = coe_control[se_colnames]
  assert(M_control.shape == (num_cpg_locations, num_cell_types))
  assert(M_control_stdev.shape == (num_cpg_locations, num_cell_types))

  # Make the disease matrix (control + change).
  M_disease = M_control + coe_change[estimate_colnames]
  M_disease_stdev = M_control_stdev + coe_change[se_colnames]
  assert(M_disease.shape == (num_cpg_locations, num_cell_types))
  assert(M_disease_stdev.shape == (num_cpg_locations, num_cell_types))

  return M_control, M_control_stdev, M_disease, M_disease_stdev


def report_significant_cpgs(cell_types, coe_change, p_value_thresh=0.05):
  """
  Print out the number of significant CpG locations that were found for each cell type.
  """
  for cell_type in cell_types:
    adjP_colname = "{}.adjP".format(cell_type)
    p_values = coe_change[coe_change[adjP_colname] <= p_value_thresh]
    num_signif = p_values.shape[0]
    print("Found {} significant CpGs for {}".format(cell_type, num_signif))


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


def run_martino2015():
  coe_control, coe_change = load_epidish_results_martino2015()

  # NOTE(milo): Had to throw away Neutro in R because it was causing errors.
  cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Eosino"]

  # print(coe_control.columns)
  # print(coe_change.columns)
  coe_control = rename_control_cols(coe_control, cell_types_m2015)
  # print(coe_control.columns)

  # report_significant_cpgs(cell_types_m2015, coe_change)

  Mc, Mc_stdev, Md, Md_stdev = cell_methylation_matrices(coe_control, coe_change, cell_types_m2015)


if __name__ == "__main__":
  run_martino2015()
