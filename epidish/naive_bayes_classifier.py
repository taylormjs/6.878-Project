import os
import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal
from matplotlib import pyplot as plt

from naive_bayes_utils import *

# NOTE(milo): Make sure to change this once we figure out where sensitized belongs!
MARTINO2015_LABEL_MAP = {"nonallergic": 0, "allergic": 1, "sensitized": 2}
MARTINO2018_LABEL_MAP = {"control": 0, "allergic": 1, "resolved": 2}


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

  return likely_ratios, signif_cpg, coe_change.loc[signif_cpg]


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

  return likely_ratios, signif_cpg, coe_change.loc[signif_cpg]


def cs_2015_classifier(p_value_thresh, likelihood_ratio, write_results=True):
  """
  Cell-specific Naive-Bayes classifier for Martino 2015.
  """
  # analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_all/" 
  # analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_neutro/"
  # analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_eosino/"
  analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_only_pbmc/"
  # NOTE(milo): Only use PBMC cell types here.
  cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono"]
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro"]
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Eosino"]
  # cell_types_m2015 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]

  likely_ratios, signif_cpgs, coe_change_signif = \
      run_cell_specific_classifier(analysis_folder, cell_types_m2015, use_mvalues=True, p_value_thresh=p_value_thresh)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "../phenotypes.csv"), index_col=0)

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

  if write_results:
    write_classifier_results(analysis_folder, p_value_thresh, likely_ratios, labels, precision,
                             recall, signif_cpgs, coe_change_signif)

  return precision, recall


def bulk_2015_classifier(p_value_thresh, likelihood_ratio, write_results=True):
  """
  Bulk Naive-Bayes classifier for Martino 2015.
  """
  test_patients = list(pd.read_csv("../analysis/martino2015/test_set.txt", header=None)[0])
  analysis_folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_bulk/"

  likely_ratios, signif_cpgs, coe_change_signif \
      = run_bulk_classifier(analysis_folder, use_mvalues=True, p_value_thresh=p_value_thresh)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "../phenotypes.csv"), index_col=0)

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

  if write_results:
    write_classifier_results(analysis_folder, p_value_thresh, likely_ratios, labels, precision,
                             recall, signif_cpgs, coe_change_signif)

  # print(labels)
  return precision, recall


def bulk_2018_classifier(p_value_thresh, likelihood_ratio, write_results=True):
  """
  Bulk Naive-Bayes classifier for Martino 2018.

  NOTE(milo): Because Martino2018 data is supposedly only CD4T cells, "bulk" means CD4T only in this
  case.
  """
  test_patients = list(pd.read_csv("../analysis/martino2018/test_set.txt", header=None)[0])
  analysis_folder = "../analysis/martino2018/Mvalues_control_vs_allergic_bulk/"

  likely_ratios, signif_cpgs, coe_change_signif = \
      run_bulk_classifier(analysis_folder, use_mvalues=True, p_value_thresh=p_value_thresh)
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "../phenotypes.csv"), index_col=0)

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

  if write_results:
    write_classifier_results(analysis_folder, p_value_thresh, likely_ratios, labels, precision,
                             recall, signif_cpgs, coe_change_signif)

  # print(labels)
  return precision, recall


def cs_2018_classifier(p_value_thresh, likelihood_ratio, write_results=True):
  """
  Cell-specific Naive-Bayes classifier for Martino 2018.
  """
  test_patients = list(pd.read_csv("../analysis/martino2018/test_set.txt", header=None)[0])
  analysis_folder = "../analysis/martino2018/Mvalues_control_vs_allergic/"
  cell_types_m2018 = ["CD4T", "CD8T"]

  likely_ratios, signif_cpgs, coe_change_signif = run_cell_specific_classifier(
      analysis_folder, cell_types_m2018, use_mvalues=True, p_value_thresh=p_value_thresh)
  
  phenotypes = pd.read_csv(os.path.join(analysis_folder, "../phenotypes.csv"), index_col=0)

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

  if write_results:
    write_classifier_results(analysis_folder, p_value_thresh, likely_ratios, labels, precision,
                             recall, signif_cpgs, coe_change_signif)

  return precision, recall


def write_classifier_results(analysis_folder, p_value_thresh, likely_ratios, labels, precision,
                             recall, signif_cpgs, coe_change_signif):
  """
  Write results from running a classifier to disk.
  """
  results_folder = os.path.join(analysis_folder, "results_{}".format(p_value_thresh))

  if not os.path.exists(os.path.abspath(results_folder)):
    print("NOTE: Making path {}".format(os.path.abspath(results_folder)))
    os.mkdir(os.path.abspath(results_folder))

  with open(os.path.join(results_folder, "likelihood_ratios.txt"), "w") as f:
    f.write("patient likelihood_ratio\n")
    for patient in likely_ratios:
      f.write("{} {}\n".format(patient, likely_ratios[patient]))

  with open(os.path.join(results_folder, "predicted_pheno.txt"), "w") as f:
    f.write("patient predicted_pheno true_pheno\n")
    for patient in labels:
      f.write("{} {} {}\n".format(patient, labels[patient][0], labels[patient][1]))

  with open(os.path.join(results_folder, "signif_cpgs.txt"), "w") as f:
    for cpg in signif_cpgs:
      f.write(str(cpg) + "\n")

  coe_change_signif.to_csv(os.path.join(results_folder, "coe_change_signif.txt"))

  with open(os.path.join(results_folder, "meta.txt"), "w") as f:
    f.write("key value\n")
    f.write("p_value_thresh {}\n".format(p_value_thresh))
    f.write("precision {}\n".format(precision))
    f.write("recall {}\n".format(recall))
    f.write("signif_cpgs {}\n".format(len(signif_cpgs)))


def precision_recall_vs_pvalue_thresh(function_to_run):
  """
  Look at the classifier performance as we vary the p_value_threshold for choosing CpG features.
      
  Outputs a text file where each line has:
  p_value_thresh    precision    recall
  """
  likelihood_ratio = 1

  # for p_value_thresh in np.arange(0.05, 0.55, 0.05):
  for p_value_thresh in np.arange(0.005, 0.1, 0.005):
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

    for lr in np.arange(0.1, 100, 0.1):
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


def report_significant_cpgs_main():
  analysis_folders = [
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_bulk/",
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_all/" ,
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_only_pbmc/",
    "../analysis/martino2018/Mvalues_control_vs_allergic_bulk/",
    "../analysis/martino2018/Mvalues_control_vs_allergic/"
  ]

  cell_types = [
    None,
    ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"],
    ["B", "NK", "CD4T", "CD8T", "Mono"],
    None,
    ["CD4T", "CD8T"]
  ]

  p_value_thresh = [
    0.70,
    0.06,
    0.01,
    0.01,
    0.25
  ]

  for i in range(len(analysis_folders)):
    folder = analysis_folders[i]
    ct = cell_types[i]
    print("=============== DATASET: {} ================".format(folder))
    print(">> p_value_thresh={}".format(p_value_thresh[i]))

    coe_change = pd.read_csv(os.path.join(folder, "coe_change.csv"), index_col=0)

    if ct is None:
      report_significant_cpgs_bulk(coe_change, p_value_thresh[i])
    else:
      report_significant_cpgs(ct, coe_change, p_value_thresh[i])


if __name__ == "__main__":
  #==================== RUN CLASSIFIERS WITH BEST PARAMS =====================
  # bulk_2015_classifier(0.60, 1) # Do for BULK.
  # cs_2015_classifier(0.01, 1) # Do for PBMC.
  # cs_2015_classifier(0.06, 1) # Do for ALL.

  # bulk_2018_classifier(0.01, 1)
  # cs_2018_classifier(0.25, 1)

  #=================== PRECISION-RECALL PLOTTING ============================
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
  # 2018 CD4: 0.01
  # 2018 CD4 and CD8: 0.25

  # lr_files = {
  #   "2015_all": ("../analysis/likelihood_ratios/2015_all.txt", 2015),
  #   "2015_pbmc": ("../analysis/likelihood_ratios/2015_pbmc.txt", 2015),
  #   "2018_cd4": ("../analysis/likelihood_ratios/2018_cd4.txt", 2018),
  #   "2018_cd4_cd8": ("../analysis/likelihood_ratios/2018_cd4_cd8.txt", 2018)
  # }
  # precision_recall_vs_cutoff(lr_files)

  report_significant_cpgs_main()
