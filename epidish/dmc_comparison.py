import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict
from naive_bayes_utils import *


def plot_pvalue_histogram(analysis_folder, cell_types=None):
  coe_change = pd.read_csv(os.path.join(analysis_folder, "coe_change.csv"), index_col=0)

  if cell_types is not None:
    p_colnames = ["{}.p".format(cellname) for cellname in cell_types]
    adjp_colnames = ["{}.adjP".format(cellname) for cellname in cell_types]
  else:
    p_colnames = ["p"]
    adjp_colnames = ["adjP"]

  for i, adjp_colname in enumerate(adjp_colnames):
    cellname = adjp_colname.split("adjP")[0][:-1]
    adjP = coe_change[adjp_colname]
    p = coe_change[p_colnames[i]]

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)

    ax1.hist(adjP, bins=100)
    ax1.set_title("Adjusted p-values for {}".format(cellname))
    ax1.set_xlabel("p-value")
    ax1.set_ylabel("CpG count")

    ax2.hist(p, bins=100)
    ax2.set_title("Raw p-values for {}".format(cellname))
    ax2.set_xlabel("p-value")
    ax2.set_ylabel("CpG count")

    fig.set_size_inches(15, 8)

    save_path = os.path.join(analysis_folder, "pvalue_hist_{}.png".format(cellname))
    fig.savefig(save_path, dpi=100)
    # plt.show()


def compare_dmcs(dataset_dmc_sets):
  """
  Report the DMCs that are (1) unique to each dataset and (2) shared between each pair of datasets.
  """
  print("---------------------------------------------------------------------")
  names = sorted(list(dataset_dmc_sets.keys()))
  summary = np.zeros((len(names), len(names)), dtype=np.int32)
  unique_dmcs = dataset_dmc_sets.copy()
  print(names)

  print("\n============ SHARED DMCs ============")
  for i in range(len(names)-1):
    for j in range(i+1, len(names)):
      name = names[i]
      other = names[j]
      # Take the intersect of CpG sets.
      shared_dmcs = dataset_dmc_sets[name] & dataset_dmc_sets[other]

      summary[i, j] = len(shared_dmcs)
      summary[j, i] = len(shared_dmcs)

      # Remove any shared DMCs from each dataset's unique set.
      unique_dmcs[name] = unique_dmcs[name].difference(dataset_dmc_sets[other])
      unique_dmcs[other] = unique_dmcs[other].difference(dataset_dmc_sets[name])

      print("\n==> Datasets {} and {} share {} DMCs:".format(name, other, len(shared_dmcs)))
      for dmc in shared_dmcs:
        print(" >>", dmc)

  # print("\n============ UNIQUE DMCs =============")
  for i, name in enumerate(names):
    # print("\n==> Dataset {} has {} unique DMCs:".format(name, len(unique_dmcs[name])))
    # for dmc in unique_dmcs[name]:
    #   print(" >>", dmc)
    summary[i, i] = int(len(unique_dmcs[name]))
  
  print("\n============== NUMBER OF CpGs SHARED ==============")
  summary = pd.DataFrame(summary)
  summary.columns = names
  summary.index = names
  print(summary)

  summary.to_csv("./shared_cpg_matrix.csv")
  return summary

  # for i in range(len(names)):
  #   for j in range(len(names)):
  #     summary.iloc[i,j] /= min(len(dataset_dmc_sets[names[i]]), len(dataset_dmc_sets[names[i]]))

  # print("\n============= PERCENTAGE CpGs SHARED ===============")
  # print(summary)


def save_signif_cpg_files(folders, cell_types_for_each, p_value_thresh):
  for i, analysis_folder in enumerate(folders):
    cell_types = cell_types_for_each[i]
    coe_control, coe_change, cell_fracs, phenotypes, beta = \
        load_epidish_results(analysis_folder, mvalues=True, has_cellfrac=False)

    if cell_types is not None:
      signif_cpg = report_significant_cpgs(cell_types, coe_change, p_value_thresh=p_value_thresh)
    else:
      signif_cpg = report_significant_cpgs_bulk(coe_change, p_value_thresh=p_value_thresh)

    print("==> All significant CpG locations:")
    print(signif_cpg)

    with open(os.path.join(analysis_folder, "signif_cpg.txt"), "w") as f:
      for cpg in signif_cpg:
        f.write(cpg + "\n")


def load_signif_cpg_files(folders, dataset_names):
  dataset_dmc_sets = defaultdict(lambda: set())
  for i, analysis_folder in enumerate(folders):
    with open(os.path.join(analysis_folder, "signif_cpgs.txt"), "r") as f:
      for l in f:
        l = l.replace("\n", "")
        dataset_dmc_sets[dataset_names[i]].add(l)

  return dataset_dmc_sets


def load_signif_cpg_files_full_path(paths, dataset_names):
  dataset_dmc_sets = defaultdict(lambda: set())
  for i, path in enumerate(paths):
    with open(path, "r") as f:
      for l in f:
        l = l.replace("\n", "")
        dataset_dmc_sets[dataset_names[i]].add(l)

  return dataset_dmc_sets


def compare_feature_dmcs_main():
  dmc_feature_folders = [
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_all/results_0.06/",
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_bulk/results_0.7/",
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_only_pbmc/results_0.01/",
    "../analysis/martino2018/Mvalues_control_vs_allergic/results_0.25/",
    "../analysis/martino2018/Mvalues_control_vs_allergic_bulk/results_0.01/",
    "../analysis/martino2018/shap/",
    "../analysis/martino2015/shap/"
  ]

  names = [
    "2015_ALL",
    "2015_BULK",
    "2015_PBMC",
    "2018_CD4_CD8",
    "2018_CD4",
    "2018_XGB_SHAP",
    "2015_XGB_SHAP"
  ]

  dataset_dmc_sets = load_signif_cpg_files(dmc_feature_folders, names)
  compare_dmcs(dataset_dmc_sets)


def compare_signif_dmcs_main():
  dmc_cutoff_01 = [
    "../analysis/signif_cpgs/2015_all.txt",
    "../analysis/signif_cpgs/2015_bulk.txt",
    "../analysis/signif_cpgs/2015_pbmc_0.1.txt",
    "../analysis/signif_cpgs/2018_cd4_cd8_0.1.txt",
    "../analysis/signif_cpgs/2018_bulk_0.1.txt",
    "../analysis/signif_cpgs/shap_2018.txt",
    "../analysis/signif_cpgs/shap_2015.txt"
  ]

  names = [
    "NB_2015_ALL",
    "NB_2015_BULK",
    "NB_2015_PBMC",
    "NB_2018_CD4_CD8",
    "NB_2018_CD4",
    "XGB_2018_BULK",
    "XGB_2015_BULK"
  ]

  dataset_dmc_sets = load_signif_cpg_files_full_path(dmc_cutoff_01, names)
  names = sorted(dataset_dmc_sets.keys())
  dmcs = compare_dmcs(dataset_dmc_sets)

  fig, ax = plt.subplots()
  im = ax.imshow(dmcs.to_numpy() > 0, cmap="coolwarm")

  ax.set_xticks(np.arange(len(names)))
  ax.set_yticks(np.arange(len(names)))
  ax.set_xticklabels(names)
  ax.set_yticklabels(names)

  # Rotate the tick labels and set their alignment.
  plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
          rotation_mode="anchor")

  # Loop over data dimensions and create text annotations.
  for i in range(len(names)):
    for j in range(len(names)):
      text = ax.text(j, i, dmcs.to_numpy()[i, j], ha="center", va="center", color="w")

  ax.set_title("Comparison of DMCs from Naive Bayes and Random Forest")
  fig.tight_layout()
  plt.show()


def pvalue_histogram_main():
  folders = [
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_eosino/",
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_neutro/",
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_all/",
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_only_pbmc/",
    "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_bulk/",
    "../analysis/martino2018/Mvalues_control_vs_allergic/",
    "../analysis/martino2018/Mvalues_control_vs_allergic_bulk/"
  ]

  cell_types_for_each = [
    ["B", "NK", "CD4T", "CD8T", "Mono", "Eosino"],
    ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro"],
    ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"],
    ["B", "NK", "CD4T", "CD8T", "Mono"],
    None,
    ["CD4T", "CD8T"],
    None
  ]

  for i in range(len(folders)):
    folder = folders[i]
    print(">> Plotting histograms for folder {}".format(folder))
    cell_types = cell_types_for_each[i]
    plot_pvalue_histogram(folder, cell_types)


if __name__ == "__main__":
  # pvalue_histogram_main()

  # Compare shared DMCs across all datasets, using a p-value threshold of 0.1.
  compare_signif_dmcs_main()

  # Compare shared DMC features across all CLASSIFIERS.
  # compare_feature_dmcs_main()
