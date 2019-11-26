import numpy as np
from collections import defaultdict
from classifier_utils import *


def compare_dmcs(dataset_dmc_sets):
  """
  Report the DMCs that are (1) unique to each dataset and (2) shared between each pair of datasets.
  """
  print("---------------------------------------------------------------------")
  names = list(dataset_dmc_sets.keys())
  summary = np.zeros((len(names), len(names)))
  unique_dmcs = dataset_dmc_sets.copy()

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

  print("\n============ UNIQUE DMCs =============")
  for i, name in enumerate(names):
    print("\n==> Dataset {} has {} unique DMCs:".format(name, len(unique_dmcs[name])))
    for dmc in unique_dmcs[name]:
      print(" >>", dmc)
    summary[i, i] = len(unique_dmcs[name])
  
  print("\n============== SUMMARY ==============")
  print(summary)

  for i in range(len(names)):
    for j in range(len(names)):
      summary[i,j] /= min(len(dataset_dmc_sets[names[i]]), len(dataset_dmc_sets[names[i]]))

  print(summary)


def save_signif_cpg_files(folders, cell_types_for_each):
  for i, analysis_folder in enumerate(folders):
    cell_types = cell_types_for_each[i]
    coe_control, coe_change, cell_fracs, phenotypes, beta = load_epidish_results(analysis_folder)
    signif_cpg = report_significant_cpgs(cell_types, coe_change, p_value_thresh=0.1)
    print("==> All significant CpG locations:")
    print(signif_cpg)

    with open(os.path.join(analysis_folder, "signif_cpg.txt"), "w") as f:
      for cpg in signif_cpg:
        f.write(cpg + "\n")


def load_signif_cpg_files(folders, dataset_names):
  dataset_dmc_sets = defaultdict(lambda: set())
  for i, analysis_folder in enumerate(folders):
    with open(os.path.join(analysis_folder, "signif_cpg.txt"), "r") as f:
      for l in f:
        l = l.replace("\n", "")
        dataset_dmc_sets[dataset_names[i]].add(l)

  return dataset_dmc_sets


if __name__ == "__main__":
  folders = [
    "../analysis/martino2015/nonallergic_vs_allergic_only_pbmc/",
    "../analysis/martino2015/nonallergic_vs_allergic_with_eosino/",
    "../analysis/martino2015/nonallergic_vs_allergic_with_neutro/",
    "../analysis/martino2018/control_vs_allergic/"
  ]
  names = [
    "2015_only_pbmc",
    "2015_with_eosino",
    "2015_with_neutro",
    "2018"
  ]
  cell_types_for_each = [
    ["B", "NK", "CD4T", "CD8T", "Mono"],
    ["B", "NK", "CD4T", "CD8T", "Mono", "Eosino"],
    ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro"]
  ]

  # save_signif_cpg_files(folders, cell_types_for_each)
  dataset_dmc_sets = load_signif_cpg_files(folders, names)
  compare_dmcs(dataset_dmc_sets)
