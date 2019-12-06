from classifier_utils import *


def write_dmc_csv(folder):
  dmct = pd.read_csv(os.path.join(folder, "dmct.csv"), index_col=0)
  dmct = dmct[dmct["DMC"] != 0]
  dmct = dmct.drop(columns="DMC")
  dmct.to_csv(os.path.join(folder, "dmct_small.csv"))


if __name__ == "__main__":
  # folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_only_pbmc/"
  # folder = "../analysis/martino2015/Mvalues_nonallergic_vs_allergic_all/"
  folder = "../analysis/martino2018/Mvalues_control_vs_allergic/"

  write_dmc_csv(folder)
