import os
import sys
import glob
import argparse
import shutil
import pandas as pd
import variables as vrb
from sims_and_glide_to_csv import report_to_dataframe


def extract_top_structures_from_report(report_path, column, path_to_sims, output_folder, top_n=10, ascending=True,
         path_to_pdbs=vrb.PATH_PATTER_TO_CLUSTERS_FROM_SIMS):
    report = report_to_dataframe(report_path, separator=";")
    report_sorted = report.sort_values(by=[column], ascending=ascending)
    top10 = report_sorted.iloc[0:int(top_n)]
    compound_name = top10.iloc[0]["compound_name"]
    print(compound_name)
    out_dit_for_compound = os.path.join(output_folder, compound_name)
    top10.to_csv(os.path.join(out_dit_for_compound, "report_top_{}_{}.csv".format(top_n, compound_name)), sep=";")
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    for n, (index, row) in enumerate(top10.iterrows()):
        ID = row["ID"]
        compound_name = row["compound_name"]
        filepath = glob.glob("{}.pdb".format(os.path.join(path_to_sims, compound_name, compound_name,
                                                          path_to_pdbs, ID)))
        out_dit_for_compound = os.path.join(output_folder, compound_name)
        if not os.path.exists(out_dit_for_compound):
            os.mkdir(out_dit_for_compound)
        shutil.copy(filepath[0], os.path.join(output_folder, out_dit_for_compound, "{}_{}_{}_{}.pdb".format(
            compound_name, ID, column, n)))


def main(reports_folder, column, path_to_sims, output_folder, top_n=10, ascending=True,
         path_to_pdbs=vrb.PATH_PATTER_TO_CLUSTERS_FROM_SIMS):
    csvs_in_path = glob.glob("{}/*.csv".format(reports_folder))
    for csv in csvs_in_path:
        extract_top_structures_from_report(csv, column, path_to_sims, output_folder, top_n, ascending, path_to_pdbs)
    return print("Finished!")


main("/home/carlespl/project/Almirall/results_summary", "Binding_Energy",
     "/home/carlespl/project/Almirall/simulation_results", "/home/carlespl/project/Almirall/top_structures")

