import os
import sys
import glob
import argparse
import pandas as pd
import variables as vrb
import multiprocessing as mp
from AdaptivePELE.atomset import atomset
import AdaptivePELE.utilities.utilities as adapt_tools


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-sp", "--sim_path", required=True,
                                help="""Absolute path where all simulations that will be analyzed are stored.""")
    required_named.add_argument("-o", "--out_path", required=True,
                                help="""Absolute path where the final csvs will be stored.""")

    parser.add_argument("-r", "--report_pattern", default="report_",
                        help="""PELE's reports prefix pattern.""")
    parser.add_argument("-t", "--traject_pattern", default="trajectory_",
                        help="""Prefix pattern of trajectory files.""")
    parser.add_argument("-ap", "--adaptive_results_pattern", default="obc_adaptive_output*",
                        help="""Pattern of folders to simulation results.""")
    parser.add_argument("-p", "--pdb", default=False,
                        help="If True, it converts all .xtc files to pdb format")

    args = parser.parse_args()

    return args.sim_path, args.out_path, args.report_pattern, args.traject_pattern, args.adaptive_results_pattern, \
           args.pdb


def trajectory_and_snapshot_to_pdb(trajectory_path, snapshot, output_path):
    """
    Given an absolute path to a trajectory of Adaptive and a snapshot (MODEL) in xtc format, the function transform it
    into a PDB format.
    :param trajectory_path: Absolute path to a trajectory from Adaptive, in xtc format.
    :type trajectory_path:str
    :param snapshot: model of a trajectory that you want to transform.
    :type snapshot: int
    :param output_path: output path of the new pdb file.
    :type output_path: str
    :return: Creates a PDB file.
    """
    topology_path_splited = trajectory_path.split("/")[0:-2]
    topology_path = os.path.join("/".join(topology_path_splited), "topology.pdb")
    topology_contents = adapt_tools.getTopologyFile(topology_path)
    trajectory = adapt_tools.getSnapshots(trajectory_path, topology=topology_path)
    try:
        single_model = trajectory[snapshot]
        PDB = atomset.PDB()
        PDB.initialise(single_model, topology=topology_contents)
    except IndexError:
        exit("You are selecting the model {} for a trajectory that has {} models, please, reselect the model index "
             "(starting from 0).".format(snapshot, len(trajectory)))
    with open(output_path, "w") as fw:
        fw.write("MODEL     %4d\n" % (snapshot + 1))
        fw.write(PDB.pdb)
        fw.write("ENDMDL\n")
        fw.write("END\n")


def get_pdb_from_xtc(row, pdbs_output_path, column_file="file_from"):
    """
    Given a row of a pandas.DataFrame this function search the column with the structure in xtc format that corresponds
    to this row's information and extract a pdb file in an output folder.
    :param row: row of a dataframe in pandas.
    :type row: pandas.DataFrame row
    :param pdbs_output_path: path to the folder where the user wants to store the pdb file.
    :type pdbs_output_path: str
    :param column_file: name of the column that contains the path to the xtc file.
    :type column_file: str
    :return: creates the pdb into the pdbs_output_path and prints the path to the pdb file.
    """
    foldername = row[column_file]
    filepath = glob.glob(foldername)[0]
    epoch = filepath.split("/")[-2]
    snapshot = row["numberOfAcceptedPeleSteps"]
    new_file_name = os.path.basename(foldername.split("/")[-1])
    new_file_name = new_file_name.split(".")[0]
    trajectory_and_snapshot_to_pdb(filepath, snapshot, os.path.join(pdbs_output_path, "{}_epoch_{}_snap_{}.pdb".format(
        new_file_name, epoch, snapshot)
                                                                    ))
    print(os.path.join(pdbs_output_path, "{}_epoch_{}_snap_{}.pdb".format(new_file_name, epoch, snapshot)))


def get_pdbs_from_df_in_xtc(df, pdbs_output_path, processors=4, column_file="file_from"):
    """
    It uses the function "get_pdb_from_xtc" for a whole dataframe using multiprocessing.
    :param df: Dataframe object (Pandas)
    :type df: pandas.DataFrame
    :param pdbs_output_path: Output path for PDB files.
    :type pdbs_output_path: str
    :param processors: Number of processes to do with multiprocessing.
    :type processors: int
    :param column_file: Column name of the dataframe that contains the path to the trajectory file.
    :type column_file: str
    :return:
    """
    pool = mp.Pool(processes=processors)
    multiprocessing_list = []
    for index, row in df.iterrows():
        multiprocessing_list.append(pool.apply_async(get_pdb_from_xtc,
                                                     (row, pdbs_output_path, column_file)))
    for process in multiprocessing_list:
        process.get()


def concat_reports_in_csv(adaptive_results_path, output_file_path, report_prefix="report_",
                          trajectory_prefix="trajectory_"):
    dataframe_lists = []
    for adaptive_epoch in range(0, 2000):
        folder = os.path.join(adaptive_results_path, str(adaptive_epoch))
        if os.path.exists(folder):
            report_list = glob.glob("{}/*{}*".format(folder, report_prefix))
            report_list = sorted(report_list, key=lambda x: int(x.split("_")[-1]))
            for n, report in enumerate(report_list):
                pandas_df = pd.read_csv(report, sep="    ", engine="python", index_col=False, header=0)
                pandas_df["epoch"] = adaptive_epoch
                pandas_df["file_from"] = glob.glob("{}/{}/*{}{}.*".format(adaptive_results_path, adaptive_epoch,
                                                                          trajectory_prefix, n + 1))[0]
                dataframe_lists.append(pandas_df)
        else:
            break
    dataframe = pd.concat(dataframe_lists, ignore_index=True)
    dataframe.to_csv(output_file_path, sep=';', index=False)


def main(simulations_path, output_file_path, report_prefix="report_", trajectory_prefix="trajectory_",
         path_to_adaptive_results=vrb.PATH_PATTER_TO_ADAPTIVE_RESULTS, extract_pdbs=True):
    path_list = glob.glob("{}/[0-9]_*".format(simulations_path))
    if not os.path.exists(output_file_path):
        os.mkdir(output_file_path)
    for path in path_list:
        path_to_result = glob.glob("{}/{}".format(path, path_to_adaptive_results))[0]
        print("Analysing {}".format(path_to_result))
        output_name = os.path.join(output_file_path, "{}_adaptive_results_summary.csv".format(os.path.basename(path)))
        concat_reports_in_csv(path_to_result, output_name, report_prefix, trajectory_prefix)
        print("Results saved in {}".format(output_name))
        if extract_pdbs:
            df = pd.read_csv(output_name, sep=";")
            get_pdbs_from_df_in_xtc(df, output_file_path, column_file="file_from",
                                    processors=4)
    return print("Finished successfully!")


if __name__ == '__main__':
    sim_path, out_path, report_pattern, traject_pattern, adaptive_results_pattern, pdb = parse_arguments()
    main(sim_path, out_path, report_pattern, traject_pattern, adaptive_results_pattern, extract_pdbs=pdb)
