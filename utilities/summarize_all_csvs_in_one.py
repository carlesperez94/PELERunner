import os
import sys
import glob
import argparse
import shutil
import pandas as pd
from AdaptivePELE.atomset import atomset
import AdaptivePELE.utilities.utilities as adapt_tools

def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""This script do two different things: 1) selects the best (by certain
    criteria and column) row for all csv files found in the input folder and 2) extract the structure correspondent to 
    this row.
    """)
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("path_results",
                                help="""Absolute path to the folder with all csvs that will be analyzed.""")
    parser.add_argument("-c", "--column",
                        help="""Column of the csv file that will be used to select the best row.""")
    parser.add_argument("-o", "--output_filepath",
                        help="""Output absolut path to the final csv.""")
    parser.add_argument("-cr", "--criteria", choices=['max', 'min'], default="min",
                        help="""Criteria that will be applied to the column to select the best row: choice
                        between 'min' and 'max'.""")
    parser.add_argument("-sp", "--simulation_path", default=None,
                        help="""If the user wants to get the structures, this variable has to be filled with the 
                        absolute path where the simulations are stored.""")
    parser.add_argument("-p2p", "--path_to_pdbs", default="*/*/",
                        help="""Pattern of folders to get the structures.""")
    parser.add_argument("-id", "--id_name", default="ID",
                        help="""Column name of the csv where is the name of the pdb structure for each row.""")
    parser.add_argument("-f", "--folder_name", default="filename",
                        help="""Column name of the csv where is the name of the file/folder where structures
                        are stored.""")
    parser.add_argument("-po", "--pdbs_output_path", default=None,
                        help="""Folder's absolute path where the structures will be stored.""")
    parser.add_argument("-x", "--xtc", default=False,
                        help="""Flag that has to be set when the structures are in 'xtc' format instead of pdb.""")

    args = parser.parse_args()

    return args.path_results, args.column, args.criteria, args.output_filepath, args.simulation_path, \
           args.path_to_pdbs, args.id_name, args.folder_name, args.pdbs_output_path, args.xtc


def read_all_dataframes_of_a_folder(path_to_results):
    """
    Given a path read all dataframes inside the folder.
    :param path_to_results: path to the folder that contains different csv files.
    :type path_to_results: str
    :return: list of pandas' dataframes
    """
    dataframes = []
    if path_to_results.endswith("/"):
        path_to_results = path_to_results.rstrip("/")
    csv_files = glob.glob("{}/*.csv".format(path_to_results))
    for csv in csv_files:
        dataframe = pd.read_csv(csv, sep=";", engine="python")
        add_filename_to_csv_as_column(dataframe, csv)
        dataframes.append(dataframe)
    return dataframes


def select_row_of_dataframe_by_criteria(dataframe, column, criteria):
    """
    Select a row of a dataframe by the minimum or maximum value of a column.
    :param dataframe: Dataframe to analyze.
    :type dataframe: pandas.DataFrame.
    :param column: Name of the column of the dataframe that we want to use to select the minimum or maximum value.
    :type column: str
    :param criteria: Criteria that will be use to select the row: 'min' or 'max' (minimum or maximum).
    :type criteria: str
    :return: row of a pandas.DataFrame
    """
    if criteria == "min":
        minimum = dataframe.loc[dataframe[column] == dataframe[column].min()]
        if len(minimum) > 1:
            return minimum.iloc[[0]]
        else:
            return minimum
    if criteria == "max":
        maximum = dataframe.loc[dataframe[column] == dataframe[column].max()]
        if len(maximum) > 1:
            return maximum.iloc[[0]]
        else:
            return maximum


def add_filename_to_csv_as_column(dataframe, filepath, column_name="origin_name"):
    """
    Adds a column with the name of the file to a pandas.DataFrame
    :param dataframe: dataframe where we want to add the column.
    :type dataframe: pandas.DataFrame
    :param filepath: path to the file that we want to add in the column.
    :type filepath: str
    :param column_name: name of the new column added.
    :type column_name: str
    :return: pandas.DataFrame with the new column incorporated.
    """
    filename = os.path.basename(filepath).split(".")[0]
    dataframe[column_name] = filename
    return dataframe


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


def compute_mean_quantile(dataframe, column, quantile_value):
    """
    It computes the mean of the sample selected into the specific quantile of a column for a pandas.DataFrame.
    :param dataframe: dataframe that contains the data.
    :type dataframe: pandas.Dataframe
    :param column: Name of the column of the variable that will be used to compute the mean of the quantile.
    :type column: str
    :param quantile_value: value of the quantile to subset the data of the column selected.
    :type quantile_value: float
    :return: value of the mean
    """
    dataframe = dataframe[dataframe[column] < dataframe[column].quantile(quantile_value)]
    dataframe = dataframe[column]
    mean_subset = dataframe.mean()
    return mean_subset


def main(path_to_results, column="Binding Energy", criteria="min", output_filepath="csv_summary.csv",
         simulation_results_path=None, path_to_pdbs="*/*/", id_name="ID", folder_name="filename",
         pdbs_output_path=None, xtc=False):
    """
    1) selects the best (by certain criteria and column) row for all csv files found in the input folder and 2) extract
    the structure correspondent to this row.
    :param path_to_results: Path of the folder where all results are stored.
    :type path_to_results: str
    :param column: Name of the column of the dataframe that we want to use to select the minimum or maximum value, and
    compute the mean of the quartile.
    :type column:str
    :param criteria: Criteria that will be use to select the row: 'min' or 'max' (minimum or maximum).
    :type criteria: str
    :param output_filepath: Path to the output csv file.
    :type output_filepath: str
    :param simulation_results_path: Path of the folder where the simulation results are stored.
    :type simulation_results_path: str
    :param path_to_pdbs: Path to the folder where the trajectories are stored.
    :type path_to_pdbs: str
    :param id_name: Name of the column of the CSV file that will contain the compound name (extracted from the folder).
    :type id_name: str
    :param folder_name: Name of the column of the CSV file that will contain the path to the structure that has been
    selected.
    :type folder_name: str
    :param pdbs_output_path: Path to the folder that will contain the pdb structures.
    :type pdbs_output_path: str
    :param xtc: If it is set, indicates that the structures of the simulation are in 'xtc' format.
    :type xtc: bool
    :return:
    """
    dataframe_list = read_all_dataframes_of_a_folder(path_to_results)
    list_of_selected = []
    for df in dataframe_list:
        mean_quartile = compute_mean_quantile(df, column, quantile_value=0.10)
        df["mean_quartile"] = mean_quartile
        selected_row = select_row_of_dataframe_by_criteria(df, column, criteria)
        list_of_selected.append(selected_row)
    single_df = pd.concat(list_of_selected, ignore_index=True)
    single_df.to_csv(output_filepath, sep=";", index=False)
    if pdbs_output_path:
        if not os.path.exists(pdbs_output_path):
            os.mkdir(pdbs_output_path)
    if simulation_results_path:
        for index, row in single_df.iterrows():
            foldername = row[folder_name]
            if xtc:
                filepath = glob.glob(foldername)[0]
                snapshot = row["numberOfAcceptedPeleSteps"]
                epoch = filepath.split("/")[-2]
                new_file_name = os.path.basename(foldername.split("/")[-1])
                new_file_name = new_file_name.split(".")[0]
                trajectory_and_snapshot_to_pdb(filepath, snapshot, os.path.join(os.path.join(pdbs_output_path,
                                                                                             "{}_epoch_{}_snap_{}.pdb".format(
                                                                                                 new_file_name, epoch,
                                                                                                 snapshot)
                                                                                             )))
                print(os.path.join(pdbs_output_path, "{}_epoch_{}_snap_{}.pdb".format(new_file_name, epoch, snapshot)))
            else:
                ID = row[id_name]
                if "recomputed" in foldername:
                    print("CSV recomputed found!")
                    foldername = foldername.split("_")[0:-1]
                    foldername = "_".join(foldername)
                filepath = glob.glob("{}.pdb".format(os.path.join(simulation_results_path, foldername, foldername,
                                                              path_to_pdbs, ID)))
                print(filepath)
                shutil.copy(filepath[0], os.path.join(pdbs_output_path, "{}_{}_{}_{}.pdb".format(foldername, ID, column,
                                                                                                 criteria)))


if __name__ == '__main__':
    path_results, column, criteria, output_filepath, simulation_path, path_to_pdbs, id_name, folder_name, \
    pdbs_output_path, xtc = parse_arguments()
    main(path_results, column, criteria, output_filepath, simulation_path, path_to_pdbs, id_name, folder_name,
         pdbs_output_path, xtc=xtc)

