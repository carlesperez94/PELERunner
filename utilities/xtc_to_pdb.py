import sys
import multiprocessing as mp
import os
import glob
import pandas as pd
from AdaptivePELE.atomset import atomset
import AdaptivePELE.utilities.utilities as adapt_tools


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


def get_pdb_from_xtc(row, pdbs_output_path, column_file="trajectory"):
    """
    Given a row of a dataframe (expected to come from a csv report) and a column name (that must contain the path to
    its correspondent trajectory), this function extract the file in PDB format in an output file.
    :param row: row of a dataframe (Pandas object).
    :type row: pandas.DataFrame
    :param pdbs_output_path: output path for the PDB file.
    :type pdbs_output_path: str
    :param column_file: Column name of the dataframe that contains the path to the trajectory file.
    :type column_file: str
    :return:
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


def concat_reports_in_csv(adaptive_results_path, output_file_path, report_prefix="report_",
                          trajectory_prefix="trajectory_", separator_out=";"):
    """
    It search report files in Adaptive's result folder and creates a csv file with everything concatenated, adding the
    epoch and trajectory information.
    :param adaptive_results_path: Path to the results folder of Adaptive.
    :type adaptive_results_path: str
    :param output_file_path: Path of the output file.
    :type output_file_path: str
    :param report_prefix: Prefix of PELE's reports.
    :type report_prefix: str
    :param trajectory_prefix: Prefix of PELE's trajectories.
    :type trajectory_prefix: str
    :param separator_out: Separator string used in the csv file.
    :type separator_out: str
    :return: Creates a csv file.
    """
    dataframe_lists = []
    for adaptive_epoch in range(0, 2000):
        folder = os.path.join(adaptive_results_path, str(adaptive_epoch))
        if os.path.exists(folder):
            report_list = glob.glob("{}/*{}*".format(folder, report_prefix))
            report_list = sorted(report_list, key=lambda x: int(x.split("_")[-1]))
            for n, report in enumerate(report_list):
                pandas_df = pd.read_csv(report, sep="    ", engine="python", index_col=False, header=0)
                pandas_df["epoch"] = adaptive_epoch
                pandas_df["trajectory"] = glob.glob("{}/{}/*{}{}.*".format(adaptive_results_path, adaptive_epoch,
                                                                          trajectory_prefix, n + 1))[0]
                dataframe_lists.append(pandas_df)
        else:
            break
    dataframe = pd.concat(dataframe_lists, ignore_index=True)
    dataframe.to_csv(output_file_path, sep=separator_out, index=False)


def get_pdbs_from_df_in_xtc(df, pdbs_output_path, processors=4, column_file="trajectory"):
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


def main(adaptive_results_path, output_folder_path, summary_done=False, separator=";",
         report_pref="report_", trajectory_pref="trajectory_", processors=4, column_file="trajectory"):
    if not os.path.exists(output_folder_path):
        os.mkdir(output_folder_path)
    summary_csv_filename = os.path.join(output_folder_path, "summary.csv")
    if not summary_done:
        concat_reports_in_csv(adaptive_results_path=adaptive_results_path, output_file_path=summary_csv_filename,
                              report_prefix=report_pref, trajectory_prefix=trajectory_pref, separator_out=separator)
    dataframe = pd.read_csv(summary_csv_filename, sep=separator, engine='python', header=0)
    multiprocessing = []
    get_pdbs_from_df_in_xtc(dataframe, output_folder_path, processors=processors, column_file=column_file)
    for process in multiprocessing:
        process.get()


if __name__ == '__main__':
    adaptive_results_path = sys.argv[1]
    output_folder_path = sys.argv[2]
    main(adaptive_results_path, output_folder_path)

