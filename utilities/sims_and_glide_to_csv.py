import os
import sys
import glob
import argparse
import pandas as pd
import variables as vrb


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""This program collects information from simulations and glide reports
                                                    and creates a single csv with everything merged for each folder.""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-o", "--out_path", required=True,
                                help="""Absolute path where the final reports will be stored.""")
    required_named.add_argument("-sp", "--sim_path", required=True,
                                help="""Absolute path where the final reports will be stored.""")
    required_named.add_argument("-gp", "--glide_path", required=True,
                                help="""Absolute path where the final reports will be stored.""")

    parser.add_argument("-spt", "--sims_pattern", default=vrb.PATH_PATTER_TO_CLUSTERS_FROM_SIMS,
                        help="""""")
    parser.add_argument("-gpt", "--glide_pattern", default=vrb.GLIDE_CSV_PATTERN,
                        help="""""")
    parser.add_argument("-rpt", "--report_pattern", default=vrb.REPORT_FILE_PATTERN,
                        help="""""")

    args = parser.parse_args()

    return args.out_path, args.sim_path, args.glide_path, args.sims_pattern, args.glide_pattern, args.report_pattern


def get_folders_in_directory(abs_path_to_folder):
    folders = glob.glob("{}/*".format(abs_path_to_folder))
    return folders


def find_common_folders(simulation_path_results, glide_results_path):
    """
    Given two paths: 1) path where the simulations and 2) path where the glide's results are stored, it returns which
    folders are common in order to relate each simulation with its glide.
    :param simulation_path_results: path of the simulations results
    :param glide_results_path: path with glide's results
    :return: name of folders that are common in both folders (that came from the same ligand)
    """
    sim_folders = get_folders_in_directory(simulation_path_results)
    glide_folders = get_folders_in_directory(glide_results_path)
    sim_folder_names = [os.path.basename(folder) for folder in sim_folders]
    glide_folder_names = [os.path.basename(folder) for folder in glide_folders]
    if set(sim_folder_names) != set(glide_folder_names):
        print("WARNING: SOME FOLDERS ARE DIFFERENT!")
        common_folders = list(set(sim_folder_names).intersection(glide_folder_names))
    else:
        common_folders = sim_folder_names
    return common_folders


def find_simulations_and_glide_reports(simulation_path_results, glide_results_path,
                                       simulations_pattern_path=vrb.PATH_PATTER_TO_CLUSTERS_FROM_SIMS,
                                       glide_pattern=vrb.GLIDE_CSV_PATTERN,
                                       report_cluster_name_pattern=vrb.REPORT_FILE_PATTERN):
    """
    This function finds for each ligand or simulation folder its glide results, returning two list; one list which
    contains a the path to the report and another path to the glide's report, and the second list that contains the main
    folder names.
    :param simulation_path_results: path where the simulations are stored.
    :param glide_results_path: path where the glide's results are stored.
    :param simulations_pattern_path: pattern that will be used to find the report file starting from the
    "simulation_path_results".
    :param glide_pattern: pattern to find the file where the report of glide with the results.
    :param report_cluster_name_pattern: suffix or prefix (pattern) that follows the report file of the simulations.
    :return: list with the first element; a list with lists of two elements (element1: path to simulation's report,
    element2: path to glide's report), and the second element that is a list with the main folders names (or ligands).
    """
    files_to_analyze = []
    ligands_to_analyze = []
    print("SIMULATIONS FOLDER:{}\nGLIDE's FOLDER:{}".format(simulation_path_results, glide_results_path))
    foldernames = find_common_folders(simulation_path_results, glide_results_path)
    for folname in foldernames:
        simulation_pattern_to_clusters_report = os.path.join(simulation_path_results, folname, folname,
                                                          simulations_pattern_path,
                                                          report_cluster_name_pattern)
        simulation_path_to_clusters_report = glob.glob(simulation_pattern_to_clusters_report)
        if len(simulation_path_to_clusters_report) > 1:
            exit("""More than one folder is being detected for simulations {}, check if patters are correct. Folders found:\n{}""".format(folname, simulation_pattern_to_clusters_report))

        glide_pattern_to_report = os.path.join(glide_results_path, folname, glide_pattern)
        glide_path_to_report = glob.glob(glide_pattern_to_report)
        if len(glide_path_to_report) > 1:
            exit("""More than one folder is being detected for glides {}, check if patters are correct. Folders found:\n{}""".format(folname, glide_path_to_report))

        files_to_analyze.append([simulation_path_to_clusters_report, glide_path_to_report])
        ligands_to_analyze.append(folname)
    print("FOLDERS TO ANALYZE:\n")
    for name in ligands_to_analyze:
        print("{}".format(name))
    return files_to_analyze, ligands_to_analyze


def report_to_dataframe(report_path, separator="\s+"):
    """
    Given a path to a report file it returns a pandas dataframe object.
    :param report_path: path to the report file
    :return: dataframe pandas object
    """
    dataframe = pd.read_csv(report_path, sep=separator, engine='python', header=0)
    return dataframe


def glide_to_dataframe(glide_results_path):
    """
    Given a path to a glide results csv it returns a pandas dataframe object.
    :param glide_results_path: path to the csv file with glide scores
    :return: dataframe pandas object
    """
    dataframe = pd.read_csv(glide_results_path, sep=',', engine='python')
    return dataframe


def get_index_of_pdb_clusterized(pdb_name):
    """
    Given the name of a pdb resulting of the clustering process with the following format name: "cluster_0_processed" it
    returns the number between underscores in order to get the real index.
    :param pdb_name: string with the name of the pdb file
    :return: index between underscores
    """
    name_to_list_splited_by_underscores = pdb_name.split("_")
    index = name_to_list_splited_by_underscores[1]
    return index


def add_index_column_to_glide_df(glide_dataframe):
    """
    This function adds to a glide's dataframe a column with the index obtained of the ID column in order to latter trace
    the information of the simulation's report file.
    :param glide_dataframe: dataframe of pandas with the information of glide scores
    :return: the same dataframe with the 'index' column added
    """
    index_list = []
    for ID in glide_dataframe['ID']:
        index = get_index_of_pdb_clusterized(ID)
        index_list.append(int(index))
    glide_dataframe['index'] = index_list
    return glide_dataframe


def add_glide_information_to_report(report_dataframe, glide_dataframe):
    """
    Given a dataframe with the information from simulations and another with the information from glide this function
    merge both into a single dataframe with all the information together.
    :param report_dataframe: dataframe pandas object from simulations
    :param glide_dataframe: dataframe pandas object from glide results
    :return: single dataframe with all information together
    """
    report_header = list(report_dataframe)
    glide_header = list(glide_dataframe)
    joined_headers = report_header + glide_header
    sorted_glide_df = glide_dataframe.sort_values('index')
    sorted_glide_df = sorted_glide_df.reset_index(drop=True)  # Sort glide df by the index of the report df
    result = pd.concat([report_dataframe, sorted_glide_df], axis=1, ignore_index=True)
    result.columns = joined_headers
    return result


def write_dataframes_to_csv_file(dataframes, foldernames, output_folder):
    """
    Given a list of dataframes it writes csv files and store them into a output_folder naming them using foldernames.
    :param dataframes: list of pandas dataframes
    :param foldernames: lisf ot folder names (or ligand names)
    :param output_folder: path to the output file
    :return: writes the csv files
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    for foldername, dataframe in zip(foldernames, dataframes):
        dataframe.to_csv(os.path.join(output_folder, "{}.csv".format(foldername)), sep=';', index=False)


def main(simulation_path_results, glide_results_path, simulations_pattern_path, glide_pattern,
         report_cluster_name_pattern, output_path):
    dataframes = []
    paths_to_simulations_and_glide, ligands = find_simulations_and_glide_reports(simulation_path_results,
                                                                                 glide_results_path,
                                                                                 simulations_pattern_path,
                                                                                 glide_pattern,
                                                                                 report_cluster_name_pattern)
    for paths in paths_to_simulations_and_glide:
        simulation_report_path, glide_report_path = paths
        sim_dataframe = report_to_dataframe(simulation_report_path[0])  # We will use the index 0 because it is a list
        glide_dataframe = glide_to_dataframe(glide_report_path[0])
        glide_dataframe = add_index_column_to_glide_df(glide_dataframe)
        merged_dataframe = add_glide_information_to_report(sim_dataframe, glide_dataframe)
        dataframes.append(merged_dataframe)
    write_dataframes_to_csv_file(dataframes, ligands, output_path)


if __name__ == '__main__':
    out_path, sim_path, glide_path, sims_pattern, glide_pattern, report_pattern = parse_arguments()
    main(sim_path, glide_path, sims_pattern, glide_pattern, report_pattern, out_path)
