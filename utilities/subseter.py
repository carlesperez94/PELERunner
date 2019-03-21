import os
import sys
import glob
import argparse
import multiprocessing as mp
import pandas as pd
import mdtraj as md
import matplotlib.pyplot as plt


def find_minimum_in_area(df, x_threshold=2.5, column="Binding Energy", x_column="hbond_H_val_690"):
    df = df.loc[df[x_column] < x_threshold ]
    minimum = df.loc[df[column] == df[column].min()]
    if len(minimum) > 1:
       return minimum.iloc[[0]]
    else:
       return minimum


def trajectory_and_snapshot_to_pdb(trajectory_path, snapshot, output_path):
    topology_path_splited = trajectory_path.split("/")[0:-2]
    topology_path = os.path.join("/".join(topology_path_splited), "topology.pdb")
    trajectory = md.load(trajectory_path, top=topology_path)
    try:
        single_model = trajectory[snapshot]
    except IndexError:
        exit("You are selecting the model {} for a trajectory that has {} models, please, reselect the model index (starting from 0).".format(snapshot, len(trajectory)))
    single_model.save_pdb(output_path)


def get_pdb_from_xtc(row, pdbs_output_path, column_file="file_from"):
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


def get_pdbs_from_df_in_xtc(df, pdbs_output_path, column_file="file_from"):
    pool = mp.Pool(processes=4)
    multiprocessing_list = []
    for index, row in df.iterrows():
        multiprocessing_list.append(pool.apply_async(get_pdb_from_xtc,
                                                     (row, pdbs_output_path, column_file)))
    for process in multiprocessing_list:
        process.get()


def subset_csv(csv_file, pdb_out_folder, x_column="hbond_H_val_690", y_column="Binding Energy", x_threshold=4,
               y_threshold=0, sep=";", definition=0.05, get_pdbs=False):
    df_original = pd.read_csv(csv_file, sep=sep)
    df = df_original.loc[df_original[x_column] <= x_threshold] # Selecting subset cutting off by x
    df = df.loc[df[y_column] <= y_threshold] # Selecting subset cutting off by y on the previous selection
    minimum = find_minimum_in_area(df, x_threshold=2.5, column=y_column, x_column=x_column)
    x_min = minimum[x_column].min()
    y_min = minimum[y_column].min()
    df_until_minimum = df.loc[df[x_column] <= x_min]
    m = (y_threshold - y_min)/(x_threshold - x_min) # Pendent of the function from the minimum to the (x_treshold, y_threshold)
    iterations = (x_threshold - x_min) / definition
    for n in range(1, int(round(iterations))+1):
        y_lim = y_min+n*m*definition
        x_lim = x_min+(definition*n)
        x = x_min+(definition*(n-1))
        df_new = df[(df[x_column] > x) & (df[x_column] <= x_lim)]
        df_new = df_new[(df_new[y_column] >= y_lim) & (df_new[y_column] < y_threshold)]
        df_until_minimum = df_until_minimum.append(df_new, ignore_index=True)
    plt.scatter(x=df_original[x_column], y=df_original[y_column])
    plt.scatter(x=df_until_minimum[x_column], y=df_until_minimum[y_column])
    plt.ylim(-70, 0)
    plt.xlabel(x_column)
    plt.ylabel(y_column)
    plt.title(os.path.basename(csv_file))
    plt.savefig("{}_subset.png".format(csv_file.split(".csv")[0]))
    df_until_minimum.to_csv("{}_subset.csv".format(csv_file.split(".csv")[0]), sep=sep, index=False)
    if get_pdbs:
        get_pdbs_from_df_in_xtc(df_until_minimum, pdbs_output_path=pdb_out_folder)


#subset_csv("/home/carlespl/project/Almirall/crystals/with_waters/summary_results/8_PKD_TYK2_SD006_with_waters_adaptive_results_summary.csv",
#           pdb_out_folder="/home/carlespl/project/Almirall/crystals/with_waters/summary_results/subset_to_machine_learning/8_PKD_TYK2_SD006",
#           y_threshold=-28, get_pdbs=True)
