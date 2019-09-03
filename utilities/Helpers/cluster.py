import pandas as pd
import numpy as np
import shutil
import os
import matplotlib.pyplot as plt
import glob
import multiprocessing as mp
from sklearn.cluster import KMeans


def cluster_df(dataframe, cluster_by_list):
    df = dataframe.filter(cluster_by_list, axis=1)
    kmean = KMeans(n_clusters=500).fit(df)
    col = df.columns.values
    new_df_km = pd.DataFrame(kmean.cluster_centers_, columns=col)
    return new_df_km


def find_minimum_for_row_in_df(dataframe, row_c, cluster_by_columns):
    vector_A = np.array([row_c[name] for name in cluster_by_columns])
    df_distances = []
    for index, row in dataframe.iterrows():
        vector_B = np.array([row[name] for name in cluster_by_columns])
        distance = np.linalg.norm(np.linalg.norm(vector_A) - np.linalg.norm(vector_B))
        df_distances.append([distance, index])
    minimum = min(df_distances)
    print(dataframe.iloc[minimum[1]])
    return dataframe.iloc[minimum[1]]


def detect_row_close_to_centroid(dataframe, centroids, cluster_by_columns):
    pool = mp.Pool(processes=4)
    multi = []
    new_df = pd.DataFrame()
    for index_c, row_c in centroids.iterrows():
        multi.append(pool.apply_async(find_minimum_for_row_in_df, [dataframe, row_c, cluster_by_columns]))
    for result in multi:
        new_df = new_df.append(result.get())
    pool.terminate()
    return new_df


def extract_structure_from_pdb_folder(path_to_pdbs, dataframe, column_trajectory="file_from",
                                      model_column="numberOfAcceptedPeleSteps"):
    pdbs = glob.glob("{}/*".format(path_to_pdbs))
    print(pdbs)
    new_path = os.path.join(path_to_pdbs, "clustering")
    print(new_path)
    if not os.path.exists(new_path):
        os.mkdir(new_path)
    for index, row in dataframe.iterrows():
        filename = row[column_trajectory].split("/")[-1].split(".xtc")[0]
        epoch = row[column_trajectory].split("/")[-2]
        model = row[model_column]
        file = os.path.join(path_to_pdbs, "{}_epoch_{}_snap_{}.pdb".format(filename, int(epoch), int(model)))
        if file in pdbs:
            shutil.copy(file, os.path.join(new_path, file.split("/")[-1]))
        else:
            print("File {} does not exist!".format(file))


def main(csv_path, out_path=None, sep=";", path_to_pdbs=None,
         clusterization_by_columns=("Binding Energy", "RMSD_ligand", "currentEnergy", "crystal_com_dist",
                                    "total_hbonds", "spec_hbonds")):
    df = pd.read_csv(csv_path, sep=sep)
    df = df.fillna(0)
    df = df.reset_index()
    print(df.columns)
    print(clusterization_by_columns)
    km = cluster_df(df, clusterization_by_columns)
    if not out_path:
        out_path = csv_path.split(".csv")[0]+"_clustering.csv"
    print(out_path)
    new_df = detect_row_close_to_centroid(df, km, clusterization_by_columns)
    if path_to_pdbs:
        extract_structure_from_pdb_folder(path_to_pdbs, new_df, column_trajectory="file_from")
    new_df.to_csv(out_path, sep=sep, index=False)


csv_list = glob.glob("/home/carlespl/project/Almirall/validation_compounds/validation_compounds_2/conf_12/results_summary/*.csv")
pdbs_paths = glob.glob("/home/carlespl/project/Almirall/validation_compounds/validation_compounds_2/conf_12/temporary/*/")

for csv in csv_list:
    print(csv)
    csv_id = csv.split("/")[-1].split("_adaptive_results_summary.csv")[0]
    print(csv_id)
    for pdb in pdbs_paths:
        print(pdb)
        pdb_id = pdb.split("/")[-2]
        print(pdb_id)
        if pdb_id == csv_id:
            main(csv_path=csv, path_to_pdbs=pdb)


