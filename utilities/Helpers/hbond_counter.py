import os
import glob
import traceback
import pandas as pd
import multiprocessing
from collections import Counter


def count_total_hbonds_for_model(csv_path, special_bonds_to_count, sep=",", model_col="model",
                                 trajectory_col="trajectory", hbonds_col="hbonds"):
    dataframe = pd.read_csv(csv_path, sep=sep, header=0)
    list_of_traj_and_models = []
    list_of_traj_and_models_with_spec_hbonds = []
    for index, row in dataframe.iterrows():
        trajectory_name = row[trajectory_col]
        model = row[model_col]
        hbond = row[hbonds_col]
        print(trajectory_name, model, hbond)
        list_of_traj_and_models.append("{}={}".format(trajectory_name, model))
        for hbond_spec in special_bonds_to_count:
            if hbond_spec in hbond:
                list_of_traj_and_models_with_spec_hbonds.append("{}={}".format(trajectory_name, model))
    dictionary = dict(Counter(list_of_traj_and_models))
    dictionary_spec = dict(Counter(list_of_traj_and_models_with_spec_hbonds))
    final_df = pd.DataFrame(columns=["trajectory", "model", "total_hbonds", "spec_hbonds"])
    for key, value in dictionary.items():
        trajectory_path, model = key.split("=")
        try:
            df_tmp = pd.DataFrame([[trajectory_path, model, value, int(dictionary_spec[key])]],
                                  columns=["trajectory", "model", "total_hbonds", "spec_hbonds"])
        except KeyError:
            df_tmp = pd.DataFrame([[trajectory_path, model, value, int(0)]],
                                  columns=["trajectory", "model", "total_hbonds", "spec_hbonds"])
        final_df = final_df.append(df_tmp, ignore_index=True)
    output_path = os.path.join("/".join(csv_path.split("/")[0:-1]), "hbond_summary.csv")
    final_df.to_csv(output_path, index=False)


def add_report_column_to_csv(csv_to_add_path, report_col_name="Binding Energy"):
    simulation_path = "/".join(csv_to_add_path.split("/")[0:-1])
    simulation_reports = glob.glob("{}/[0-9]*/*report_*".format(simulation_path))
    hbond_df = pd.read_csv(csv_to_add_path, sep=",", header=0, engine="python")
    hbond_df[report_col_name] = "NaN"
    for report in simulation_reports:
        trajectory_file = report.replace("report", "trajectory")+".xtc"
        traj_index = hbond_df.index[hbond_df["trajectory"] == trajectory_file]
        df = pd.read_csv(report, sep="    ", header=0, engine="python")
        for n, index in enumerate(traj_index):
            hbond_df.loc[index, report_col_name] = df[report_col_name][n]
    hbond_df.to_csv(os.path.join(simulation_path, "hbond_summary_completed.csv"), index=False)


def add_csv_column_to_csv(csv_original, csv_to_add_path, sep_o=",", sep_a=",", col_names_to_add=("total_hbonds",
                                                                                                 "spec_hbonds"),
                          col_name_to_id_original="trajectory", col_name_to_model_original="numberOfAcceptedPeleSteps",
                          col_name_to_id_add="file_from", col_name_to_model_add="model"):
    dt_original = pd.read_csv(csv_original, sep=sep_o, engine="python")
    print(csv_to_add_path)
    dt_to_add = pd.read_csv(csv_to_add_path, sep=sep_a, engine="python")
    for index, row in dt_original.iterrows():
        ID = row[col_name_to_id_original].split("/")[-2:]
        model_original = row[col_name_to_model_original]
        for index_add, row_add in dt_to_add.iterrows():
            ID_add = row_add[col_name_to_id_add].split("/")[-2:]
            model_add = row_add[col_name_to_model_add]
            if ID == ID_add and model_add == model_original:
                for col in col_names_to_add:
                    print(ID, model_add)
                    data_to_add = dt_to_add.iloc[index_add][col]
                    if not data_to_add:
                        data_to_add = 0
                    print(data_to_add)
                    dt_original.loc[index, col] = data_to_add
    dt_original.to_csv(csv_original, sep_o, index=False)




#folder_list = glob.glob("/home/carlespl/project/Almirall/validation_compounds/validation_compounds_2/conf_12/7_BMS986165_no_wat_2*/"
#                        "obc_adaptive*/hbonds_analysis.csv")
#
#for folder in folder_list:
#    count_total_hbonds_for_model(folder, special_bonds_to_count=["690-N", "VAL690-O", "GLU598-O"], sep=",",
#                                 model_col="model", trajectory_col="trajectory", hbonds_col="hbonds")
#
folder_list = glob.glob("/home/carlespl/project/Almirall/validation_compounds/validation_compounds_1/simulation_results_with_exp_data/*/obc_adaptive*/hbond_summary.csv")
print(folder_list)
csv_list = glob.glob("/home/carlespl/project/Almirall/validation_compounds/to_ML/*/selection_report.csv")
print(csv_list)
for folder in folder_list:
   #print(folder)
    folder_id = folder.split("/")[-3].split("_")[-1]
    print(folder_id)
    for csv in csv_list:
        csv_id = csv.split("/")[-2]
        print(csv_id)
        if folder_id == csv_id:
            add_csv_column_to_csv(csv_original=csv, csv_to_add_path=folder, col_name_to_id_original="trajectory",
                                  col_name_to_id_add="trajectory", sep_o=";")

            print(folder, csv)

