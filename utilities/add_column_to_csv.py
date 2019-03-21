import os
import sys
import glob
import argparse
import pandas as pd
from sims_and_glide_to_csv import report_to_dataframe


def add_filename_to_csv_as_column(dataframe, filepath, column_name="origin_name"):
    filename = os.path.basename(filepath).split(".")[0]
    dataframe[column_name] = filename
    return dataframe


def add_column_to_csvs_in_folder(in_folder, column_added_name="compound_name"):
    files = glob.glob("{}/*.csv".format(in_folder))
    for file in files:
        compound_name = os.path.basename(file).split(".")[0]
        dataframe = report_to_dataframe(file, separator=";")
        dataframe_with_column = add_filename_to_csv_as_column(dataframe, compound_name, column_name=column_added_name)
        dataframe_with_column.to_csv(file, sep=';', index=False)
        print("{} completed!".format(file))
    return print("Done")


add_column_to_csvs_in_folder("/home/carlespl/project/Almirall/crystals/with_waters/filtered_results")