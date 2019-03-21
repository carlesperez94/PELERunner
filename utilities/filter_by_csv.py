import os
import sys
import argparse
import pandas as pd
from sims_and_glide_to_csv import report_to_dataframe
from summarize_all_csvs_in_one import add_filename_to_csv_as_column


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""Given a csv file it returns a new csv filtering the dataset by
                                                    a certain column following the selected criteria.""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-i", "--input_path", required=True,
                                help="""Absolute path to the csv file.""")
    required_named.add_argument("-cr", "--criteria", required=True, choices=['greater', 'greater-equal', 'lower',
                                                                             'lower-equal', 'equal', 'different'],
                                help="""Criteria that will be used in order to do a subset of the dataset. Example: if 
                                you would like to get all values greater than 3, select: '>' [value=3].""")
    required_named.add_argument("-v", "--value", required=True,
                                help="""Value that will set the threshold to filter by.""")
    required_named.add_argument("-c", "--column", required=True,
                                help="""Column of the dataframe that will be used as variable to do the subset.""")
    required_named.add_argument("-o", "--output_folder", required=True,
                                help="""Path to the output folder.""")

    args = parser.parse_args()

    return args.input_path, args.criteria, args.value, args.column, args.output_folder


def filter_by(csv_file, criteria, value, column_to_filter, output_folder):
    df = report_to_dataframe(csv_file, separator=";")
    if criteria == "greater":
        df = df[df[column_to_filter] > float(value)]
    if criteria == "greater-equal":
        df = df[df[column_to_filter] >= float(value)]
    if criteria == "lower":
        df = df[df[column_to_filter] < float(value)]
    if criteria == "lower-equal":
        df = df[df[column_to_filter] <= float(value)]
    if criteria == "equal":
        df = df[df[column_to_filter] == value]
    if criteria == "different":
        df = df[df[column_to_filter] != value]
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    add_filename_to_csv_as_column(df, os.path.basename(csv_file).split(".")[0], column_name="filename")
    if len(df) > 0:
        df.to_csv(os.path.join(output_folder, "{}_filtered_by_{}_{}.csv".format(os.path.basename(csv_file).split(".")[0],
                                                                          column_to_filter, value)), sep=';',
                                                                          index=False)
    else:
        print("{} has no values of {} {} than {}".format(csv_file, column_to_filter, criteria, value))


if __name__ == '__main__':
    input_path, criteria, value, column, output_folder = parse_arguments()
    filter_by(input_path, criteria, value, column, output_folder)
