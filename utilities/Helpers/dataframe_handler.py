import sys
import pandas as pd


def add_column_to_dataframe(data, column_name, dataframe_to_be_added):
    dataframe_to_be_added[column_name] = data
    return dataframe_to_be_added
