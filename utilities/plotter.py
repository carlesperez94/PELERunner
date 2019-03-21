import os
import sys
import argparse
import glob
import mdtraj as md
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import multiprocessing as mp
from sims_and_glide_to_csv import report_to_dataframe
from AdaptivePELE.atomset import atomset
import AdaptivePELE.utilities.utilities as adapt_tools


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    desc = "Generates a scatterplot of a csv given two or three columns (X, Y, and Z if set).\n" \
           "This plot has a 'selection mode', where the user can set a selection of desired points by drawing.\n" \
           "Structures will be selected and stored into an output folder. Additionally, a report file of this selected\n" \
           "structures will be created. To be run for example like: \n" \
           "\">python plotter.py -i /home/usr/csvfile.csv -x 'Binding Energy' -y epoch -z NumberOfAcceptedStep -s -o figures_to_save\"\n"\
           "Exemple of selection: \n" \
           "\">python plotter.py -i /home/usr/csvfile.csv -x 'Binding Energy' -y epoch -z NumberOfAcceptedStep -slc -pdbc file_from\"\n"
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description=desc)
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-i", "--input_file", required=True, type=str,
                                help="""Absolute path of the csv file""")
    required_named.add_argument("-x", "--col_x", required=True, type=str,
                                help="""Column of the csv file that will be used in the X-axis.""")
    required_named.add_argument("-y", "--col_y", required=True, type=str,
                                help="""Column of the csv file that will be used in the Y-axis.""")
    required_named.add_argument("-z", "--col_z", required=True, type=str,
                                help="""Column of the csv that will be used as Z axis (gradient of color).""")
    
    parser.add_argument("-t", "--title", default=None, type=str,
                        help="""Title of the graph. If None, the title will be the basename of the csv file.""")
    parser.add_argument("-xl", "--x_lab", default=None, type=str,
                        help="""Name of the X axis label. If None, put the name of the col_x by default.""")
    parser.add_argument("-yl", "--y_lab", default=None, type=str,
                        help="""Name of the Y axis label. If None, put the name of the col_y by default.""")
    parser.add_argument("-zl", "--z_lab", default=None, type=str,
                        help="""Name of the Z axis label. If None, put the name of the col_z by default.""")
    parser.add_argument("-s", "--save", action="store_true",
                        help="""If True it saves a file in the output folder.""")
    parser.add_argument("-o", "--out_fold", default=None, type=str,
                        help="""Folder where the images will be saved.""")
    parser.add_argument("-slc", "--select", action="store_true",
                        help="If set you can select structures in the graph")
    parser.add_argument("-os", "--out_sel", default="", type=str,
                        help="Folder were the selection will be stored.")
    parser.add_argument("-pdbc", "--pdbs_column", default=None, type=str,
                        help="Column of the dataframe with the path of the structures.")

    args = parser.parse_args()

    return args.input_file, args.col_x, args.col_y, args.col_z, args.title, args.x_lab, args.y_lab, args.z_lab, \
           args.save, args.out_fold, args.select, args.out_sel, args.pdbs_column


class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    """

    def __init__(self, ax, collection, alpha_other=0.1):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


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


def create_scatter_plot(csv_file, column_to_x, column_to_y, column_to_z, title=None, x_label=None,
                        y_label=None, z_label=None, save_image=False, output_folder=None, select_mode=False,
                        output_selection_folder="", select_pdbs_column=None, processors=4, separator=";"):
    """
    Generates a scatterplot of a csv file given two or three columns (X, Y, and Z if set).
    This plot allows the selection of desired points by drawing. Structures will be selected and
    stored into an output folder. Additionally, a report file of this selected structures
    will be created.
    :param csv_file: Path to csv file that will be analyzed.
    :type csv_file: str
    :param column_to_x: Column name of the report file that will be used in the X axis.
    :type column_to_x: str
    :param column_to_y: Column name of the report file that will be used in the Y axis.
    :type column_to_y: str
    :param column_to_z: If set, column name of the report file that will be used in the Z axis (colorbar).
    :type column_to_z: str
    :param title: Title of the graph. If None, the title will be the basename of the file.
    :type title: str
    :param output_selection_folder: path to the output's folder.
    :type output_selection_folder: str
    :param x_label: Name of the X axis label. If None, put the name of the colum_to_x by default.
    :type x_label: str
    :param y_label: Name of the Y axis label. If None, put the name of the colum_to_y by default.
    :type y_label: str
    :param z_label: Name of the Z axis label. If None, put the name of the colum_to_z by default.
    :type z_label: str
    :param processors: Number of processors that you want to use in order to save time.
    :type processors: int
    :param save_image: If the user wants to save the graph this boolean must be set.
    :type save_image: bool
    :param output_folder: Path to the folder that to store the graphs.
    :type output_folder: str
    :param separator: Separator string that will be used in the CSV files.
    :type separator: str
    :param select_pdbs_column: Column name of the dataframe that contains the path to the trajectory file.
    :type select_pdbs_column: str
    :param select_mode: If it is set, the user will be able to draw lines on the plot to select specific structures.
    :type select_mode: bool
    :return:
    """
    df = report_to_dataframe(csv_file, separator=";")
    if select_mode:
        fig, ax = plt.subplots()
        pts = ax.scatter(df[column_to_x], df[column_to_y], c=df[column_to_z], s=20)
        selector = SelectFromCollection(ax, pts)

        def accept(event):
            if event.key == "enter":
                print("Selected points:")
                df_select = df.loc[selector.ind]
                print(df_select)
                if not os.path.exists(output_selection_folder):
                    os.mkdir(output_selection_folder)
                df_select.to_csv(os.path.join(output_selection_folder, "selection_report.csv"), sep=separator,
                                 index=False)
                if select_pdbs_column:
                    get_pdbs_from_df_in_xtc(df_select, output_selection_folder, column_file=select_pdbs_column,
                                            processors=processors)
                selector.disconnect()
                ax.set_title("")
                fig.canvas.draw()
        fig.canvas.mpl_connect("key_press_event", accept)
        ax.set_title("Press enter to accept selected points.")
        ax.set_xlabel(column_to_x)
        ax.set_ylabel(column_to_y)
        plt.show()
    else:
        if not title:
            title = os.path.basename(csv_file).split(".")[0]
        if not x_label:
            x_label = column_to_x
        if not y_label:
            y_label = column_to_y
        if not z_label:
            z_label = column_to_z
        plot_df = df.plot.scatter(column_to_x, column_to_y, c=column_to_z, title=title, colormap='viridis')
        plot_df.set_xlabel(x_label)
        plot_df.set_ylabel(y_label)
        if save_image:
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)
            plt.savefig(os.path.join(output_folder, "{}-{}-{}-{}.png".format(title, x_label, y_label, z_label)))
        else:
            plt.show()


if __name__ == '__main__':
    input_file, col_x, col_y, col_z, title, x_lab, y_lab, z_lab, save, output_folder, select, out_sel, pdbs_col = \
        parse_arguments()
    create_scatter_plot(input_file, col_x, col_y, col_z, title, x_lab, y_lab, z_lab, save, output_folder, select,
                        out_sel, pdbs_col)

