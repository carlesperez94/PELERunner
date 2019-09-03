import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np
import variables as vrb
import mdtraj as md
import prody
import multiprocessing as mp
from Helpers.structure_handler import open_structure_from_dataframe
from Helpers.structure_handler import read_structure_pdb
from Helpers.structure_handler import read_structure_xtc
from Helpers.structure_handler import select_atom_pair_from_ligand_and_residue_mdtraj
from Helpers.structure_handler import select_atom_pair_from_ligand_and_residue
from Helpers.dataframe_handler import add_column_to_dataframe


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""This program compute atom atom distances
                                     for AdaptivePELE simulations.""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-csv_file",
                                help="""""")
    required_named.add_argument("-path_to_simulations",
                                help="""""")
    required_named.add_argument("-resid",
                                help="""""")
    required_named.add_argument("-atom_residue",
                                help="""""")
    required_named.add_argument("-atom_ligand",
                                help="""""")
    parser.add_argument("-cn", "--column_name", default="distance",
                                help="""""")
    parser.add_argument("-p2c", "--path_to_clusters", default=vrb.PATH_PATTER_TO_CLUSTERS_FROM_SIMS,
                        help="""""")
    parser.add_argument("-lr","--ligand_resname", default="LIG",
                                help="""""")
    parser.add_argument("-fc","--file_column", default="file_from",
                                help="""""")
    parser.add_argument("-s", "--separator", default=";",
                        help=""".""")
    parser.add_argument("-lc", "--ligand_chain", default="L",
                        help="""""")
    parser.add_argument("-xtc", "--read_xtc", default=False,
                        help="""""")

    args = parser.parse_args()

    return args.csv_file, args.path_to_simulations, args.resid, args.atom_residue, args.atom_ligand, args.column_name, \
           args.path_to_clusters, args.file_column, args.separator, args.ligand_chain, args.ligand_resname, \
           args.read_xtc


def calculate_distance_of_to_atoms(atom1, atom2):
    coords_atom1 = atom1.getCoords()
    coords_atom2 = atom2.getCoords()
    distance = prody.calcDistance(coords_atom1, coords_atom2)
    return distance[0]


def select_atom_of_ligand_and_residue_and_get_distance(pdb_file, resid, name_residue, name_ligand, ligand_chain="Z"):
    atom_ligand, atom_residue = select_atom_pair_from_ligand_and_residue(pdb_file, resid, name_residue, name_ligand,
                                                                         ligand_chain)
    distance = calculate_distance_of_to_atoms(atom_ligand, atom_residue)
    return distance


def compute_distances_mdtraj(trajectory, resid, name_residue, name_ligand, ligand_resname="LIG"):

    atom_ligand, atom_residue = select_atom_pair_from_ligand_and_residue_mdtraj(trajectory, resid, name_residue,
                                                                                name_ligand, ligand_resname)
    atoms = np.ndarray(buffer=np.array([atom_ligand, atom_residue]), shape=(1, 2), dtype=int)
    distances = md.compute_distances(trajectory, atoms)  # distances in nm
    distances_amstrongs = [float(x)*10 for x in distances]
    return distances_amstrongs


def compute_distances_mdtraj_from_file(file_path, resid, name_residue, name_ligand, ligand_resname="LIG"):
    topology_path_splited = file_path.split("/")[0:-2]
    topology_path = os.path.join("/".join(topology_path_splited), "topology.pdb")
    trajectory = read_structure_xtc(file_path, topology_path)
    distances = compute_distances_mdtraj(trajectory, resid, name_residue, name_ligand, ligand_resname)
    for distance in distances:
        print(distance)
    return distances


def main(csv_file, path_to_simulations, resid, name_residue, name_ligand, column_name="distance",
         path_to_clusters=vrb.PATH_PATTER_TO_CLUSTERS_FROM_SIMS, file_column="file_from", separator=";",
         ligand_chain="Z", ligand_resname="LIG", read_xtc=False):
    files_to_analyze = open_structure_from_dataframe(csv_file=csv_file, path_to_simulations=path_to_simulations,
                                                     path_to_clusters=path_to_clusters, file_column=file_column,
                                                     separator=separator, read_xtc=read_xtc)
    distances_list = []
    multiprocessing_list = []
    pool = mp.Pool(processes=4)
    for n, file in enumerate(files_to_analyze):
        if not read_xtc:
            distance = select_atom_of_ligand_and_residue_and_get_distance(file, resid, name_residue, name_ligand, ligand_chain)
            distances_list.append(distance)
        else:
            multiprocessing_list.append(pool.apply_async(compute_distances_mdtraj_from_file,
                                                 (file, resid, name_residue, name_ligand, ligand_resname)))
    for process in multiprocessing_list:
        distances_list.extend(process.get())
    dataframe = pd.read_csv(csv_file, sep=separator)
    add_column_to_dataframe(distances_list, column_name, dataframe)
    dataframe.to_csv(os.path.splitext(csv_file)[0] + "recomputed.csv", sep=separator, index=False)
    return dataframe


if __name__ == '__main__':
    csv_file, path_to_simulations, resid, atom_residue, atom_ligand, column_name, path_to_clusters, file_column, \
    separator, ligand_chain, ligand_resname, read_xtc = parse_arguments()
    main(csv_file, path_to_simulations, resid, atom_residue, atom_ligand, column_name, path_to_clusters, file_column,
         separator, ligand_chain, ligand_resname, read_xtc)
