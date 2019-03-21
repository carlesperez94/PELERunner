import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np
import mdtraj as md
import prody

import variables as vrb


def open_structure_from_dataframe(csv_file, path_to_simulations, path_to_clusters=vrb.PATH_PATTER_TO_CLUSTERS_FROM_SIMS,
                                  file_column="ID", separator=";", read_xtc=False):
    files_to_analyze = []
    dataframe = pd.read_csv(csv_file, sep=separator)
    if not read_xtc:
        for index, row in dataframe.iterrows():
            file_to_analyse = row[file_column]
            compound_name = os.path.basename(csv_file).split(".")[0]
            path_to_file = os.path.join(path_to_simulations, compound_name, compound_name, path_to_clusters, "{}.*".format(file_to_analyse))
            file_in = glob.glob(path_to_file)[0]
            files_to_analyze.append(file_in)
        return files_to_analyze
    else:
        for index, row in dataframe.iterrows():
            file_to_analyse = row[file_column]
            file_in = glob.glob(file_to_analyse)[0]
            if not file_in in files_to_analyze:
                files_to_analyze.append(file_in)
        return files_to_analyze


def read_structure_pdb(file_path):
    molecule = prody.parsePDB(file_path)
    return molecule


def read_structure_xtc(file_path, topology_file_path=None, pdb=False):
    if not pdb:
        trajectory = md.load(file_path, top=topology_file_path)
    else:
        trajectory = md.load(file_path)
    return trajectory


def select_atom_pair_from_ligand_and_residue_mdtraj(trajectory, resid, name_residue, name_ligand,
                                                    ligand_resname="LIG"):
    topology = trajectory.topology
    selection_string_residue = "resSeq {} and name {}".format(resid, name_residue)
    selection_string_ligand = "resname {} and name {}".format(ligand_resname, name_ligand)
    selection_residue = topology.select(selection_string_residue)
    selection_ligand = topology.select(selection_string_ligand)
    return selection_ligand, selection_residue


def select_atom_pair_from_ligand_and_residue(pdb_file, resid, name_residue, name_ligand, ligand_chain="Z"):
    molecule = read_structure_pdb(pdb_file)
    selection_ligand = molecule.select("chain {} and name {}".format(ligand_chain, name_ligand))
    selection_residue = molecule.select("resid {} and name {}".format(resid, name_residue))
    return selection_ligand, selection_residue


def select_ligand_mdtraj(trajectory, ligand_resname="LIG", heavy=True):
    topology = trajectory.topology
    if heavy:
        selection_string_ligand = "resname '{}' and (not type H)".format(ligand_resname)
    else:
        selection_string_ligand = "resname '{}'".format(ligand_resname)
    selection_ligand = topology.select(selection_string_ligand)
    return selection_ligand


def select_ligand_and_residues(trajectory, list_of_residues_ids, ligand_resname="LIG", heavy=True):
    topology = trajectory.topology
    if heavy:
        selection_string = "(resname '{}' or resSeq {}) and (not type H)".format(ligand_resname, list_of_residues_ids)
    else:
        selection_string = "resname '{}' or resSeq {}".format(ligand_resname, list_of_residues_ids)
    selection_ligand_and_residues = topology.select(selection_string)
    return selection_ligand_and_residues

