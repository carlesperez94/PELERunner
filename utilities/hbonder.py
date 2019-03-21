import os
import sys
import argparse
import traceback

import pandas as pd
import mdtraj as md
import numpy as np
import multiprocessing as mp
import md_traj_h_bond_modified as mdm



def extract_single_model_from_structure(trajectory_file_path, topology_file_path, model_number):
    trajectory = md.load(trajectory_file_path, top=topology_file_path)
    try:
        single_model = trajectory[model_number]
    except IndexError:
        exit("You are selecting the model {} for a trajectory that has {} models, please, reselect the model index (starting from 0).".format(model_number, len(trajectory)))
    return single_model


def select_residue(trajectory_object, residue):
    topology = trajectory_object.topology
    selection_string = "resSeq {}".format(residue)
    selection = topology.select(selection_string)
    return selection


def get_residues_given_dictionary(trajectory_object, residues_dictionary):
    residues = residues_dictionary.keys()
    results = []
    for residue in residues:
        selected_residue = select_residue(trajectory_object, residue)
        results.append(selected_residue)
    return results


def select_ligand(trajectory_object, ligand_resname="LIG"):
    topology = trajectory_object.topology
    selection_string = "resname {}".format(ligand_resname)
    selection = topology.select(selection_string)
    return selection


def print_atoms_of_selection(selection, trajectory):
    topology = trajectory.topology
    for atom in selection:
        print(topology.atom(atom))


def hunt_h_bonds_of_model(trajectory, cutoff=0.25):
    h_bonds = mdm.baker_hubbard(trajectory, cutoff)
    return h_bonds


def select_specific_h_bond_with_ligand(trajectory, residues_dict, lig_resname="LIG", cutoff=0.25):
    ligand = select_ligand(trajectory, lig_resname)
    residues_list = get_residues_given_dictionary(trajectory, residues_dict)
    residues_list = [j for i in residues_list for j in i]
    total_selection = np.concatenate([ligand, residues_list])
    subset_selected = trajectory.atom_slice(total_selection)
    total_selection = sorted(total_selection)
    h_bonds = hunt_h_bonds_of_model(subset_selected, cutoff)
    ligand_residue = trajectory.topology.atom(ligand[0]).residue
    h_bonds_with_ligand = []
    for hbond in h_bonds:
        h_b = [total_selection[x] for x in hbond]
        donor_id, hydrogen_id, acceptor_id = h_b
        residue_donor = trajectory.topology.atom(donor_id)
        residue_acceptor = trajectory.topology.atom(acceptor_id)
        if donor_id in ligand:
            ligand_atom = residue_donor.name
            residue_acceptor_resseq = residue_acceptor.residue.resSeq
            if (residue_acceptor.name == residues_dict[residue_acceptor_resseq]) and (
                residue_acceptor_resseq in residues_dict.keys()):
                residue_to_bond = residue_acceptor
                print(ligand_atom, residue_to_bond)
        if acceptor_id in ligand:
            ligand_atom = residue_acceptor.name
            residue_donor_resseq = residue_donor.residue.resSeq
            if (residue_donor.name == residues_dict[residue_donor_resseq]) and (
                residue_donor_resseq in residues_dict.keys()):
                residue_to_bond = residue_donor
                print(ligand_atom, residue_to_bond)
        try:
            h_bonds_with_ligand.append((ligand_atom, residue_to_bond))
        except:
            traceback.print_exc()

    return h_bonds_with_ligand


def get_bonds_from_row_in_csv(row, residues_dictionary, cutoff=0.25, lig_resname="LIG", file_column="file_from",
                              path_to_topology="topology.pdb"):
    trajectory_file = row[file_column]
    model = row["numberOfAcceptedPeleSteps"]
    topology_path_splited = trajectory_file.split("/")[0:-2]
    topology_path = os.path.join("/".join(topology_path_splited), path_to_topology)
    structure = extract_single_model_from_structure(trajectory_file, topology_path, model)
    result = select_specific_h_bond_with_ligand(structure, residues_dictionary, lig_resname, cutoff)
    return result


def load_trajectories_of_csvs(csv_file, residues_dictionary, output_folder, lig_resname="LIG", cutoff=0.25, file_column="file_from", separator=";", path_to_topology="topology.pdb"):
    dataframe = pd.read_csv(csv_file, sep=separator)
    pool = mp.Pool(processes=4)
    hbonds_results = []
    for index, row in dataframe.iterrows():
        hbonds_results.append(pool.apply_async(get_bonds_from_row_in_csv,
                                                 (row, residues_dictionary, cutoff, lig_resname, file_column,
                                                  path_to_topology)))
    for n, hbonds in enumerate(hbonds_results):
        hbonds = hbonds.get()
        if hbonds:
            for hb in hbonds:
                ligand_atom, residue = hb
                residue = str(residue)
                if residue not in dataframe.columns:
                    dataframe[residue] = ""
                else:
                    if not dataframe.loc[n, residue]:
                        atom = dataframe.loc[n, residue]
                        dataframe.loc[n, residue] = atom + "" + ligand_atom
                    else:
                        dataframe.loc[n, residue] = ligand_atom
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    dataframe.to_csv(os.path.join(output_folder, os.path.basename(csv_file)), sep=separator)
    return print("FINISHED SUCCESSFULLY")

#load_trajectories_of_csvs("/home/carlespl/project/Almirall/validation_compounds/"
#                          "testing_smalldataset.csv",
#                          {642: "NZ", 688: "N", 689: "N", 690: "N"},
#                                "/home/carlespl/project/Almirall/testing_hbonder",
#                                cutoff=0.27)

