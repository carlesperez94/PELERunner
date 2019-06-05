import os
import sys
import argparse
import traceback
import glob
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
    residues = residues_dictionary.items()
    results = []
    for residue, name in residues:
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
    residues_list = [j for i in residues_list for j in i]  # Delete arrays and join them
    total_selection = np.concatenate([residues_list, ligand])
    total_selection.sort()  # Sorting is important, otherwise the structure is aberrant
    subset_selected = trajectory.atom_slice(total_selection)
    subset_selected.save_pdb("/home/carlespl/project/Almirall/selection.pdb")
    h_bonds = hunt_h_bonds_of_model(subset_selected, cutoff)
    total_selection_reindexed = {n: i for n, i in enumerate(total_selection)}
    h_bonds_with_ligand = []
    for hbond in h_bonds:
        donor = hbond[0]
        acceptor = hbond[2]
        donor_name = trajectory.topology.atom(total_selection_reindexed[donor]).name
        acceptor_name = trajectory.topology.atom(total_selection_reindexed[acceptor]).name
        donor_seq = trajectory.topology.atom(total_selection_reindexed[donor]).residue.resSeq
        acceptor_seq = trajectory.topology.atom(total_selection_reindexed[acceptor]).residue.resSeq
        if total_selection_reindexed[donor] in ligand or total_selection_reindexed[acceptor] in ligand:
            # Residue = acceptor, ligand = donor
            if acceptor_seq in residues_dict.keys():
                if acceptor_name in residues_dict[acceptor_seq]:
                    hbond_to_save = [(trajectory.topology.atom(total_selection_reindexed[acceptor]).residue.resSeq,
                                     trajectory.topology.atom(total_selection_reindexed[acceptor]).residue.name,
                                     acceptor_name),
                                     (trajectory.topology.atom(total_selection_reindexed[donor]).residue.resSeq,
                                     trajectory.topology.atom(total_selection_reindexed[donor]).residue.name,
                                     donor_name)]

                    print("The acceptor is the residue {}{} {}".format(trajectory.topology.atom(
                          total_selection_reindexed[acceptor]).residue.resSeq,
                                                                       trajectory.topology.atom(
                                                                           total_selection_reindexed[
                                                                               acceptor]).residue.name,
                                                                       trajectory.topology.atom(
                                                                           total_selection_reindexed[acceptor]).element))
                    print("The donor is the ligand {}{} {}{}".format(trajectory.topology.atom(
                        total_selection_reindexed[donor]).residue.resSeq,
                                                                     trajectory.topology.atom(
                                                                         total_selection_reindexed[donor]).residue.name,
                                                                     trajectory.topology.atom(
                                                                         total_selection_reindexed[donor]).element,
                                                                     trajectory.topology.atom(
                                                                         total_selection_reindexed[donor]).name))
                    h_bonds_with_ligand.append(hbond_to_save)
            if donor_seq in residues_dict.keys():
                if donor_name in residues_dict[donor_seq]:
                    hbond_to_save = [(trajectory.topology.atom(total_selection_reindexed[donor]).residue.resSeq,
                                      trajectory.topology.atom(total_selection_reindexed[donor]).residue.name,
                                      donor_name),
                                      (trajectory.topology.atom(total_selection_reindexed[acceptor]).residue.resSeq,
                                      trajectory.topology.atom(total_selection_reindexed[acceptor]).residue.name,
                                      acceptor_name)]
                    print("The donor is the residue atom {}{} {}".format(trajectory.topology.atom(
                        total_selection_reindexed[donor]).residue.resSeq,
                                                                         trajectory.topology.atom(
                                                                             total_selection_reindexed[donor]).residue.name,
                                                                         trajectory.topology.atom(
                                                                             total_selection_reindexed[donor]).element))
                    print("The acceptor is the ligand {}{} {}{}".format(trajectory.topology.atom(
                        total_selection_reindexed[acceptor]).residue.resSeq,
                                                                        trajectory.topology.atom(
                                                                            total_selection_reindexed[
                                                                                acceptor]).residue.name,
                                                                        trajectory.topology.atom(
                                                                            total_selection_reindexed[acceptor]).element,
                                                                        trajectory.topology.atom(
                                                                            total_selection_reindexed[acceptor]).name))
                    h_bonds_with_ligand.append(hbond_to_save)

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


def load_trajectories_of_csvs(csv_file, residues_dictionary, output_folder, lig_resname="LIG", cutoff=0.25,
                              file_column="file_from", separator=";", path_to_topology="topology.pdb"):
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


def check_bonds_of_folder(folder_to_analyze, residues_dictionary, lig_resname="LIG", cutoff=0.27, processors=4):
    pool = mp.Pool(processes=processors)
    hbonds = []
    trajectories = glob.glob("{}/*.xtc".format(folder_to_analyze))
    for trajectory in trajectories:
        topology_path_splited = trajectory.split("/")[0:-2]
        topology_path = os.path.join("/".join(topology_path_splited), "topology.pdb")
        traj = md.load(trajectory, top=topology_path)
        multi = []
        for n, model in enumerate(traj):
            multi.append(pool.apply_async(select_specific_h_bond_with_ligand,((model, residues_dictionary, lig_resname, cutoff))))
        for n, process in enumerate(multi):
            hbonds_lig = process.get()
            hbonds.append([trajectory, n, hbonds_lig])
            print(trajectory, n, len(hbonds_lig))
    return hbonds

trajectory_path = "/home/carlespl/project/Almirall/validation_compounds/simulation_results_with_exp_data/73_wo14670_13/obc_adaptive_output_conf_11/0/73_wo14670_13_trajectory_1.xtc"
topology_path = "/home/carlespl/project/Almirall/validation_compounds/simulation_results_with_exp_data/73_wo14670_13/obc_adaptive_output_conf_11/topology.pdb"
traj = md.load(trajectory_path, top=topology_path)
for model in traj[0]:
    print(model)
    select_specific_h_bond_with_ligand(model, residues_dict={688: ["O"], 690: ["N", "O"]}, cutoff=0.27)

#dataframe = pd.DataFrame(columns=["trajectory_file", "model", "hbonds"])
#paths = glob.glob("/home/carlespl/project/Almirall/validation_compounds/simulation_results_with_exp_data/73_wo14670_13/obc_adaptive_output_conf_11/[0-9]*")
#for path in paths:
#    results = check_bonds_of_folder(path, residues_dictionary={688: ["O"], 690: ["N", "O"]})
#    for hbond in results:
#        trajectory, model, hbonds = hbond
#        amount_hbonds = len(hbonds)
#        df = pd.DataFrame([[trajectory, model, amount_hbonds]], columns=["trajectory_file", "model", "hbonds"])
#        dataframe = pd.concat([dataframe, df])
#        print(dataframe)
#dataframe.to_csv("/home/carlespl/project/Almirall/csv_testing.csv")


#load_trajectories_of_csvs("/home/carlespl/project/Almirall/validation_compounds/"
#                          "testing_smalldataset.csv",
#                          {642: "NZ", 688: "N", 689: "N", 690: "N"},
#                                "/home/carlespl/project/Almirall/testing_hbonder",
#                                cutoff=0.27)

