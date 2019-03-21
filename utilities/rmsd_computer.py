import multiprocessing as mp
import os
import numpy as np
import mdtraj as md
import pandas as pd
from Helpers.dataframe_handler import add_column_to_dataframe
from Helpers.structure_handler import open_structure_from_dataframe
from Helpers.structure_handler import read_structure_xtc
from Helpers.structure_handler import select_ligand_mdtraj
from Helpers.structure_handler import select_ligand_and_residues


def calculate_rmsd_of_selection(selection, trajectory_target, trajectory_reference):
    rmsd = np.sqrt(3 * np.mean((trajectory_target.xyz[:, selection, :] - trajectory_reference.xyz[:, selection, :]) ** 2,
                               axis=(1, 2)))
    rmsd_amstrongs = [float(x) * 10 for x in rmsd]
    return rmsd_amstrongs


def compute_ligand_rmsd(trajectory_target, trajectory_reference, ligand_resname="LIG", heavy=True):
    ligand = select_ligand_mdtraj(trajectory_target, ligand_resname, heavy)
    rmsd_amstrongs = calculate_rmsd_of_selection(selection=ligand, trajectory_target=trajectory_target,
                                                 trajectory_reference=trajectory_reference)
    for result in rmsd_amstrongs:
        print(result)
    return rmsd_amstrongs


def compute_ligand_and_close_rmsd(trajectory_target, trajectory_reference, residues_ids, ligand_resname="LIG",
                                  heavy=True):
    selection = select_ligand_and_residues(trajectory_target, list_of_residues_ids=residues_ids,
                                           ligand_resname=ligand_resname, heavy=heavy)
    rmsd_amstrongs = calculate_rmsd_of_selection(selection=selection, trajectory_target=trajectory_target,
                                                 trajectory_reference=trajectory_reference)
    for result in rmsd_amstrongs:
        print(result)
    return rmsd_amstrongs


def rmsd_calculator(file_path, file_reference, ligand_resname="LIG", pdb=True, heavy=True, residues_ids=False):
    topology_path_splited = file_path.split("/")[0:-2]
    topology_path = os.path.join("/".join(topology_path_splited), "topology.pdb")
    trajectory_target = read_structure_xtc(file_path, topology_path, pdb=False)
    trajectory_reference = read_structure_xtc(file_path=topology_path, topology_file_path=None, pdb=pdb)
    if not residues_ids:
        rmsd = compute_ligand_rmsd(trajectory_target, trajectory_reference, ligand_resname, heavy)
    else:
        rmsd = compute_ligand_and_close_rmsd(trajectory_target, trajectory_reference, residues_ids, ligand_resname,
                                             heavy)
    return rmsd


def main(csv_file, reference, path_to_simulations, output_filepath, file_column="file_from", separator=";", read_xtc=True,
         path_to_clusters=None, ligand_resname="LIG", pdb=True, column_name="RMSD_recomputed", heavy=True,
         residues_ids=False):
    data = open_structure_from_dataframe(csv_file, path_to_simulations, path_to_clusters, file_column,
                                         separator, read_xtc)
    rmsd_list = []
    multiprocessing_list = []
    pool = mp.Pool(processes=4)
    for n, file in enumerate(data):
        multiprocessing_list.append(pool.apply_async(rmsd_calculator,
                                                     (file, reference, ligand_resname, pdb, heavy, residues_ids)))
    for process in multiprocessing_list:
        rmsd_list.extend(process.get())

    dataframe = pd.read_csv(csv_file, sep=separator)
    add_column_to_dataframe(rmsd_list, column_name, dataframe)
    dataframe.to_csv(output_filepath, sep=separator, index=False)
    return dataframe


main(csv_file="/home/carlespl/project/Almirall/crystals/with_waters/summary_results/8_PKD_TYK2_SD006_with_waters_adaptive_results_summary.csv",
     reference="/home/carlespl/project/Almirall/crystals/with_waters/simulation_results/8_PKD_TYK2_SD006_with_waters/8_PKD_TYK2_SD006_with_waters_processed.pdb",
     output_filepath="/home/carlespl/project/Almirall/crystals/with_waters/summary_results/8_PKD_TYK2_SD006_with_waters_adaptive_results_summary.csv",
     path_to_simulations="/home/carlespl/project/Almirall/crystals/with_waters/simulation_results",
     ligand_resname="INX",
     column_name="RMSD_close_area",
     residues_ids="595 596 597 598 601 603 640 642 662 647 672 673 687 688 689 690 691 692 693 694 697 738 739 741 758 759")
