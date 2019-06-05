import glob
import multiprocessing as mp
import os
import traceback

import mdtraj as md
import prody
import numpy as np
import pandas as pd
from xtc_to_pdb import trajectory_and_snapshot_to_pdb


def select_ligand(pdb, ligand_chain="L"):
    ligand = pdb.select("chain {}".format(ligand_chain))
    return ligand


def select_ligand_md(model, lig_resname="LIG"):
    topology = model.topology
    ligand = topology.select("resname '{}'".format(lig_resname))
    #ligand_struct = model.atom_slice(ligand)
    return ligand


def select_prot_and_waters(model):
    topology = model.topology
    protein = topology.select("protein or water")
    protein_struct = model.atom_slice(protein)
    return protein_struct


def select_possible_donors_and_acceptors_protein(pdb, cutoff=3.5, ligand_chain="L"):
    acceptors_and_donors = pdb.select("protein and (element N C O H within {} of chain {})".format(cutoff, ligand_chain))
    return acceptors_and_donors


def find_donors(selection, cutoff=1.2):
    pdb_atom_names = selection.getNames()
    residues = selection.getResnums()
    elements = selection.getElements()
    donors = []
    for atom, id, element in zip(pdb_atom_names, residues, elements):
        if element == "H":
            donor_group = selection.select("(element N C O within {} of (name {} and resnum {})) or (element H and name"
                                           " {} and resnum {})".format(cutoff, atom, id, atom, id))
            if len(donor_group) > 1:
                donors.append(donor_group)
                #prody.writePDB("/home/carlespl/project/Almirall/donor_{}_{}".format(atom, id), donor_group)
    return donors


def find_donors_md(selection):
    topology = selection.topology
    hydrogen_atoms = topology.select("element H")
    return hydrogen_atoms


def find_acceptors(selection, cutoff=1.2):
    pdb_atom_names = selection.getNames()
    residues = selection.getResnums()
    elements = selection.getElements()
    acceptors = []
    for atom, id, element in zip(pdb_atom_names, residues, elements):
        if element == "N":
            acceptor_group = selection.select("(name {} and resnum {}) or (element H within {} of (name {} and "
                                               "resnum {}))".format(atom, id, cutoff, atom, id))
            if len(acceptor_group) < 2:
                acceptors.append(acceptor_group)
        if element == "O":
            acceptor_group = selection.select("(name {} and resnum {}) or (element H within {} of (name {} and "
                                              "resnum {}))".format(atom, id, cutoff, atom, id))
            if len(acceptor_group) <= 2:
                acceptors.append(acceptor_group)
            #prody.writePDB("/home/carlespl/project/Almirall/acceptor_{}_{}".format(atom, id), acceptor_group)
    return acceptors


def find_acceptors_md(selection):
    topology = selection.topology
    possible_acceptors_atoms = topology.select("element N or element O")
    acceptors_banned = []
    for atom in possible_acceptors_atoms:
        atoms_bonded_to_possible_acceptors = md.compute_neighbors(selection, 0.165, [atom])[0]
        if len(atoms_bonded_to_possible_acceptors) > 0:  # If bond check if has a H atom bonded to ban them
            checker = False
            for atom_bonded in atoms_bonded_to_possible_acceptors:
                if str(topology.atom(atom_bonded).element) == "hydrogen" and str(topology.atom(atom).element) != "oxygen":
                    checker = True
            if checker:
                acceptors_banned.append(atom)
    return list(set(possible_acceptors_atoms) - set(acceptors_banned))


def merge_selections(list_of_selections):
    i = 0
    while i < len(list_of_selections):
        if i == 0:
            selections = list_of_selections[i]
        else:
            selections = selections + list_of_selections[i]
        i += 1
    return selections


def sort_triplet(triplet):
    sorted_triplet = []
    if len(triplet) == 3:
        for n, atom in enumerate(triplet):
            if atom.getElement() != "H":
                sorted_triplet.append(atom)
            else:
                sorted_triplet.insert(1, atom)
        return sorted_triplet


def check_h_bond_cond_accomplished(acceptors, donors, hbond_cutoff=2.5, angle_cutoff=130):
    hbonds = []
    for acceptor in acceptors:
        acceptor_name = acceptor.getNames()[0]
        acceptor_id = acceptor.getResnums()[0]
        acceptor_resname = acceptor.getResnames()[0]
        if len(acceptor) > 1:
            acceptor = acceptor[0]
            acceptor_name = acceptor.getName()
            acceptor_id = acceptor.getResnum()
            acceptor_resname = acceptor.getResname()
        for donor in donors:
            triplet = acceptor + donor
            h_at_distance = triplet.select("hydrogen within {} of name {}".format(hbond_cutoff, acceptor_name))
            try:
                name_h = h_at_distance.getNames()[0]
                id_h = h_at_distance.getResnums()[0]
                id_resname = h_at_distance.getResnames()[0]
                triplet_sorted = sort_triplet(triplet)
                angle = prody.calcAngle(triplet_sorted[0], triplet_sorted[1], triplet_sorted[2])
                if abs(angle) > angle_cutoff:
                    print("H-bond detected between {}{}-{} and {}{}-{}". format(acceptor_id, acceptor_resname,
                                                                                acceptor_name, id_h, id_resname,
                                                                                name_h))
                    hbonds.append(triplet_sorted)
            except AttributeError:
                pass
    return hbonds


def get_hbonds_from_pdb(pdb_file, ligand_chain="L", residues_cutoff=3.5, hbond_cutoff=2.5, angle_cutoff=130):
    hbonds = []
    pdb = prody.parsePDB(pdb_file)
    ligand = select_ligand(pdb=pdb, ligand_chain=ligand_chain)
    residues = select_possible_donors_and_acceptors_protein(pdb=pdb, cutoff=residues_cutoff, ligand_chain=ligand_chain)
    donors_protein = find_donors(residues)
    donors_ligand = find_donors(ligand)
    acceptors_protein = find_acceptors(residues)
    acceptors_ligand = find_acceptors(ligand)
    hbond_acceptor_lig = check_h_bond_cond_accomplished(acceptors=acceptors_ligand, donors=donors_protein,
                                                        hbond_cutoff=hbond_cutoff, angle_cutoff=angle_cutoff)
    hbond_acceptor_prot = check_h_bond_cond_accomplished(acceptors=acceptors_protein, donors=donors_ligand,
                                                         hbond_cutoff=hbond_cutoff, angle_cutoff=angle_cutoff)
    for triplet_lig in hbond_acceptor_lig:
        hbonds.append(triplet_lig)
    for triplet_prot in hbond_acceptor_prot:
        hbonds.append(triplet_prot)
    return hbonds


def get_hbonds_from_xtc(trajectory_path, ligand_resname="LIG", hbond_cutoff=0.27, angle_cutoff=130):
    print("TRAJECTORY: {}".format(trajectory_path))
    # Create an empty dataframe
    results = []
    topology_path_splited = trajectory_path.split("/")[0:-2]
    topology_path = os.path.join("/".join(topology_path_splited), "topology.pdb")
    trajectory = md.load(trajectory_path, top=topology_path)
    ligand = select_ligand_md(trajectory, lig_resname=ligand_resname)
    ligand_traj = trajectory.atom_slice(ligand)
    h_donors_ligand = find_donors_md(ligand_traj)
    h_donors_ligand = [ligand[index] for index in h_donors_ligand]
    donors_ligand = []
    for hydrogen in h_donors_ligand:
        donor = md.compute_neighbors(trajectory, cutoff=0.145, query_indices=[hydrogen])
        donor = donor[0][0]
        donors_ligand.append(donor)
    acceptors_ligand = find_acceptors_md(ligand_traj)
    acceptors_ligand = [ligand[index] for index in acceptors_ligand]
    binding_sites = md.compute_neighbors(trajectory, cutoff=hbond_cutoff, query_indices=ligand)
    angle_cut_in_rad = angle_cutoff * (np.pi / 180)
    for n, bs in enumerate(binding_sites):
        hbonds = []
        bs = [atom for atom in bs if not atom in ligand]
        binding_sites_traj = trajectory.atom_slice(bs)
        acceptors_protein = find_acceptors_md(binding_sites_traj)
        acceptors_protein = [bs[index] for index in acceptors_protein]
        h_donors_protein = find_donors_md(binding_sites_traj)
        h_donors_protein = [bs[index] for index in h_donors_protein]
        donors_protein = []
        for hydrogen in h_donors_protein:
            donor = md.compute_neighbors(trajectory, cutoff=0.145, query_indices=[hydrogen])
            if len(donor[0]) == 0:
                donor = md.compute_neighbors(trajectory, cutoff=0.16, query_indices=[hydrogen])
            donor = donor[0][0]
            donors_protein.append(donor)
        print("Model {}".format(n))
        for acceptor in acceptors_protein:
            for h, donor in zip(h_donors_ligand, donors_ligand):
                distance = md.compute_distances(trajectory, [[acceptor, h]])
                distance = distance[n][0]  # Get only the distance for the correct model
                if distance <= hbond_cutoff:
                    angle = md.compute_angles(trajectory, [[acceptor, h, donor]])
                    angle = angle[n][0]
                    if angle_cut_in_rad < angle:
                        print("Hbond between {} and {}".format(trajectory.topology.atom(donor),
                                                               trajectory.topology.atom(acceptor)))
                        hbonds.append([trajectory.topology.atom(donor), trajectory.topology.atom(acceptor)])
        for acceptor in acceptors_ligand:
            for h, donor in zip(h_donors_protein, donors_protein):
                distance = md.compute_distances(trajectory, [[acceptor, h]])
                distance = distance[n][0]  # Get only the distance for the correct model
                if distance <= hbond_cutoff:

                    angle = md.compute_angles(trajectory, [[acceptor, h, donor]])
                    angle = angle[n][0]
                    if angle_cut_in_rad < angle:
                        print("Hbond between {} and {}".format(trajectory.topology.atom(donor),
                                                               trajectory.topology.atom(acceptor)))
                        hbonds.append([trajectory.topology.atom(donor), trajectory.topology.atom(acceptor)])
        hbond_info = ["{}_{}".format(donor, acceptor) for donor, acceptor in hbonds]
        data = pd.DataFrame({"trajectory": trajectory_path, "model": n, "hbonds": hbond_info})
        results.append(data)

    return results


def process_snapshot_from_trajectory_file(trajectory, snapshot, lig_chain="L", residues_cutoff=3.5, hbond_cutoff=2.5,
                                          angle_cutoff=130):
    print("TRAJECTORY: {}".format(trajectory))
    results = []
    if trajectory.endswith(".xtc"):
        output_path = "{}_snapshot_{}.pdb".format(trajectory.split(".xtc")[0], snapshot)
    else:
        output_path = "{}_snapshot_{}.pdb".format(trajectory.split(".pdb")[0], snapshot)
    trajectory_and_snapshot_to_pdb(trajectory, snapshot, output_path)
    hbonds = get_hbonds_from_pdb(pdb_file=output_path, ligand_chain=lig_chain, residues_cutoff=residues_cutoff,
                                hbond_cutoff=hbond_cutoff, angle_cutoff=angle_cutoff)
    for hbond in hbonds:
        results.append(["{}{}-{}".format(hbond[0].getResnum(), hbond[0].getResname(), hbond[0].getName()),
                        "{}{}-{}".format(hbond[1].getResnum(), hbond[1].getResname(), hbond[1].getName()),
                        "{}{}-{}".format(hbond[2].getResnum(), hbond[2].getResname(), hbond[2].getName())])
    return trajectory, snapshot, results


def get_hbonds_from_trajectory_file(trajectory, lig_chain="L", residues_cutoff=3.5, hbond_cutoff=2.5, angle_cutoff=130,
                                    processors=4):
    topology_path_splited = trajectory.split("/")[0:-2]
    topology_path = os.path.join("/".join(topology_path_splited), "topology.pdb")
    trajectory_len = len(md.load(trajectory, top=topology_path))
    pool = mp.Pool(processors)
    hbonds = []
    multi = []
    for snapshot in range(0, trajectory_len):
        multi.append(pool.apply_async(process_snapshot_from_trajectory_file, [trajectory, snapshot, lig_chain,
                                                                              residues_cutoff, hbond_cutoff,
                                                                              angle_cutoff]))
    for process in multi:
        hbond = process.get()
        hbonds.append(hbond)
    pool.terminate()
    return hbonds


def main(simulation_results_path, traj_prefix="*trajectory_*.xtc", ligand_resname="LIG",
         hbond_cutoff=0.27, angle_cutoff=130, processors=4, output_path=None):
    paths_to_analyze = glob.glob("{}/[0-9]*/{}".format(simulation_results_path, traj_prefix))
    results = []
    pool = mp.Pool(processors)
    multi = []
    for trajectory in paths_to_analyze:
        multi.append(pool.apply_async(get_hbonds_from_xtc, [trajectory, ligand_resname, hbond_cutoff, angle_cutoff]))
    for process in multi:
        for df in process.get():
            results.append(df)
    pool.terminate()
    final_df = pd.concat(results)
    if not output_path:
        final_df.to_csv(os.path.join(simulation_results_path, "hbonds_analysis.csv"), index=False)
    else:
        final_df.to_csv(output_path, index=False)

folder_list = glob.glob("/home/carlespl/project/Almirall/validation_compounds/simulation_results_with_exp_data/*/obc_adaptive*")
for folder in folder_list:
    print(folder)
    try:
        main(folder, ligand_resname="LIG")
    except Exception:
        traceback.print_exc()
