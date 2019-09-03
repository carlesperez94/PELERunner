import os
import argparse
import subprocess
import glob
import traceback

import prody

import pdb_prepare
import prepare_controls
import prepare_folders
import configuration as c

FilePath = os.path.abspath(__file__)
package_path = os.path.dirname(FilePath)


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""This program takes a PDB of a receptor and a PDB of a ligand and
    merge both into the same PDB file, making all necessary changes to run PELE afterward (TER additions, water molecules
    renaming, etc.) """)
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-pdb_receptor",
                                help="""Absolute path to the PDB file used as receptor. Remember to use Protein
                                Preparation Wizard of maestro before running this program. Additionally, take in
                                mind that you must ensure that the PDB atom names are correct to then run PELE.""")
    required_named.add_argument("-pdb_ligand", type=str,
                                help="""Absolute path the PDB file of the desired ligand. Take into account that
                                this software does not perform superimpositions, so ensure that the ligand is
                                correctly placed in the binding site (previous docking required).""")
    required_named.add_argument("-adapt_template",
                                help="""""")
    required_named.add_argument("-pele_template",
                                help="""""")

    parser.add_argument("-c", "--chain_ligand", default="L",
                        help="""Chain name for the ligand.""")
    args = parser.parse_args()

    return args.pdb_receptor, args.pdb_ligand, args.adapt_template, args.pele_template, args.chain_ligand


def compute_center_of_chain(pdb_object, chain="L"):
    molecule_to_center = pdb_object.select("chain {}".format(chain))
    center = prody.calcCenter(molecule_to_center)
    return center


def get_resnum(pdb_object, chain="L"):
    molecule = pdb_object.select("chain {}".format(chain))
    resnum = molecule.getResnums()[0]
    return resnum


def get_receptor(pdb_file, chain_to_center="L"):
    pdb = prody.parsePDB(pdb_file)
    ligand_or_metal = pdb.select("hetero and not water")
    if chain_to_center:
        center = compute_center_of_chain(pdb, chain_to_center)
    if ligand_or_metal:
        answer = input("WARNING! Ligand or metal detected. This part of the pdb will be deleted. Do you want to \n"
                       "continue? [y/n]: ")
        while True:
            if answer == "y":
                receptor = pdb.select("protein or water")
                break
            elif answer == "n":
                answer2 = input("Select the ligand chain: ")
                ligand_chain = answer2
                receptor = pdb.select("(protein or water) or (hetero and not chain {})".format(ligand_chain))
                break
            else:
                print("Ligand or metal detected. This part of the pdb will be deleted. Do you want to \n"
                       "continue? [y/n]: ")

    return receptor.copy(), center


def add_ligand_to_receptor(receptor_object, ligand_object, ligand_chain="L"):
    chains_in_receptor = list(set(receptor_object.getChids()))
    while True:
        if ligand_chain in chains_in_receptor:
            ligand_chain = input("The receptor already contains the chain '{}'. "
                                 "Please, select another name: ".format(ligand_chain))
        else:
            break

    ligand_object.setChids(ligand_chain)
    complex_structure = receptor_object + ligand_object
    return complex_structure


def prepare_ligands_templates(sch_path, plop_relative_path, pdb_lig, rotamers, out_temp, out_rot):
    sch_python = os.path.join(sch_path, "run")
    cmd = "{} {} {} {} {} {}".format(sch_python, plop_relative_path, pdb_lig, rotamers, out_temp, out_rot)
    print(cmd)
    try:
        subprocess.call(cmd.split())
    except OSError:
        print("WARNING: You are running an SCHRODINGER version older than 2019!")
        try:
            sch_python = os.path.join(sch_path, "python")
            cmd = "{} {} {} {} {} {}".format(sch_python, plop_relative_path, pdb_lig, rotamers, out_temp, out_rot)
            subprocess.call(cmd.split())
        except OSError:
            raise OSError("Path {} not foud. Change schrodinger path under configuration.py".format(sch_python))


def prepare_obc_parameters(sch_python, obc_param_path, template_file, folder):
    if not sch_python.endswith("python"):
    	cmd = "{}/python {} {}".format(sch_python, obc_param_path, template_file)
    else:
        try:
            cmd = "{} {} {}".format(sch_python, obc_param_path, template_file)
        except FileNotFoundError:
            raise FileNotFoundError("OBC does not found python to run the script, try to found schrodinger python"
                          " and add the path properly!")
    subprocess.call(cmd.split())
    obc_path_in_data = "{}/Data/OBC/solventParamsHCTOBC.txt".format(folder)
    with open(obc_path_in_data) as obc_data:
        obc = obc_data.read()
    obc_template = "{}_OBCParams.txt".format(template_file)
    with open(obc_template) as obc_tmpl:
        obc_lig = obc_tmpl.read()
    with open(os.path.join(folder, "DataLocal/OBC/solventParamsHCTOBC.txt"), "w") as final_obc:
        final_obc.write(obc + obc_lig)


def prepare_single_pdb(pdb_receptor, pdb_ligand, adaptive_template, pele_template, documents_path=c.DOCUMENTS,
                       data_path=c.DATA, ligand_chain="L", license=c.LICENSE, rotamers=c.ROTAMERS, sch_path=c.SCHRODINGER,
                       obc_script_path=c.OBC_PATH, radius=4, ligname="LIG", chain_dist_1=None, chain_dist_2=None,
                       resnum_dist_1=None, resnum_dist_2=None, atom_dist_1=None, atom_dist_2=None):

    receptor, center = get_receptor(pdb_file=pdb_receptor, chain_to_center=ligand_chain)
    ligand = prody.parsePDB(pdb_ligand)
    complex_stru = add_ligand_to_receptor(receptor_object=receptor, ligand_object=ligand, ligand_chain=ligand_chain)
    resnum = get_resnum(ligand, ligand_chain)
    out_folder = os.path.splitext(pdb_receptor)[0] + "_" + os.path.basename(os.path.splitext(pdb_ligand)[0])
    complex_filename = os.path.splitext(pdb_receptor)[0] + "_" + os.path.basename(os.path.splitext(pdb_ligand)[0]) + \
                       ".pdb"
    os.mkdir(out_folder)
    prody.writePDB(os.path.join(out_folder, complex_filename), complex_stru)
    pdb_corrected = pdb_prepare.add_ters_to_pdb(os.path.join(out_folder, complex_filename))
    with open(os.path.join(out_folder, complex_filename), "w") as pdb_out:
        pdb_out.write(pdb_corrected)
    center_string = ", ".join([str(center[0]), str(center[1]), str(center[2])])
    prepare_controls.main(pele_control_file=pele_template, pele_fileout=os.path.join(out_folder, pele_template),
                          adaptive_control_file=adaptive_template, adaptive_fileout=os.path.join(out_folder,
                                                                                                 adaptive_template),
                          license=license, chain=ligand_chain,
                          center=center_string, radius=radius, resnum=resnum, native=complex_filename,
                          pdb_in=complex_filename, ligname=ligname, chain_dist_1=chain_dist_1, chain_dist_2=chain_dist_2,
                          resnum_dist_1=resnum_dist_1, resnum_dist_2=resnum_dist_2, atom_dist_1=atom_dist_1,
                          atom_dist_2=atom_dist_2)
    prepare_folders.main(destination_path=out_folder, documents_path=documents_path, data_path=data_path)
    prody.writePDB(os.path.join(out_folder, "{}.pdb".format(ligname)), ligand)
    out_temp = os.path.join(out_folder, "DataLocal/Templates/OPLS2005/HeteroAtoms")
    out_rot = os.path.join(out_folder, "DataLocal/LigandRotamerLibs")
    plop_relative_path = os.path.join(package_path, "PlopRotTemp_S_2017/ligand_prep.py")
    prepare_ligands_templates(sch_path, plop_relative_path, pdb_ligand, rotamers, out_temp, out_rot)
    if obc_script_path:
        prepare_obc_parameters(sch_path, obc_script_path, os.path.join(out_temp, ligname.lower() + "z"), out_folder)


def main(pdb_receptor, ligand_folder, adaptive_template, pele_template, documents_path=c.DOCUMENTS,
         data_path=c.DATA, ligand_chain="L", license=c.LICENSE, rotamers=c.ROTAMERS, sch_path=c.SCHRODINGER,
         obc_script_path=c.OBC_PATH, radius=4, ligname="LIG", distances_instructions=c.DISTANCE_ATOMS):
    ligands = glob.glob(ligand_folder)
    print("LIGANDS found: {}".format(ligands))
    ligands = sorted(ligands)
    if distances_instructions:
        for pdb_ligand, distance_ins in zip(ligands, distances_instructions):
            reference_1, reference_2 = distance_ins
            chain_1, residue_1, atom_1 = reference_1
            chain_2, residue_2, atom_2 = reference_2
            try:
                print("Preparing {}... ".format(pdb_ligand))
                prepare_single_pdb(pdb_receptor, pdb_ligand, adaptive_template, pele_template, documents_path, data_path,
                                   ligand_chain, license, rotamers, sch_path, obc_script_path, radius, ligname, chain_1,
                                   chain_2, residue_1, residue_2, atom_1, atom_2)
            except Exception:
                traceback.print_exc()
    else:
        for pdb_ligand in ligands:
            try:
                print("Preparing {}... ".format(pdb_ligand))
                prepare_single_pdb(pdb_receptor, pdb_ligand, adaptive_template, pele_template, documents_path,
                                   data_path, ligand_chain, license, rotamers, sch_path, radius, ligname)
            except Exception:
                traceback.print_exc()


if __name__ == '__main__':
    pdb_receptor, ligand_folder, adapt_template, pele_template, ligand_chain = parse_arguments()
    main(pdb_receptor=pdb_receptor, ligand_folder=ligand_folder, ligand_chain=ligand_chain,
         pele_template=pele_template, adaptive_template=adapt_template)

