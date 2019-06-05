import os
import argparse
import prody
import pdb_prepare


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
    required_named.add_argument("pdb_receptor",
                                help="""Absolute path to the PDB file used as receptor. Remember to use Protein
                                Preparation Wizard of maestro before running this program. Additionally, take in
                                mind that you must ensure that the PDB atom names are correct to then run PELE.""")
    required_named.add_argument("pdb_ligand",
                                help="""Absolute path the PDB file of the desired ligand. Take into account that
                                this software does not perform superimpositions, so ensure that the ligand is
                                correctly placed in the binding site (previous docking required).""")

    parser.add_argument("-c", "--chain_ligand", default="L",
                        help="""Chain name for the ligand.""")
    args = parser.parse_args()

    return args.pdb_receptor, args.pdb_ligand, args.chain_ligand


def compute_center_of_chain(pdb_object, chain="L"):
    molecule_to_center = pdb_object.select("chain {}".format(chain))
    center = prody.calcCenter(molecule_to_center)
    return center


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


def main(pdb_receptor, pdb_ligand, ligand_chain="L"):
    receptor, center = get_receptor(pdb_file=pdb_receptor, chain_to_center=ligand_chain)
    ligand = prody.parsePDB(pdb_ligand)
    complex_stru = add_ligand_to_receptor(receptor_object=receptor, ligand_object=ligand, ligand_chain=ligand_chain)
    complex_filename = os.path.splitext(pdb_receptor)[0] + "_" + os.path.basename(os.path.splitext(pdb_ligand)[0]) + \
                       ".pdb"
    prody.writePDB(complex_filename, complex_stru)
    pdb_corrected = pdb_prepare.add_ters_to_pdb(complex_filename)
    with open(complex_filename, "w") as pdb_out:
        pdb_out.write(pdb_corrected)
    return center


if __name__ == '__main__':
    pdb_receptor, pdb_ligand, ligand_chain = parse_arguments()
    main(pdb_receptor=pdb_receptor, pdb_ligand=pdb_ligand, ligand_chain=ligand_chain)

