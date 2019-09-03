import glob
import os


def open_pdb_file_as_lines(path_to_pdb):
    with open(path_to_pdb) as pdb:
        content = pdb.readlines()
    return content


def find_ligand_in_complex(complex, lig_chain="L"):
    pdb = open_pdb_file_as_lines(complex)
    ligand_lines = []
    for line in pdb:
        if line[21:22] == lig_chain:
            ligand_lines.append(line)
    ligand = "".join(ligand_lines)
    return ligand


def replace_ligand1_by_ligand2(ligand1, ligand2, chain_lig1, chain_lig2):
    ligands = []
    for ligand, chain in zip([ligand1, ligand2], [chain_lig1, chain_lig2]):
        ligand_in_pdb = find_ligand_in_complex(complex=ligand, lig_chain=chain)
        print("LIGAND:\n{}\nCHAIN:{}\n".format(ligand_in_pdb, chain))
        ligands.append(ligand_in_pdb)
    with open(ligand1) as pdb:
        pdb_to_replace = pdb.read()
    old_ligand, new_ligand = ligands[0], ligands[1]
    pdb_replaced = pdb_to_replace.replace(old_ligand, new_ligand)
    return pdb_replaced


def main(complex, docked_pdbs_folder, chain_ligand_in_complex="L", chain_ligand_in_docked_files="L"):
    docked_files = glob.glob("{}/*.pdb".format(docked_pdbs_folder))
    print(docked_files)
    file_names = []
    for docking in docked_files:
        new_pdb = replace_ligand1_by_ligand2(ligand1=complex, ligand2=docking, chain_lig1=chain_ligand_in_complex,
                                             chain_lig2=chain_ligand_in_docked_files)
        name = "{}_complex.pdb".format(docking.split(".pdb")[0])
        print(name)
        with open(name, "w") as out_pdb:
            out_pdb.write(new_pdb)
        file_names.append(name)
    return file_names
