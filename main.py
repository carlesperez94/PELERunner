# General imports
import os
import shutil
import subprocess
import glob
import argparse
# Local imports
from Helpers import folder_handler, constraints, center_of_mass, templatize
import configuration as c

FilePath = os.path.abspath(__file__)
PackagePath = os.path.dirname(FilePath)
curr_dir = os.path.abspath(os.path.curdir)


def parse_arguments():
    """
        Parse user arguments
        Output: list with all the user arguments
    """
    parser = argparse.ArgumentParser(description="""Description: This program will prepare several folders of PELE ready
    to run given several pdb complexes and templatized control files.
    """)
    parser.add_argument("folder_to_analyze", type=str,
                        help="""Path to the folder that will be analyzed, it must contain all complexes in PDB 
                        format.""")
    parser.add_argument("control_template", type=str,
                        help="""Path to the control file templatized.""")
    parser.add_argument("--obc", type=bool, default=False,
                        help="If it is set, OBC parameters will be computed")

    args = parser.parse_args()

    return args.folder_to_analyze, args.control_template, args.obc


def control_file_modifier(control_template, pdb, license, overlap, results_path="/growing_output", steps=100,
                          chain="L", constraints=" ", center="", temperature=1000):
    """
    This function creates n control files for each intermediate template created in order to change
    the logPath, reportPath and trajectoryPath to have all control files prepared for PELE simulations.
    """

    ctrl_fold_name = "control_folder"
    controlfile_path = control_template
    control_template = os.path.basename(control_template)

    # As pdb is a list of complexes, we have to create one line per complex in pdb and then put them all in the control
    list_of_lines_complex = []
    for complex in pdb:
        control_file_complex = '{"files" : [{"path": "%s" }] }' % (complex)
        list_of_lines_complex.append(control_file_complex)
    lines_complex = ",\n".join(list_of_lines_complex)

    # Definition of the keywords that we are going to substitute from the template
    keywords = {"LICENSE": license,
                "RESULTS_PATH": results_path,
                "CHAIN": chain,
                "CONSTRAINTS": constraints,
                "CENTER": center,
                "PDB": lines_complex,
                "STEPS": steps,
                "OVERLAP": overlap,
                "TEMPERATURE": temperature,
                }
    # Creation of a folder where we are going to contain our control files, just if needed
    if not os.path.exists(ctrl_fold_name):
        os.mkdir(ctrl_fold_name)

    # Create a copy of the control template in the control folder, because templatize.py replace the original template
    if not os.path.exists(os.path.join(ctrl_fold_name, control_template)):
        shutil.copyfile(controlfile_path, control_template)

    # Else, if has been created this means that we already have a template in this folder, so we will need a copy of the
    # file in the main folder to then replace the template for a real control file
    else:
        shutil.copyfile(controlfile_path, control_template)

    # Modifying the control file template
    templatize.TemplateBuilder(control_template, keywords)
    # Make a copy in the control files folder
    print("{} has been created successfully!".format(control_template))


def simulation_runner(path_to_pele, control_in, cpus=4):
    """
    Runs a PELE simulation with the parameters described in the input control file.

    Input:

    path_to_pele --> Complete path to PELE folder

    control_in --> Name of the control file with the parameters to run PELE
    """
    if cpus:
        cpus = int(cpus)
        if cpus < 2:
            exit("Sorry, to run mpi PELE you need at least 2 CPUs!")
        else:
            print("Starting PELE simulation. You will run mpi PELE with {} cores.".format(cpus))
            cmd = "mpirun -np {} {} {}".format(cpus, path_to_pele, control_in)
            print("Running {}".format(cmd))
            subprocess.call(cmd.split())
    else:
        print("Starting PELE simulation. You will run serial PELE.")
        cmd = "{} {}".format(path_to_pele, control_in)
        print("Running {}".format(cmd))
        subprocess.call(cmd.split())


def extract_ligand_of_template(in_pdb_file, out_path, lig_chain="L"):
    with open(in_pdb_file) as pdb:
        pdb_as_list = pdb.readlines()
    ligand_lines = []
    for line in pdb_as_list:
        if line[21:22] == lig_chain:
            ligand_lines.append(line)
    try:
        ligand = "".join(ligand_lines)
    except:
        print("Ligand Chain empty. Check if the chain {} contains the ligand or not".format(lig_chain))
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    resname = ligand[17:20]
    out_prefix = os.path.join(out_path, resname)
    out_pdb = "{}.pdb".format(out_prefix)
    counter = 0
    while True:
        try:
            with open(out_pdb, "w") as out:
                out.write(ligand)
            break
        except FileExistsError:
            counter += 1
            out_pdb = out_prefix + "_" + str(counter) + ".pdb"
    return out_pdb

def run_plop_from_pdb(sch_python, plop_relative_path, pdb_file, py2_env):
    cmd = "{} {} {}".format(sch_python, plop_relative_path, pdb_file)
    new_env = os.environ.copy()
    new_env["PYTHONPATH"] = py2_env
    subprocess.call(cmd.split(), env=new_env)


def prepare_obc_parameters(sch_python, obc_param_path, template_file, folder):
    cmd = "{} {} {}".format(sch_python, obc_param_path, template_file)
    subprocess.call(cmd.split())
    obc_path_in_data = "Data/OBC/solventParamsHCTOBC.txt"
    with open(obc_path_in_data) as obc_data:
        obc = obc_data.read()
    obc_template = "{}_OBCParams.txt".format(template_file)
    with open(obc_template) as obc_tmpl:
        obc_lig = obc_tmpl.read()
    with open(os.path.join(folder,"DataLocal/OBC/solventParamsHCTOBC.txt"), "w") as final_obc:
        final_obc.write(obc + obc_lig)


def prepare_pele_simulation(pdb_complex, control_template, obc=False, plop_path=c.PLOP_PATH, out_ligands=c.PATH_OUTPUT_LIGANDS,
                           sch_python=c.SCHRODINGER_PY_PATH, py2_env=c.PYTHON2_SCH_ENV, results_folder=c.RESULTS_PATH,
                           license_path=c.LICENSE, overlap=c.OVERLAP, pele_steps=c.STEPS, chain=c.CHAIN,
                           temp=c.TEMPERATURE):
    # Path definition
    plop_relative_path = os.path.join(PackagePath, plop_path)
    # Creation of output folder
    folder_handler.check_and_create_DataLocal()
    # Creating symbolic links
    folder_handler.create_symlinks(c.PATH_TO_PELE_DATA, 'Data')
    folder_handler.create_symlinks(c.PATH_TO_PELE_DOCUMENTS, 'Documents')
    # Creating constraints
    const = "\n".join(constraints.retrieve_constraints(pdb_complex, {}, {}, 0.2, 0.2, 10))
    # Creating templates
    out_pdb = extract_ligand_of_template(in_pdb_file=pdb_complex, out_path=out_ligands, lig_chain=chain)
    template_name = out_pdb.split("/")[-1].split(".pdb")[0].lower()+"z"
    run_plop_from_pdb(sch_python=sch_python, plop_relative_path=plop_relative_path, pdb_file=out_pdb, py2_env=py2_env)
    if obc:
        prepare_obc_parameters(sch_python, obc_param_path=c.OBC_PATH,
                               template_file="DataLocal/Templates/OPLS2005/HeteroAtoms/{}".format(template_name), folder=".")
    # Control file preparation
    center = center_of_mass.center_of_mass(out_pdb)
    control_file_modifier(control_template=control_template, pdb=[pdb_complex], license=license_path, overlap=overlap,
                          results_path=results_folder, steps=pele_steps, chain=chain, constraints=const, center=center,
                          temperature=temp)
    # Creating results folder
    folder_handler.check_and_create_folder(results_folder)


def main(folder_to_analyze, control_template, obc=False, plop_path=c.PLOP_PATH, out_ligands=c.PATH_OUTPUT_LIGANDS,
        sch_python=c.SCHRODINGER_PY_PATH, py2_env=c.PYTHON2_SCH_ENV, results_folder=c.RESULTS_PATH,
        license_path=c.LICENSE, overlap=c.OVERLAP, pele_steps=c.STEPS, chain=c.CHAIN,
        temp=c.TEMPERATURE):
    pdbs_in_folder = glob.glob("{}/*.pdb".format(folder_to_analyze))
    for pdb in pdbs_in_folder:
        new_folder = os.path.join(folder_to_analyze, os.path.basename(pdb).split(".")[0])
        if not os.path.exists(new_folder):
             os.mkdir(new_folder)
        os.chdir(new_folder)
        prepare_pele_simulation(pdb, control_template, obc, plop_path=plop_path, out_ligands=out_ligands,
                                sch_python=sch_python, py2_env=py2_env, results_folder=results_folder,
                                license_path=license_path, overlap=overlap, pele_steps=pele_steps, chain=chain,
                                temp=temp)


if __name__ == '__main__':
    fol2analyze, cntr_temp, obc = parse_arguments()
    main(folder_to_analyze=fol2analyze, control_template=cntr_temp, obc=obc)

