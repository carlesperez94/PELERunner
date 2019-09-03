import os


def check_and_create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)


def check_and_create_datalocal(destination_path):
    check_and_create_folder((os.path.join(destination_path, "DataLocal")))
    check_and_create_folder((os.path.join(destination_path, "DataLocal/LigandRotamerLibs")))
    check_and_create_folder((os.path.join(destination_path, "DataLocal/OBC")))
    check_and_create_folder((os.path.join(destination_path, "DataLocal/Templates")))
    check_and_create_folder((os.path.join(destination_path, "DataLocal/Templates/OPLS2005")))
    check_and_create_folder((os.path.join(destination_path, "DataLocal/Templates/OPLS2005/HeteroAtoms")))


def create_symbolic_link(path_to, name):
    if not os.path.islink(path_to):
        os.symlink(path_to, name)


def main(destination_path, data_path, documents_path):
    check_and_create_datalocal(destination_path)
    create_symbolic_link(data_path, os.path.join(destination_path, "Data"))
    create_symbolic_link(documents_path, os.path.join(destination_path, "Documents"))
