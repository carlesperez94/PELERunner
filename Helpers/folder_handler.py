import sys
import os

def create_symlinks(src, dst):
    if not os.path.islink(dst):
        os.symlink(src, dst)


def check_and_create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)


def check_and_create_DataLocal():
    # Get current path
    curr_dir = os.path.abspath(os.path.curdir)
    check_and_create_folder((os.path.join(curr_dir, "DataLocal")))
    check_and_create_folder((os.path.join(curr_dir, "DataLocal/LigandRotamerLibs")))
    check_and_create_folder((os.path.join(curr_dir, "DataLocal/Templates")))
    check_and_create_folder((os.path.join(curr_dir, "DataLocal/Templates/OPLS2005")))
    check_and_create_folder((os.path.join(curr_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")))


