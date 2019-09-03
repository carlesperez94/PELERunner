import sys
import os
import socket
# PATHS Definitions

machine = socket.getfqdn()
if "bsc.mn" in machine:
    # PELE parameters
    PATH_TO_PELE = "/gpfs/projects/bsc72/PELE++/mniv/rev12536/bin/Pele_mpi"
    PATH_TO_PELE_DATA = "/gpfs/projects/bsc72/PELE++/data/rev12360/Data"
    PATH_TO_PELE_DOCUMENTS = "/gpfs/projects/bsc72/PELE++/Documents/rev12360"
    LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    # PlopRotTemp parameters
    SCHRODINGER_PY_PATH = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC/utilities/python"
    PYTHON2_SCH_ENV = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC/internal/lib/python2.7/site-packages/:/home/bsc72/bsc72292/projects/PELERunner/PELERunner"
    OBC_PATH = "/gpfs/projects/bsc72/PELE++/scripts/solventOBCParamsGenerator.py"
else:  # If you run the program in another machine fill this parameters
    # PELE parameters
    PATH_TO_PELE = ""
    PATH_TO_PELE_DATA = ""
    PATH_TO_PELE_DOCUMENTS = ""
    LICENSE = ""
    # PlopRotTemp parameters
    SCHRODINGER_PY_PATH = ""
    OBC_PATH = ""


# Control file parameters to replace
STEPS = 1000
RESULTS_PATH = "simulation_results"
CHAIN = "L"
TEMPERATURE = 1500
OVERLAP = 0.7
PROCESSORS = 48
# PlopRotTemp configuration
PLOP_PATH = "PlopRotTemp_S_2017/ligand_prep.py"
ROTRES = 10

PATH_OUTPUT_LIGANDS = "."

