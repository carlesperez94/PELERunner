import socket

machine = socket.getfqdn()
print(machine)

if "bsc.mn" in machine:
    LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    DATA = "/gpfs/projects/bsc72/PELE++/mniv/V1.6/Data"
    DOCUMENTS = "/gpfs/projects/bsc72/PELE++/mniv/V1.6/Documents"
    SCHRODINGER = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC/utilities"
    OBC_PATH = "/gpfs/projects/bsc72/PELE++/scripts/solventOBCParamsGenerator.py"
else:
    LICENSE = "/home/carlespl/repos/PELE-repo/licenses"
    DATA = "/home/carlespl/repos/PELE-repo/Data"
    DOCUMENTS = "/home/carlespl/repos/PELE-repo/Documents"
    SCHRODINGER = "/home/carlespl/schrodinger2019-2"
    OBC_PATH = ""
ROTAMERS = "10.0"
RADIUS = 6
RESNUM = 1
DISTANCE_ATOMS = None
#DISTANCE_ATOMS = [
#                 [["A", 690, "_H__"] ,["L", 900 ,"_N5_"]],
#                 [["A", 690, "_H__"] ,["L", 900 ,"_N2_"]],
#                 [["A", 690, "_H__"] ,["L", 900 ,"_N2_"]],
#                 [["A", 690, "_H__"] ,["L", 900 ,"_N1_"]],
#                 [["A", 690, "_H__"] ,["L", 900 ,"_N4_"]],
#                  ]
