import sys
import os
from string import Template


class TemplateBuilder(object):
    """
        Description: Base clase to substitute de keywords
        int the already templetize file.
        e.g.
        If we have a file like that --> dededea $RESIDUENAME $CHAIN
        to replace $RESIDUENAME to ASL and $CHAIN to Z we should do:
        file= /path/to/file
        keywords= { "RESIDUENAME" : "ASL", "CHAIN": "Z" }
        TemplateBuilder(file, keywords)
    """

    def __init__(self, file, keywords, fileout):

        self.file = file
        self.keywords = keywords
        self.fileout = fileout
        self.fill_in()

    def fill_in(self):
        """
        Fill the control file in
        """
        with open(os.path.join(self.file), 'r') as infile:
            confile_data = infile.read()

        confile_template = Template(confile_data)

        confile_text = confile_template.safe_substitute(self.keywords)

        with open(os.path.join(self.fileout), 'w') as outfile:
            outfile.write(confile_text)


def main(pele_control_file, pele_fileout, adaptive_control_file, adaptive_fileout, license, chain, center, radius,
         resnum, native, pdb_in, ligname, chain_dist_1=None, chain_dist_2=None, resnum_dist_1=None, resnum_dist_2=None,
         atom_dist_1=None, atom_dist_2=None):
    keywords = {"LICENSE": license,
                "CHAIN": chain,
                "CENTER": center,
                "RADIUS": radius,
                "RESNUM": resnum,
                "NATIVE": native,
                # Identify atoms to compute distances:
                "CHAIN_1": chain_dist_1,
                "CHAIN_2": chain_dist_2,
                "RESNUM_1": resnum_dist_1,
                "RESNUM_2": resnum_dist_2,
                "ATOM_1": atom_dist_1,
                "ATOM_2": atom_dist_2
                }

    keywords_adaptive = {"PDB_IN": pdb_in,
                         "CONTROL_FILE": pele_control_file,
                         "LIGAND_NAME": ligname,
                         }
    TemplateBuilder(pele_control_file, keywords, fileout=pele_fileout)
    TemplateBuilder(adaptive_control_file, keywords_adaptive, fileout=adaptive_fileout)
