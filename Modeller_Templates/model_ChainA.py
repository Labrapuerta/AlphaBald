# MODELLER script for ChainA
# Run this script from the Modeller_Templates directory
from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()

# Template PDB should be in current directory
env.io.atom_files_directory = ['.']

a = AutoModel(env,
              alnfile='ChainA_alignment.pir',
              knowns='1nkp',
              sequence='ChainA',
              assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 5  # Generate 5 models

a.make()
