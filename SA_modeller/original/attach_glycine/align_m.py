#!/usr/bin/env python

from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

# directories for input atom files
env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)
env.io.atom_files_directory = ['.', 'final_mirror.pdb']

a = automodel(env, alnfile = 'alignment_m.ali',
              knowns = 'final_mirror', sequence = 'output_m')

a.starting_model = 1
a.ending_model   = 5

# Very thorough VTFM optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 300

# Thorough MD optimization:
a.md_level = refine.slow

a.make()
