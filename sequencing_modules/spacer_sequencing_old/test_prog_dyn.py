import os
import sys
import inspect
import prog_dyn

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0, parentdir+"/synthesis_simulation")

import utils.dna_file_reader as dfr

fragments = dfr.read_fasta("2_fragmented_seq_file.fasta")

OVLP_MIN = 20

for namei, fi in fragments.items():
    for namej, fj in fragments.items():
        w = prog_dyn.compareFrag(fi[-OVLP_MIN:],fj[:OVLP_MIN],10)
        print(fi[-OVLP_MIN:],fj[:OVLP_MIN], w)