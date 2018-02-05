#!/usr/bin/env python
"""
    This script takes a list of sequences and finds primers that go off of
    either end. There are eight values you must enter for it to work:
    1. Lowest acceptable Tm
    2. Highest acceptable Tm
    3. How far in from the end of the sequence the primer can be
    4. Minimum number of mismatches between the primer and the assembly
    5. Maximum times the primer can have a sequence in the assembly with fewer
       mismatches than #4
    6. Fasta file of sequences that you want primers for
    7. Blast database of the assembly
    8. Output file
    It tries all possible primers of from 18-27 nucleotides long.
"""

import sys
import subprocess
import re
from Bio.Seq import Seq

if len(sys.argv) != 9:
    print "Usage: create_pcr_primers.py <lowest Tm> <highest Tm> \
    <range from ends> <minimum mismatches> <max occurrences in assembly \
    (Min = 1)> <fasta file> <blast db> <output file>\n"
    quit()

low_tm = int(sys.argv[1])
hi_tm = int(sys.argv[2])
seq_range = int(sys.argv[3])
min_mismatches = int(sys.argv[4])
max_occurrences = int(sys.argv[5])
out = open(sys.argv[8], "w")

seqs = open(sys.argv[6], "r")
seqs_dict = {}
for line in seqs:
    if line.startswith(">"):
        header = line.strip()
        seqs_dict[header] = ""
    else:
        seqs_dict[header] += line.strip()


def check_clamp(primer, reverse):
    """Function to test if the putative primer has a GC clamp"""

    if reverse == True:
        clamp = primer[-2:]
    else:
        clamp = primer[0:2]
    if clamp == "GG" or clamp == "GC" or clamp == "CG" or clamp == "CC":
        return True


def check_Tm(GC, primer):
    """Function to calculate Tm and determine if it is acceptable"""

    Tm = 64.9 + ((41 * (GC - 16.4)) / (len(primer)))
    if Tm >= low_tm and Tm <= hi_tm:
        return Tm, True
    else:
        return None, False

def check_targets(primer):
    """Function to determine how many times the putative primer occurs in the
       full assembly"""

    temp_blastn = open("temp_blastn", "r")
    targets = 0
    for line in temp_blastn:
        if line.startswith(" Identities ="):
            regex = re.match(" Identities = ([0-9]*)/([0-9]*) .*?\n", line)
            if len(primer) - int(regex.group(1)) <= min_mismatches:
                targets += 1
    temp_blastn.close()
    return targets

# This runs through all possible primers and runs the above functions to see
# if they are good
for item in seqs_dict:
    forward = seqs_dict[item][0:seq_range].upper()
    reverse = seqs_dict[item][-seq_range:].upper()
    print "Finding primers for %s" % item
    for y in range(18, 28):
        print "Trying primers of length %s" % y
        for x in range(0, (seq_range - y)):
            f_primer = forward[x:x + y]
            if check_clamp(f_primer, False) == True:
                GC = f_primer.count("G") + f_primer.count("C")
                Tm, good_tm = check_Tm(GC, f_primer)
                if good_tm == True:
                    temp_query = open("temp_query", "w")
                    temp_query.write(">primer\n")
                    temp_query.write(f_primer)
                    temp_query.close()
                    blastn = "blastn -query temp_query -db %s -word_size 7 -out\
                              temp_blastn" % sys.argv[7]
                    subprocess.call(blastn, stdin=None, stdout=None,
                                            stderr=None, shell = True)
                    targets = check_targets(f_primer)
                    if targets <= max_occurrences:
                        f_primer_rc = Seq(f_primer)
                        f_primer_rc = f_primer_rc.reverse_complement()
                        out.write("%s_%sF Tm: %s GC: %s Occurrences in\
                                  assembly: %s\n%s\n" % (item, x, Tm,
                                  (float(GC) / y) * 100, targets, f_primer_rc))

        for x in range(0, (seq_range - y)):
            r_primer = reverse[x:x + y]
            if check_clamp(r_primer, True) == True:
                GC = r_primer.count("G") + r_primer.count("C")
                Tm, good_tm = check_Tm(GC, r_primer)
                if good_tm == True:
                    temp_query = open("temp_query", "w")
                    temp_query.write(">primer\n")
                    temp_query.write(r_primer)
                    temp_query.close()
                    blastn = "blastn -query temp_query -db %s -word_size 7\
                    -out temp_blastn" % sys.argv[7]
                    subprocess.call(blastn, stdin=None, stdout=None,
                                            stderr=None, shell = True)
                    targets = check_targets(r_primer)
                    if targets <= max_occurrences:
                        out.write("%s_%sR Tm: %s GC: %s Occurrences in\
                                  assembly: %s\n%s\n" % (item, seq_range - x,
                                  Tm, (float(GC) / y) * 100, targets, r_primer))

# A bit of cleanup
rm = "rm temp_query temp_blastn"
subprocess.call(rm, stdin=None, stdout=None, stderr=None, shell = True)
