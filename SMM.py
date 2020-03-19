# python3 SMM.py input.fasta output.fasta 12 norm

import sys
import subprocess

order = sys.argv[3]
stype = sys.argv[4]
smm = './smm'

if (stype != "norm"):
    stype = "raw"

read_file = sys.argv[1]
out_file = sys.argv[2]

if (order != "12") and (stype == "norm"):
    print("WARNING: Normalization not compatible with order != 12.  Reverting to raw scores.")
    stype = "raw"

a = open(out_file,'w')

with open(read_file) as infile:
    for lines in infile:
        if lines.startswith('>'):
            lines = lines.rstrip()
            a.write('\t' + lines[1:])
a.write('\n')
            
a.close()


subprocess.run([smm,"Genomes.txt",read_file,order,stype],stdout=open(out_file,'a'),stderr=open('SMM_Errors.txt','w'))

