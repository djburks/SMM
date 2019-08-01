import glob
import sys
import os

order = sys.argv[3]
g_dir = 'Genomes/'
genomes = glob.glob(g_dir + '*')

read_file = sys.argv[1]
out_file = sys.argv[2]

a = open(out_file,'w')

with open(read_file) as infile:
    for lines in infile:
        if lines.startswith('>'):
            lines = lines.rstrip()
            a.write('\t' + lines[1:])
a.write('\n')
            
a.close()

for g in genomes:
    os.system('./wgm4 ' + g + ' ' + order + ' ' + read_file + ' >>' + out_file)
