import glob
import sys
import os

order = sys.argv[3]
g_dir = 'Genomes/'
genomes = glob.glob(g_dir + '*')
smm = './smm'


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
    
b = open('Genomes.txt','w')
for g in genomes:
    b.write(g + '\n')
b.close()

os.system(smm + ' Genomes.txt ' + read_file + ' ' + order + ' >> ' + out_file)
