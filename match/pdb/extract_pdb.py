import sys

pdb = sys.argv[1]
f = open(pdb,'r+')
clean_pdb = 'clean_' + pdb
nf=open(clean_pdb,'w+')

for line in f:
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
        nf.write(line)