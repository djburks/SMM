import sys

gfile = sys.argv[1]

metafile = {}

def fastaImporter(gfile):
    nucleo = set('ATGC')
    with open(gfile) as infile:
        for lines in infile:
            lines = lines.rstrip()
            if lines.startswith('>'):
                print(lines)
            else:
                lines = ''.join(filter(nucleo.__contains__,lines))
                print(lines)

fastaImporter(gfile)
