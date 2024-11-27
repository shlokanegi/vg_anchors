import sys
import gzip
from contextlib import contextmanager

SEQ_DV_FIELD=15
MIN_SEQ_DV=0.5
GAF_FIELDS = ["Query name","Querylength","Query start","Query end","Strand","Path","Path length","Path start","Path end","# Residue matches","Alignment_block_length","MapQ","DP Alignment Score","bq","Difference sequence","Approx per-base sequence divergence"]
#dv:f:0.0008
@contextmanager
def open_fasta(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, 'rt') as f:
            yield f
    else:
        with open(filename, 'r') as f:
            yield f

def read_alignment_file(filename):
    with open(filename, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > len(GAF_FIELDS):
                if float(fields[SEQ_DV_FIELD + 2][5:]) >= MIN_SEQ_DV:
                    print(line,end="")  # Convert y values to float
            else:
                if float(fields[SEQ_DV_FIELD][5:]) >= MIN_SEQ_DV:
                    print(line,end="")  # Convert y values to float

if __name__ == "__main__":
    read_alignment_file(sys.argv[1])