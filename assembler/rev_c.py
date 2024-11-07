import sys

def rev_c(string) -> str:
        rev_str = string[::-1]
        r_c = ""
        for el in rev_str:
            if el == 'A':
                r_c += "T"
            if el == 'C':
                r_c += "G"
            if el == 'G':
                r_c += "C"
            if el == 'T':
                r_c += "A"
        return r_c

seq_in = sys.argv[1]
print(rev_c(seq_in))