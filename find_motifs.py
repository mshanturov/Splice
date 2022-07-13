file = open("C:/Users/sizov/bioinf/splicing/splice3_50_upper.fa", 'r')
lines = file.readlines()
motif = 'TTAAT'
len_motif = len(motif)
out = open("TTAAT.fa", 'w')

header = ''
for line in lines:
    if line[0] == ">":
        header = line
    else:
        start = line.find(motif) - 49
        stop = start + len_motif - 1
        if start > -50:
            out.write("".join([header, str(start), "\t", str(stop), "\n"]))
