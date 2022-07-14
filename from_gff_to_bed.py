part = "intron"

gff = open("ex.gff", 'r')
bed = open(part + ".bed", 'w')
lines = gff.readlines()

for line in lines:
    if line[0] != "#":
        componentsLine = line.split("\t")
        strand = componentsLine[6]
        if componentsLine[2] == part:
            strand = componentsLine[6]
            if strand == '+':
                start = int(componentsLine[3]) - 1
                stop = int(componentsLine[4])
                bed.write("\t".join([componentsLine[0], str(start), str(stop), ".", ".", strand, "\n"]))
            if strand == '-':
                start = int(componentsLine[3])
                stop = int(componentsLine[4]) + 1
                bed.write("\t".join([componentsLine[0], str(start), str(stop), ".", ".", strand, "\n"]))
