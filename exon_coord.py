def list_protein_gene(filepath):
    file1 = open(filepath, "r")
    file2 = open("genes.txt", 'w')
    lines = file1.readlines()
    arr = set()
    for line in lines:
        if line[0] != "#":
            components_line = line.split("\t")
            if components_line[2] == "gene":
                if "gene_biotype=protein_coding" in line:
                    components_part_eight = components_line[8].split(";")
                    geneID = components_line[0] + "_" + components_part_eight[1].split("GeneID:")[1].split(",")[0]
                    arr.add(geneID)
                    file2.write("{0}\n".format(geneID))
    return arr


def dict2(filepath):
    file1 = open(filepath, "r")
    lines = file1.readlines()
    dict = {}
    for line in lines:
        if line[0] != "#":
            components_line = line.split("\t")
            if components_line[2] == "exon":
                components_part_eight = components_line[8].split(";")
                exonID = components_line[0] + "_" + components_part_eight[0][8::]
                start = components_line[3]
                stop = components_line[4]
                if exonID in dict:
                    dict[exonID].append([start, stop])
                else:
                    dict[exonID] = [start, stop]
    return dict


def dict3(filepath):
    file1 = open(filepath, "r")
    lines = file1.readlines()
    dict = {}
    index = 0
    for line in lines:
        if line[0] != "#":
            components_line = line.split("\t")
            if components_line[2] == "CDS":
                index += 1
                components_part_eight = components_line[8].split(";")
                CDSID = components_line[0] + "_" + components_part_eight[0][7::] + "-" + str(index)
                start = components_line[3]
                stop = components_line[4]
                if CDSID in dict:
                    dict[CDSID].append([start, stop])
                else:
                    dict[CDSID] = [start, stop]
            else:
                index = 0

    return dict