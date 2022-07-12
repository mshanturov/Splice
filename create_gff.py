import sys

# filepath - путь к файлу GRCh38_latest_genomic.gff
# filepath = "C:/Users/sizov/bioinf/splicing/GRCh38_latest_genomic.gff"
# python create_gff.py -i <fasta> -a <input_gff> -o <output_gff>
filepath = sys.argv[4]


def exon_param(exon, cds_coord, strand, phases):
    if strand == '+':
        for i in range(len(cds_coord)):
            if exon == cds_coord[i]:
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + phases[i] + "\tCDS"

        if int(exon[0]) < int(cds_coord[0][0]):
            if int(exon[1]) < int(cds_coord[0][0]):
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + "." + "\tUTR5"
            else:
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + phases[0] + "\tUTR5"

        if int(exon[1]) > int(cds_coord[-1][1]):
            if int(exon[0]) > int(cds_coord[-1][1]):
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + "." + "\tUTR3"
            else:
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + phases[-1] + "\tUTR3"

    if strand == '-':
        for i in range(len(cds_coord)):
            if exon == cds_coord[i]:
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + phases[i] + "\tCDS"

        if int(exon[1]) > int(cds_coord[0][1]):
            if int(exon[0]) > int(cds_coord[0][1]):
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + "." + "\tUTR5"
            else:
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + phases[0] + "\tUTR5"

        if int(exon[0]) < int(cds_coord[-1][0]):
            if int(exon[1]) < int(cds_coord[-1][0]):
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + "." + "\tUTR3"
            else:
                return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + phases[-1] + "\tUTR3"
    return exon[0] + "\t" + exon[1] + "\t.\t" + strand + "\t" + "." + "\tENCLOSED"


def create_gff_module(chr, exons_list, cds_coord, strand, phases, exon_id):
    introns_list = []
    res = ''
    if strand == "+":
        for i in range(len(exons_list) - 1):
            introns_list.append([exons_list[i][1], exons_list[i+1][0]])
        for i in range(len(introns_list)):
            #res += chr + "\t.\t" + "intron\t" + introns_list[i][0] + "\t" + introns_list[i][1] + "\t.\t+\t." + "\n"
            res += chr + "\t.\t" + "exon\t" + exon_param(exons_list[i], cds_coord, strand, phases) + ";" + exon_id[i] + "\n" +\
                   chr + "\t.\t" + "splice5\t" + introns_list[i][0] + "\t" + str(int(introns_list[i][0]) + 1) + "\t.\t+\t." + "\n" +\
                   chr + "\t.\t" + "intron\t" + introns_list[i][0] + "\t" + introns_list[i][1] + "\t.\t+\t." + "\n" +\
                   chr + "\t.\t" + "splice3\t" + str(int(introns_list[i][1]) - 1) + "\t" + introns_list[i][1] + "\t.\t+\t." + "\n"
        res += chr + "\t.\t" + "exon\t" + exon_param(exons_list[-1], cds_coord, strand, phases) + ";" + exon_id[-1] + "\n"

    if strand == "-":
        for i in range(len(exons_list) - 1):
            introns_list.append([exons_list[i+1][1], exons_list[i][0]])
        for i in range(len(introns_list)):
            #res += chr + "\t.\t" + "intron\t" + introns_list[i][0] + "\t" + introns_list[i][1] + "\t.\t-\t." + "\n"
            res += chr + "\t.\t" + "exon\t" + exon_param(exons_list[i], cds_coord, strand, phases) + ";" + exon_id[i] + "\n" +\
                   chr + "\t.\t" + "splice5\t" + str(int(introns_list[i][1]) - 1) + "\t" + introns_list[i][1] + "\t.\t-\t." + "\n" +\
                   chr + "\t.\t" + "intron\t" + introns_list[i][0] + "\t" + introns_list[i][1] + "\t.\t-\t." + "\n" +\
                   chr + "\t.\t" + "splice3\t" + introns_list[i][0] + "\t" + str(int(introns_list[i][0]) + 1) + "\t.\t-\t." + "\n"
        res += chr + "\t.\t" + "exon\t" + exon_param(exons_list[-1], cds_coord, strand, phases) + ";" + exon_id[-1] + "\n"

    return res


def write_in_file(file, chr, exons_coord, cds_coord, strand, phases, exons_id):
    file4 = open(file, 'a')
    file4.write(create_gff_module(chr, exons_coord, cds_coord, strand, phases, exons_id))


outfile = sys.argv[6]
exons_coord = []
cds_coord = []
exons_id = []
file1 = open(filepath, 'r')
lines = file1.readlines()
parent = ''
strand = ''
chr = ''
phases = []
for line in lines:
    if line[0] != "#":
        componentsLine = line.split("\t")
        tmp_parent = componentsLine[8].split(";")[1].split('=')[1][4::]
        if (componentsLine[2] == "exon" or componentsLine[2] == "CDS") and tmp_parent != parent:
            if len(cds_coord) > 0:
                write_in_file(outfile, chr, exons_coord, cds_coord, strand, phases, exons_id)
            parent = tmp_parent
            exons_coord.clear()
            cds_coord.clear()
            phases.clear()
            exons_id.clear()
        if componentsLine[2] == "exon":
            componentsPartEight = componentsLine[8].split(";")
            chr = componentsLine[0]
            exon_id = chr + "_" + componentsPartEight[0][8::]
            start = componentsLine[3]
            stop = componentsLine[4]
            strand = componentsLine[6]
            exons_coord.append([start, stop])
            exons_id.append(exon_id)
        if componentsLine[2] == "CDS":
            componentsPartEight = componentsLine[8].split(";")
            cds_id = componentsLine[0] + "_" + componentsPartEight[0][8::]
            start = componentsLine[3]
            stop = componentsLine[4]
            phase = componentsLine[7]
            cds_coord.append([start, stop])
            phases.append(phase)
