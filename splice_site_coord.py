def set_ss_coord(ss_coord, chrom, start, stop, strand):
    key1 = chrom + '_' + start + '_' + stop
    key2 = chrom + '_' + start + '_' + stop
    index1 = '_5'
    index2 = '_3'
    if strand == '-':
        index1, index2 = index2, index1
    ss_coord[key1+index1] = (int(start), int(start)+1)
    ss_coord[key2+index2] = (int(stop)-1, int(stop))


def write_data(ss_coord):
    output = open("ss_coord_dict.txt", "w")
    for k in ss_coord.keys():
        output.write("{}    {}  {}\n".format(k, ss_coord[k][0], ss_coord[k][1]))


if __name__ == "__main__":
    introns = open("introns.gff", "r").readlines()
    ss_coord = {}
    for intron in introns:
        chrom, _, _, start, stop, score, strand, phase = intron.split()
        set_ss_coord(ss_coord, chrom, start, stop, strand)
    write_data(ss_coord)


