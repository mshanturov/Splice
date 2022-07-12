from Bio import SeqIO


def complementary(s):
    s = s.lower()
    s.split()
    arr = s[::-1]
    dict = {'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    compl_string = ''
    for el in arr:
        compl_string += dict[el]
    return compl_string


def get_one_sub_sequence(chromosome, start, stop):
    record = SeqIO.read("chromosomes/{}.fna".format(chromosome), "fasta")
    return str(record[start - 1:stop].seq)


def set_ss_seq(ss_seq, chrom, start, stop, strand):
    key1 = chrom + '_' + start + '_' + stop
    key2 = chrom + '_' + start + '_' + stop
    index1 = '_5'
    index2 = '_3'
    s1_seq = get_one_sub_sequence(chrom, int(start), int(start)+1)
    s2_seq = get_one_sub_sequence(chrom, int(stop)-1, int(stop))
    if strand == '-':
        s1_seq = complementary(s1_seq)
        s2_seq = complementary(s2_seq)
        index1, index2 = index2, index1
    ss_seq[key1+index1] = s1_seq
    ss_seq[key2+index2] = s2_seq
    return ss_seq


def write_data(ss_seq):
    output = open("ss_seq_dict.txt", "w")
    for k in ss_seq.keys():
        output.write("{}    {}\n".format(k, ss_seq[k]))


if __name__ == "__main__":
    introns = open("introns_for_NC_000001.11_chromosome.gff", "r").readlines()
    ss_seq = {}
    for intron in introns:
        chrom, _, _, start, stop, score, strand, phase = intron.split()
        set_ss_seq(ss_seq, chrom, start, stop, strand)
    write_data(ss_seq)

