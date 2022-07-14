
genome_file = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
ggf_file = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz"

fasta_file = "C:/Users/khani/IdeaProjects/SWW/GRCh38_latest_genomic.fna"
gff_annotation_file = "C:/Users/khani/IdeaProjects/SWW/GRCh38_latest_genomic.gff"
output_fl_ss3 = open("C:/Users/khani/IdeaProjects/SWW/ss3_fl.fasta", "w")
output_fl_ss5 = open("C:/Users/khani/IdeaProjects/SWW/ss5_fl.fasta", "w")

REVCOMP_DICTIONARY = dict(zip('ATCGNatcgn~[]', 'TAGCNtagcn~]['))

def sc_iter_fasta_brute(file_name, inmem=False, lower=False):
    """ Iter over fasta file."""

    header = None
    seq = []
    with open(file_name) as fh:
        if inmem:
            data = fh.readlines()
        else:
            data = fh
        for line in data:
            if line.startswith(">"):
                if seq or header:
                    sequence = "".join(seq)
                    if lower:
                        sequence = sequence.lower()
                    yield header, sequence
                header = line.strip()
                seq = []
                continue
            seq.append(line.strip())
        if seq or header:
            sequence = "".join(seq)
            if lower:
                sequence = sequence.lower()
            yield header, sequence

#%%



def filter_name_chrmosomes_hg19(chrm):
    if not chrm.startswith("NC") or "NC_012920" in chrm:
        return True
    return False

contig2seq = {}

for header, seq in sc_iter_fasta_brute(fasta_file):
    chrm = header[1:].split()[0]
    if filter_name_chrmosomes_hg19(chrm):
        continue
    print(chrm, end=" ")
    contig2seq[chrm] = seq

#%%

from collections import defaultdict

genes = {}
gene2exons = defaultdict(list)
gene2cds = defaultdict(list)

with open(gff_annotation_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        d = line.strip().split("\t")
        d[3] = int(d[3])
        d[4] = int(d[4])
        features = dict([x.split("=") for x in d[-1].split(";")])
        if "Dbxref" in features:
            features["Dbxref"] = dict([x.split(":")[:2] for x in features["Dbxref"].split(",")])
        d[-1] = features
        if d[2] == "gene":
            if features["gene_biotype"] != "protein_coding":
                continue
            gene_id = features["Dbxref"]["GeneID"]
            genes[gene_id] = d
        if d[2] == "exon":
            gene_id = features["Dbxref"]["GeneID"]
            gene2exons[gene_id].append(d)
        if d[2] == "cds":
            gene_id = features["Dbxref"]["GeneID"]
            gene2cds[gene_id].append(d)


#%%

for gene_id in genes:
    gene2exons[gene_id].sort(key=lambda x: (x[3], x[4]))

#%%

gene2intervals = defaultdict(list)
for gene_id in genes:
    for d in gene2exons[gene_id]:
        d = (d[0], d[3], d[4], d[6])
        gene2intervals[gene_id].append(d)


#%%
def get_revcomp(sequence):
    return ''.join(REVCOMP_DICTIONARY.get(nucleotide, '') for nucleotide in reversed(sequence))

flank = 15
i = 0
exons = []
for gene_id in gene2intervals:
    if len(gene2intervals[gene_id]) < 2:
        continue
    if filter_name_chrmosomes_hg19(gene2intervals[gene_id][0][0]):
        continue
    if gene2intervals[gene_id][0][-1] == "+":
        for j in range(1, len(gene2intervals[gene_id])):
            a = gene2intervals[gene_id][j-1]
            b = gene2intervals[gene_id][j]
            pre_ss3 = contig2seq[a[0]][a[2]-flank:a[2]].upper()
            ss3 = contig2seq[a[0]][a[2]:a[2]+2].upper()
            post_ss3 = contig2seq[a[0]][a[2]+2:a[2]+2+flank].upper()
            pre_ss5 = contig2seq[b[0]][b[1]-3-flank:b[1]-3].upper()
            ss5 = contig2seq[b[0]][b[1]-3:b[1]-3+2].upper()
            post_ss5 = contig2seq[b[0]][b[1]-3+2:b[1]-3+2+flank].upper()
            data = (gene_id, "+", pre_ss3, ss3, post_ss3, pre_ss5, ss5, post_ss5)
            #
            seq_name = gene2intervals[gene_id][0][0]
            output_fl_ss3.write('>' + seq_name + ':(+)\n' + pre_ss5 + ss5 + post_ss5 + '\n')#наоборот, чтобы было правильно(????)
            output_fl_ss5.write('>' + seq_name + ':(+)\n' + pre_ss3 + ss3 + post_ss3 + '\n')
            #
            exons.append(data)
    else:
        for j in range(1, len(gene2intervals[gene_id])):
            a = gene2intervals[gene_id][j-1]
            b = gene2intervals[gene_id][j]
            pre_ss5 = contig2seq[a[0]][a[2]-flank:a[2]].upper()
            ss5 = contig2seq[a[0]][a[2]:a[2]+2].upper()
            post_ss5 = contig2seq[a[0]][a[2]+2:a[2]+2+flank].upper()
            pre_ss3 = contig2seq[b[0]][b[1]-3-flank:b[1]-3].upper()
            ss3 = contig2seq[b[0]][b[1]-3:b[1]-3+2].upper()
            post_ss3 = contig2seq[b[0]][b[1]-3+2:b[1]-3+2+flank].upper()
            data = (pre_ss3, ss3, post_ss3, pre_ss5, ss5, post_ss5)
            data = [gene_id, "-"] + [x for x in map(get_revcomp, data)]
            #
            seq_name = gene2intervals[gene_id][0][0]
            output_fl_ss3.write('>' + seq_name + ':(-)\n' + get_revcomp(pre_ss5 + ss5 + post_ss5) + '\n')
            output_fl_ss5.write('>' + seq_name + ':(-)\n' + get_revcomp(pre_ss3 + ss3 + post_ss3) + '\n')
            #
            exons.append(data)





#%%

len(exons)


#%%

from collections import Counter

for k,v in Counter([x[6] for x in exons]).items():
    print(k, round(100.*v/len(exons),2))

#%%

for k,v in Counter([x[3] for x in exons]).items():
    print(k, round(100.*v/len(exons),2))

#%%

left_seqs = []
for x in exons:
    left_seqs.append("".join(x[2:2+3]))
right_seqs = []
for x in exons:
    right_seqs.append("".join(x[2+3:2+3+3]))

#%%

N = len(right_seqs)
cons = []
for i in range(len(right_seqs[0])):
    c = [(k, round(100.*v/N,2))  for k,v in Counter([x[i] for x in right_seqs]).most_common(1)]
    cons.append((c[0][0]))
"".join(cons)

#%%

N = len(left_seqs)
cons = []
for i in range(len(left_seqs[0])):
    c = [(k, round(100.*v/N,2)) for k,v in Counter([x[i] for x in left_seqs if i < len(x)]).most_common(1)]
    cons.append((c[0][0]))
#     print(i-200, [(k, round(100.*v/N,2))  for k,v in Counter([x[i] for x in left_seqs]).most_common(1)])
"".join(cons)

#%%

motifs = Counter()
motifs_to_poses = defaultdict(Counter)
for seq in left_seqs:
    for i in range(140, 200):
        kmer = seq[i:i+5]
        motifs[kmer] += 1
        motifs_to_poses[kmer][i-200] += 1

#%%

for x in motifs.most_common(2000):
    kmer, tf = x
    x, y = motifs_to_poses[kmer].most_common(1)[0]
    max_enrich = round(100.*y/tf,1)
    if  max_enrich > 5 and x > 10:
        print(x, " ".join([f'{y[0]} : {round(100.*y[1]/tf,1)}%'
                           for y in motifs_to_poses[kmer].most_common(1)]))

#%%

motifs["CTAAC"], [f'{x[0]} : {round(100.*x[1]/motifs["CTAAC"],1)}%'
                  for x in motifs_to_poses["CTAAC"].most_common(10)]

#%%


