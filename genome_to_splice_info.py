#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 15.07.2022
# @author: Aleksey Komissarov and SWW splicing team
# @contact: ad3002@gmail.com

import argparse
import gzip
from collections import defaultdict, Counter


REVCOMP_DICTIONARY = dict(zip("ATCGNatcgn~[]", "TAGCNtagcn~]["))

output_fl_ss3 = open("C:/Users/khani/IdeaProjects/SWW/ss3.fasta", "w")
output_fl_ss5 = open("C:/Users/khani/IdeaProjects/SWW/ss5.fasta", "w")

def get_revcomp(sequence):
    """Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    """
    return "".join(
        REVCOMP_DICTIONARY.get(nucleotide, "") for nucleotide in reversed(sequence)
    )


def get_opener(file_name):
    if file_name.endswith("gz"):
        return gzip.open(file_name, "rt")
    else:
        return open(file_name)


def sc_iter_fasta_brute(file_name, inmem=False, lower=False):
    """Iter over fasta file."""

    header = None
    seq = []
    with get_opener(file_name) as fh:
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute ss3 and ss5")
    parser.add_argument(
        "-o", "--output", help="Output file", required=True, default=None
    )
    parser.add_argument("-s", "--stats", help="Stats file", required=True, default=None)
    parser.add_argument(
        "-i", "--fasta", help="Input fasta file", required=True, default=None
    )
    parser.add_argument(
        "-g", "--gff", help="Input gff file", required=True, default=None
    )
    parser.add_argument(
        "-f", "--flanks", help="Flank length", required=False, default=100
    )
    args = vars(parser.parse_args())

    fasta_file = args["fasta"]
    gff_file = args["gff"]
    output_file = args["output"]
    stats_output_file = args["stats"]
    flank = int(args["flanks"])

    contig2seq = {}
    for header, seq in sc_iter_fasta_brute(fasta_file):
        chrm = header[1:].split()[0]
        contig2seq[chrm] = seq

    genes = {}
    gene2exons = defaultdict(list)
    gene2cds = defaultdict(list)

    with get_opener(gff_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            d = line.strip().split("\t")
            d[3] = int(d[3])
            d[4] = int(d[4])
            features = dict([x.split("=") for x in d[-1].split(";")])
            if "Dbxref" in features:
                features["Dbxref"] = dict(
                    [x.split(":")[:2] for x in features["Dbxref"].split(",")]
                )
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

    for gene_id in genes:
        gene2exons[gene_id].sort(key=lambda x: (x[3], x[4]))

    gene2intervals = defaultdict(list)
    for gene_id in genes:
        for d in gene2exons[gene_id]:
            d = (d[0], d[3], d[4], d[6])
            gene2intervals[gene_id].append(d)

    i = 0
    exons = []
    for gene_id in gene2intervals:
        if len(gene2intervals[gene_id]) < 2:
            continue
        if gene2intervals[gene_id][0][-1] == "+":
            for j in range(1, len(gene2intervals[gene_id])):
                a = gene2intervals[gene_id][j - 1]
                b = gene2intervals[gene_id][j]
                pre_ss5 = contig2seq[a[0]][a[2] - flank : a[2]].upper()
                ss5 = contig2seq[a[0]][a[2] : a[2] + 2].upper()
                post_ss5 = contig2seq[a[0]][a[2] + 2 : a[2] + 2 + flank].upper()
                pre_ss3 = contig2seq[b[0]][b[1] - 3 - flank : b[1] - 3].upper()
                ss3 = contig2seq[b[0]][b[1] - 3 : b[1] - 3 + 2].upper()
                post_ss3 = contig2seq[b[0]][b[1] - 3 + 2 : b[1] - 3 + 2 + flank].upper()
                exon = (gene_id, "+", pre_ss3, ss3, post_ss3, pre_ss5, ss5, post_ss5)
                output_fl_ss5.write('>\n' + pre_ss5 + ss5 + post_ss5 + '\n')
                output_fl_ss3.write('>\n' + pre_ss3 + ss3 + post_ss3 + '\n')
                exons.append(exon)
        else:
            for j in range(1, len(gene2intervals[gene_id])):
                a = gene2intervals[gene_id][j - 1]
                b = gene2intervals[gene_id][j]
                pre_ss3 = contig2seq[a[0]][a[2] - flank : a[2]].upper()
                ss3 = contig2seq[a[0]][a[2] : a[2] + 2].upper()
                post_ss3 = contig2seq[a[0]][a[2] + 2 : a[2] + 2 + flank].upper()
                pre_ss5 = contig2seq[b[0]][b[1] - 3 - flank : b[1] - 3].upper()
                ss5 = contig2seq[b[0]][b[1] - 3 : b[1] - 3 + 2].upper()
                post_ss5 = contig2seq[b[0]][b[1] - 3 + 2 : b[1] - 3 + 2 + flank].upper()
                data = (pre_ss3, ss3, post_ss3, pre_ss5, ss5, post_ss5)
                exon = [gene_id, "-"] + [x for x in map(get_revcomp, data)]
                output_fl_ss5.write('>\n' + get_revcomp(pre_ss5 + ss5 + post_ss5) + '\n')
                output_fl_ss3.write('>\n' + get_revcomp(pre_ss3 + ss3 + post_ss3) + '\n')
                exons.append(tuple(exon))

    exons = set(exons)

    with open(output_file, "w") as fw:
        for exon in exons:
            fw.write("%s\n" % "\t".join(exon))

    with open(stats_output_file, "w") as fw:
        for k, v in Counter([x[6] for x in exons]).most_common(10):
            fw.write(f"ss5\t{k}\t{round(100.*v/len(exons),2)}\n")
        for k, v in Counter([x[3] for x in exons]).most_common(10):
            fw.write(f"ss3\t{k}\t{round(100.*v/len(exons),2)}\n")
