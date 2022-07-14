from weblogo import *
filepath = 'C:/Users/liza_/Desktop/sww/'


def make_logo(filename):
    f = open('{}{}'.format(filepath, filename))
    seqs = read_seq_data(f, alphabet='ACGTacgt')
    logo_data = LogoData.from_seqs(seqs)
    logo_options = LogoOptions()
    logo_options.color_scheme = classic
    logo_options.stacks_per_line = 60
    logo_format = LogoFormat(logo_data, logo_options)
    png = png_formatter(logo_data, logo_format)
    with open('{}logo_{}.png'.format(filepath, filename.split(".")[0]), 'wb') as f:
        f.write(png)


make_logo("seqs_in_intervals.fasta")
