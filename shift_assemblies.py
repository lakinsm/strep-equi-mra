#!/usr/bin/env python3

import sys

OLD = "ACGT"
NEW = "TGCA"
TAB = str.maketrans(OLD, NEW)


def load_blastn_outfmt7(infile):
    ret = []
    with open(infile, 'r') as f:
        data = f.read().split('\n')
        for line in data:
            if not line or line[0] == '#':
                continue
            entries = line.split('\t')
            ret.append(entries)
    return ret


def fasta_parse(infile):
    with open(infile, 'r') as fasta_file:
        # Skip whitespace
        while True:
            line = fasta_file.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            all_lines = []
            line = fasta_file.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                all_lines.append(line.rstrip())
                line = fasta_file.readline()
            yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration


def circular_shift(alignments, assembly):
    genomic_contig = ''
    for ali in alignments:
        tc = int(ali[0]) - 1  # target contig
        frame = int(ali[-1].split('/')[-1])
        if frame == 1:
            start_idx = int(ali[3]) - 1
            assembly[tc][1] = assembly[tc][1][start_idx:] + assembly[tc][1][:start_idx]
        elif frame == -1:
            start_idx = len(assembly[tc][1]) - int(ali[3])
            assembly[tc][1] = assembly[tc][1].translate(TAB)[::-1]
            assembly[tc][1] = assembly[tc][1][start_idx:] + assembly[tc][1][:start_idx]
        else:
            raise ValueError('Frame must be 1 or -1: {}'.format(frame))
        genomic_contig = assembly[tc][1]
        print('Frame: {}\t{}'.format(frame, assembly[tc][1][:150]))
    return assembly, genomic_contig


def write_shifted_assembly(assembly, genome, outfile, genomic_file, samplename):
    with open(outfile, 'w') as out:
        for header, seq in assembly:
            out.write('>{}\n{}\n'.format(header, seq))
    with open(genomic_file, 'a') as gen:
        gen.write('>{} Genomic Contig\n{}\n'.format(samplename, genome))


if __name__ == '__main__':
    alignments = load_blastn_outfmt7(sys.argv[2])
    assembly = [[k, v] for k, v in fasta_parse(sys.argv[3])]
    shifted, genome = circular_shift(alignments, assembly)
    write_shifted_assembly(shifted, genome, sys.argv[4], sys.argv[5], sys.argv[1])
