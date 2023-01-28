#!/usr/bin/env python3
#
# Convert an OGM index to BED format (for visualization)
#

import sys
import gzip
import math


def aln_to_bed(aln_fname: str, bp_per_pixel=1):
    f = gzip.open(aln_fname, 'rt')
    motifs = []
    for line in f:
        if not line.strip():
            continue
        
        if line.startswith("#motif"):
            spl = line.strip().split(':')
            num = int(spl[0].split(' ')[1])
            motif = spl[1].split(';')[0].strip()

            while len(motifs) < num:
                motifs.append(None)
            
            motifs[num-1] = motif
            continue

        cols = line.strip('\n').split('\t')
        molecule_id = cols[0]
        flags = int(cols[1])
        ref = cols[2]
        start = int(cols[3])-1
        end = int(cols[4])
        strand = cols[5]
        score = cols[6]
        cigar = cols[7]
        labels = cols[8].split(',')
        ref_pos = [int(x) for x in cols[9].split(',')]

        blocks = []

        l_idx = 0
        p_idx = 0

        for op in expand_cigar(cigar):
            # only possible options are S, M, D, I
            # the only ones that are represented in the ref_pos values are M and I
            if op == 'M' or op == 'I':
                blocks.append(ref_pos[p_idx])
                p_idx += 1

        start = math.floor(start - (bp_per_pixel/2))
        end = math.floor(end + (bp_per_pixel/2))

        out = [ref, start, end, molecule_id, score, strand, start, end, '0,0,0', len(blocks)]
        out.append(','.join([str(bp_per_pixel),] * len(blocks)))
        out.append(','.join([str(math.floor(x-1-start-(bp_per_pixel/2))) for x in blocks]))

        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in out]))


def expand_cigar(cigar):
    ret = ''
    n_buf = ''
    while cigar:
        if cigar[0] in '0123456789':
            n_buf += cigar[0]
        else:
            op = cigar[0]
            count = int(n_buf)
            for i in range(count):
                ret += op
            n_buf = ''
        cigar = cigar[1:]

    return ret
        



if __name__ == '__main__':
    fname = None
    fname = sys.argv[1]

    aln_to_bed(fname, 375)