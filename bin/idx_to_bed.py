#!/usr/bin/env python3
#
# Convert an OGM index to BED format (for visualization)
#

import sys
import ogmio


def idx_to_bed(ref_fname: str):
    idx = ogmio.RefIndex(ref_fname)
    for ref in idx.get_refs():
        for label in idx.get_ref_labels(ref):
            sys.stdout.write('%s\t%s\t%s\t%s\n' % (ref, label[0]-1, label[0]-1 + len(idx.get_motif(label[1])), idx.get_motif(label[1])))



if __name__ == '__main__':
    fname = None
    fname = sys.argv[1]

    idx_to_bed(fname)