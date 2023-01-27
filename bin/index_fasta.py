#!/usr/bin/env python3
#
# Create an OGM index for a FASTA file and a given motif
#

import sys
import ogmio


def revcomp(s):
    ret = ""
    for b in s.upper()[::-1]:
        if b == 'A':
            ret += 'T'
        elif b == 'T':
            ret += 'A'
        elif b == 'C':
            ret += 'G'
        elif b == 'G':
            ret += 'C'
        else:
            ret += 'N'
    return ret



def index_fasta(fasta: str, motifs: list[str], min_match=9, min_length=150000):
    log = ogmio.Logger()

    motif_rcs = []
    max_motif_len = len(motifs[0])
    for i, motif in enumerate(motifs):
        max_motif_len = max(max_motif_len, len(motifs[0]))
        motif_rc = revcomp(motif)
        if (motif_rc != motif):
            motif_rcs.append(motif_rc)
        else:
            motif_rcs.append(None)

        sys.stdout.write('#motif %s: %s/%s\n' % (i+1, motif, motif_rc))
        log.write('motif %s: %s/%s\n' % (i+1, motif, motif_rc))

    f = open(fasta, 'rt')
    
    cur_ref = None
    cur_pos = 0
    last_pos = 0
    buf = ""

    cur_matches = []

    last_log = 0

    while True:
        try:
            line = next(f)
        except StopIteration:
            break
        if line[0] == '>':
            if cur_ref:
                if len(cur_matches) >= min_match and cur_pos > min_length:
                    log.write('>%s\t%s\t%s\n' % (cur_ref, cur_pos, len(cur_matches)), True)
                    sys.stdout.write('>%s\t%s\t%s\n' % (cur_ref, cur_pos, len(cur_matches)))
                    for i, (pos, motif_num) in enumerate(cur_matches):
                        # write this out as a 1-based index...
                        sys.stdout.write('%s\t%s\t%s\n' % (i+1, pos+1, motif_num+1))
                else:
                    log.write('>%s\t%s\t%s\tskipped -- too short or too ref labels\n' % (cur_ref, cur_pos, len(cur_matches)), True)
                    
            cur_ref = line[1:].strip().split(' ')[0]
            cur_matches = []
            cur_pos = 0
            last_pos = 0
            last_log = 0
            buf = ""
            log.write('>%s' % cur_ref)

        else:
            buf += line.strip().upper()

        if  cur_pos - last_log > 100000:
            log.write('>%s %s' % (cur_ref, cur_pos), True)
            last_log = cur_pos

        while len(buf) >= max_motif_len:
            for i, motif in enumerate(motifs):
                motif_rc = motif_rcs[i]

                if buf[:len(motif)] == motif:
                    cur_matches.append((cur_pos, i))
                    last_pos = cur_pos
                elif motif_rc and buf[:len(motif)] == motif_rc:
                    cur_matches.append((cur_pos, i))
                    last_pos = cur_pos
            
            buf = buf[1:]
            cur_pos += 1

    if cur_ref:
        if len(cur_matches) >= min_match and cur_pos > min_length:
            log.write('>%s\t%s\t%s\n' % (cur_ref, cur_pos, len(cur_matches)), True)
            sys.stdout.write('>%s\t%s\t%s\n' % (cur_ref, cur_pos, len(cur_matches)))
            for i, (pos, motif_num) in enumerate(cur_matches):
                # write this out as a 1-based index...
                sys.stdout.write('%s\t%s\t%s\n' % (i+1, pos+1, motif_num+1))
        else:
            log.write('>%s\t%s\t%s\tskipped -- too short or too ref labels\n' % (cur_ref, cur_pos, len(cur_matches)), True)

    log.close()

if __name__ == '__main__':
    motifs = []
    fasta = None

    fasta = sys.argv[1]
    for motif in sys.argv[2:]:
        motifs.append(motif.upper())

    index_fasta(fasta, motifs)