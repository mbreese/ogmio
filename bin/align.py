#!/usr/bin/env python3

import sys
import ogmio
import math


class Aligner(object):
    def __init__(self, target: list[tuple[int, int]], bp_per_pixel: int = 1, match_score = 10, mismatch_score=-10, match_distance_weight=0.2, motif_mismatch_penalty = -10, missing_label_penalty = -10, extra_label_penalty = -10):
        self._match_score = match_score
        self._mismatch_score = mismatch_score
        self._match_distance_weight = match_distance_weight
        self._motif_mismatch_penalty = motif_mismatch_penalty
        self._missing_label_penalty = missing_label_penalty
        self._extra_label_penalty = extra_label_penalty

        # setup the target gap list
        last_pixel_t = 0
        last_motif_t = -1
        last_pos_t = 0
        self.t_gaps: list[tuple[int,int,int,int]] = []

        for i, (pos, motif) in enumerate(target):
            pixel = math.floor(pos / bp_per_pixel)

            if last_pixel_t:
                v = (pixel-last_pixel_t, last_motif_t, last_pos_t, pos)
                if len(self.t_gaps) == 0 or self.t_gaps[-1] != v:
                    self.t_gaps.append(v)

            last_pixel_t = pixel
            last_motif_t = motif
            last_pos_t = pos

        self.t_gaps.append((-1, last_motif_t, last_pos_t, -1))


    def align(self, query: list[tuple[int, int]], bp_per_pixel: int = 1, revcomp = False, stretch_factor = 1.0, qname=None, tname=None):

        # We will align based on the distance between labels/motifs. This means that we will be aligning based on 1-based indexing
        # but used 0-based indexing for the query/target positions. Effectively, we will have one less gap than label in the molecule.
        # (and then the genome posistions are 1-based by default, so... confusing)

        # q_gaps and t_gaps will be a list of tuples
        # (scaled_distance, label_num, scaled_left_pos, scaled_right_pos)

        q_gaps: list[tuple[int,int,int,int]] = []

        last_raw_pos = 0
        last_scaled_pos = 0
        last_pixel = 0
        last_motif = -1

        for i, (pos, motif) in enumerate(query):
            dist = pos - last_raw_pos
            scaled_dist = dist * stretch_factor
            scaled_pos = math.floor(last_scaled_pos + scaled_dist)
            pixel = math.floor(scaled_pos / bp_per_pixel)

            if last_pixel:
                v = (pixel-last_pixel, last_motif, last_scaled_pos, scaled_pos)
                # if len(q_gaps) == 0 or q_gaps[-1] != v:
                q_gaps.append(v)

            last_pixel = pixel
            last_raw_pos = pos
            last_scaled_pos = scaled_pos
            last_motif = motif


        if revcomp:
            # for revcomp, we just need to reverse the order of the gaps. The distances are still accurate,
            # but we are looking for the inverted order.
            q_gaps = q_gaps[::-1]

        # this marks the end of the gaps
        q_gaps.append((-1, last_motif, last_scaled_pos, -1))


        # setup matrix (query in rows[i], target in cols[j])
        m = []
        best_score = -1
        best_i = -1
        best_j = -1

        # calculate the matrix
        for i in range(len(q_gaps)):
            m.append([])
            for j in range(len(self.t_gaps)):
                left_score = 0
                up_score = 0
                diag_score = 0

                if i > 0:
                    # if the above pixel is the same (zero distance, same motif), merge, otherwise, it's an insert
                    if q_gaps[i][0] == 0 and q_gaps[i-1][1] == q_gaps[i][1]:
                        up_score = m[i-1][j][0]
                    else:
                        up_score = m[i-1][j][0] + self._extra_label_penalty
                if j > 0:
                    # if the left pixel is the same (zero distance, same motif), merge, otherwise, it's a deletion
                    if self.t_gaps[i][0] == 0 and self.t_gaps[j-1][1] == self.t_gaps[j][1]:
                        left_score = m[i][j-1][0]
                    else:
                        left_score = m[i][j-1][0] + self._missing_label_penalty
                if i > 0 and j > 0:
                    # if we have the same pixels, then merge, otherwise look for match/mismatch
                    if q_gaps[i][0] == 0 and q_gaps[i-1][1] == q_gaps[i][1] and self.t_gaps[i][0] == 0 and self.t_gaps[j-1][1] == self.t_gaps[j][1]:
                        diag_score = m[i-1][j-1][0]
                    else:
                        diag_score = m[i-1][j-1][0] + self.score_pos(q_gaps[i-1][0], self.t_gaps[j-1][0], q_gaps[i-1][1], self.t_gaps[j-1][1])

                scores = [(diag_score, 'd'), (up_score, 'u'), (left_score, 'l')]
                scores = sorted(scores)

                # take the highest score as the alignment score (up, left, diag)
                # this is an if/elif block to promote diagonals scores in the event of a tie

                if diag_score >= up_score and diag_score >=left_score:
                    # prefer match over indel
                    best = (diag_score, 'd')
                elif up_score >= left_score:
                    # prefer deletion over insertion
                    best = (up_score, 'u')
                else:
                    best = (left_score, 'l')

                if (best[0] > best_score):
                    best_score = best[0]
                    best_i = i
                    best_j = len(m[i])

                m[i].append(best)

        # backtrack the matrix (use SW local alignment, so backtrack from the best score anywhere)
        i = best_i
        j = best_j

        qmatch = [q_gaps[i]]
        tmatch = [self.t_gaps[j]]

        dir = m[i][j][1]
        cigar = ''
        if dir == 'd':
            cigar = 'M' + cigar
        elif dir == 'u':
            cigar = 'I' + cigar
        elif dir == 'l':
            cigar = 'D' + cigar

        while m[i][j][0] > 0:
            dir = m[i][j][1]
            m[i][j] = (m[i][j][0], '*')

            if dir == 'd':
                i -= 1
                j -= 1

                # yes, do this after the decrement
                cigar = 'M' + cigar
                qmatch.insert(0, q_gaps[i])
                tmatch.insert(0, self.t_gaps[j])

            elif dir == 'u':
                i -= 1
                cigar = 'I' + cigar
                qmatch.insert(0, q_gaps[i])

            elif dir == 'l':
                j -= 1
                cigar = 'D' + cigar
                tmatch.insert(0, self.t_gaps[j])
            else:
                sys.stderr.write("Bad matrix!\n")
                sys.exit(1)

        m[best_i][best_j] = (m[best_i][best_j][0], '***')

        cigar = self.simplify_cigar(cigar)


        # self.write_matrix(q_gaps, t_gaps, m, query, target, bp_per_pixel, stretch_factor)

        # sys.stdout.write('Best score: %s (%s, %s)\n' % (best_score, best_i, best_j))
        # sys.stdout.write('%s, %s, %s %s\n' % (query[best_i], target[best_j], best_score, cigar))
        # sys.stdout.write('[%s] %s\n' % (len(qmatch), ','.join([str(x) for x in qmatch])))
        # sys.stdout.write('[%s] %s\n' % (len(tmatch), ','.join([str(x) for x in tmatch])))

        # the aligned positions will be offset
        # the first query label will be 0
        # all other query positions will be the scaled position, offset by the first pos.
        # finally, the pos will be set relative to the first ref. position

        q_offset = qmatch[0][2] # the first pos
        if revcomp:
            q_offset = qmatch[-1][2]

        t_offset = tmatch[0][2] # the first pos
        
        aligned_pos = []
        for dist, label, start, end in qmatch:
            aligned_pos.append(start - q_offset + t_offset)

        # sys.stdout.write('[%s] %s\n' % (len(aligned_pos), ','.join([str(x) for x in aligned_pos])))

        return Alignment(qname, tname, query, aligned_pos, best_score, cigar)




        # sys.stdout.write('Query : %s\n' % len(q_gaps))
        # sys.stdout.write('Target: %s\n' % len(t_gaps))


        #sys.stdout.write("Query: %s\n" % (','.join([str(x) for x in query[:10]])))
        return None

    def simplify_cigar(self, cigar):
        ret = ''
        cur = None
        acc = 0
        for c in cigar:
            if c == cur:
                acc += 1
            else:
                if cur:
                    ret += '%s%s' % (acc, cur)
                cur = c
                acc = 1
        if cur:
            ret += '%s%s' % (acc, cur)

        return ret


    def score_pos(self, q_gap, t_gap, q_motif, t_motif):
        if q_motif != t_motif:
            return self._motif_mismatch_penalty
        
        return max(0, (1 - (abs(q_gap - t_gap) * self._match_distance_weight)) * self._match_score) 


    def write_matrix(self, q_gaps, t_gaps, m, query, target, bp_per_pixel, stretch_factor):
        sys.stdout.write('\t')
        last_pixel = 0
        for j in range(len(m[0])):
            if j > 0:
                pixel = math.floor(target[j-1][0] / bp_per_pixel)
                #v = math.floor((pixel-last_pixel) * stretch_factor) #you don't stretch the ref pos
                interval = math.floor(pixel-last_pixel)
                sys.stdout.write('%s/%s\t' % (math.floor(last_pixel + interval), t_gaps[j-1][0]))
                last_pixel = last_pixel + interval

            else:
                sys.stdout.write('\t')
        sys.stdout.write('\n')
        last_pixel = 0
        for i in range(len(m)):
            if i > 0:
                pixel = math.floor(query[i-1][0] / bp_per_pixel)
                interval = math.floor((pixel-last_pixel) * stretch_factor)
                sys.stdout.write('%s/%s\t' % (last_pixel + interval, q_gaps[i-1][0])) 
                last_pixel = last_pixel + interval

            else:
                sys.stdout.write('\t')

            for j in range(len(m[i])):
                sys.stdout.write(('%.1f %s\t' % m[i][j]) if m[i][j] else '') 
            sys.stdout.write('\n')


class Alignment(object):
    def __init__(self, qname, tname, q_orig_label, q_aligned_pos, score, cigar):
        self._qname = qname
        self._tname = tname
        self._q_orig_label = q_orig_label
        self._q_aligned_pos = q_aligned_pos
        self._score = score
        self._cigar = cigar

    def write(self, out=sys.stdout):
        cols = [self._qname, self._tname, 
                self._score, self._cigar, 
                ','.join(['%s:%s' % (x[0], x[1]) for x in self._q_orig_label]), 
                ','.join([str(x) for x in self._q_aligned_pos])
                ]
        out.write('%s\n' % '\t'.join([str(x) for x in cols]))


def main(ref_fname, bnx_fname, use_stretch=True):
    idx = ogmio.RefIndex(ref_fname)
    bnx = ogmio.BNXFile(bnx_fname)

    aligners: dict[tuple[str,int], Aligner] = {}

    log = ogmio.Logger()

    i = 1
    for mol in bnx.molecules():
        if mol:
            #sys.stdout.write(repr(mol))
            #sys.stdout.write('===============\n')

            best_aln = None

            mol_all =  mol.get_all_labels()

            log.write('[%s/%s] %s:' % (i, mol._header._num_of_molecules,  mol._molecule_id), True)
            #sys.stderr.write("Labels  : %s\n" % (','.join([str(x) for x in all])))

            for ref in idx.get_refs():
                k = (ref, mol.get_run()._bases_per_pixel)
                if not k in aligners:
                    aligners[k] = Aligner(idx.get_ref_labels(ref), mol.get_run()._bases_per_pixel)
                
                aligner = aligners[k]


                log.write('[%s/%s] %s: %s' % (i, mol._header._num_of_molecules, mol._molecule_id, ref), True)
                # if ref != 'chr12':
                #     continue

                aln1 = aligner.align(mol_all, bp_per_pixel=mol.get_run()._bases_per_pixel, revcomp=False, stretch_factor=mol.get_run()._stretch_factor if use_stretch else 1, qname=mol._molecule_id, tname=ref)
                aln2 = aligner.align(mol_all, bp_per_pixel=mol.get_run()._bases_per_pixel, revcomp=True, stretch_factor=mol.get_run()._stretch_factor if use_stretch else 1, qname=mol._molecule_id, tname=ref)

                if not best_aln or best_aln._score < aln1._score:
                    best_aln = aln1

                if not best_aln or best_aln._score < aln2._score:
                    best_aln = aln2

            if best_aln:
                log.write('[%s/%s] %s\t%s\t%s\t%s\n' % (i, mol._header._num_of_molecules, mol._molecule_id, best_aln._tname, best_aln._score, best_aln._cigar), True)
                best_aln.write()

        i += 1
    log.close()


if __name__ == '__main__':
    use_stretch = True
    ref_fname = None
    bnx_fname = None

    for arg in sys.argv[1:]:
        if arg == '--no-stretch':
            use_stretch = False

        elif not ref_fname:
            ref_fname = arg
        elif not bnx_fname:
            bnx_fname = arg

    main(ref_fname, bnx_fname, use_stretch=use_stretch)
