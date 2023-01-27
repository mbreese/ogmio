import gzip
import sys
import math

class BNXRun(object):
    def __init__(self, folder, instrument_serial, time, pixels_per_scan, stretch_factor, bases_per_pixel, num_scans, chip_id, flowcell, snr_filter, min_length, min_snr, run_id):
        self._folder=folder
        self._instrument_serial=instrument_serial
        self._time=time
        self._pixels_per_scan=pixels_per_scan
        self._stretch_factor=stretch_factor
        self._bases_per_pixel=bases_per_pixel
        self._num_scans=num_scans
        self._chip_id=chip_id
        self._flowcell=flowcell
        self._snr_filter=snr_filter
        self._min_length=min_length
        self._min_snr=min_snr
        self._run_id=run_id        


class BNXHeader(object):
    #
    # rh SourceFolder    InstrumentSerial        Time    NanoChannelPixelsPerScan        StretchFactor   BasesPerPixel   NumberofScans   ChipId  Flowcell \
    #    LabelSNRFilterType      MinMoleculeLength       MinLabelSNR1    RunId
    #
    def __init__(self, fname):
        f = gzip.open(fname, 'rt')
        self._runs = {}
        self._molecule_header_idx = {}
        self._molecule_header_vals = []
        self._label_count = -1
        self._labels = []
        self._quality_scores = {}
        self._num_of_molecules = -1

        header_idx = {}
        for line in f:
            if line[0] == '#':
                if line.startswith("#rh"):
                    cols = line[3:].strip().split('\t')
                    for i,v in enumerate(cols):
                        header_idx[v] = i
                    continue

                if header_idx and line.startswith("# Run Data"):
                    cols = line[10:].strip().split('\t')

                    folder = cols[header_idx["SourceFolder"]]
                    instrument_serial = cols[header_idx["InstrumentSerial"]]
                    time = cols[header_idx["Time"]]
                    pixels_per_scan = int(cols[header_idx["NanoChannelPixelsPerScan"]])
                    stretch_factor = float(cols[header_idx["StretchFactor"]])
                    bases_per_pixel = int(cols[header_idx["BasesPerPixel"]])
                    num_scans = int(cols[header_idx["NumberofScans"]])
                    chip_id = cols[header_idx["ChipId"]]
                    flowcell = int(cols[header_idx["Flowcell"]])
                    snr_filter = cols[header_idx["LabelSNRFilterType"]]
                    min_length = float(cols[header_idx["MinMoleculeLength"]])
                    min_snr = float(cols[header_idx["MinLabelSNR1"]])
                    run_id = int(cols[header_idx["RunId"]])

                    k = (chip_id, flowcell, run_id)
                    run = BNXRun(folder, instrument_serial, time, pixels_per_scan, stretch_factor, bases_per_pixel, num_scans, chip_id, flowcell, snr_filter, min_length, min_snr, run_id)

                    self._runs[k] = run

                if line.startswith("# Label Channels"):
                    self._label_count = int(line.split('\t')[1])
                    while len(self._labels) < self._label_count:
                        self._labels.append(None)

                if self._label_count > 0 and line.startswith("# Nickase Recognition Site"):
                    num = int(line.split(":")[0].split(" ")[-1])
                    label = line.split(":")[1].split(";")[0].strip()
                    label_desc = line.split(":")[1].split(";")[1].strip()
                    self._labels[num-1] = (label, label_desc)

                if line.startswith("# Number of Molecules:"):
                    self._num_of_molecules = int(line.split(":")[1].strip())

                if line.startswith("# Quality Score"):
                    scoreID = line.split(":")[0].split(" ")[-1]
                    scoreDesc = line.split(":")[1].strip()
                    self._quality_scores[scoreID] = scoreDesc

                if line.startswith("#0h"):
                    cols = line[3:].strip().split('\t')
                    for i,v in enumerate(cols):
                        self._molecule_header_idx[v] = i
                        self._molecule_header_vals.append(v)

            else:
                break
        f.close()


    def get_molecule_header_byidx(self, idx):
        return self._molecule_header_vals[idx]


    def get_run(self, chip_id: str, flowcell: int, run_id: int):
        k = (chip_id, flowcell, run_id)
        if k in self._runs:
            return self._runs[k]

        sys.stderr.write("Missing run: %s\n" % str(k))
        return None


class Molecule(object):
    def __init__(self, header: BNXHeader, molecule_id, length, ave_intesity, snr, num_labels, 
                       orig_molecule_id, scan_num, scan_dir, chip_id, flowcell, run_id, column, start_fov, 
                       start_x, start_y, end_fov, end_x, end_y, labels: dict[int, list[int]], qual_scores: dict[str, list[float]]):

        # #0h	LabelChannel	MoleculeId	Length	AvgIntensity	SNR	NumberofLabels	OriginalMoleculeId	ScanNumber	
        #       ScanDirection	ChipId	Flowcell	RunId	Column	StartFOV	StartX	StartY	EndFOV	EndX	EndY    

        self._header = header
        self._molecule_id = molecule_id
        self._length = length
        self._ave_intesity = ave_intesity
        self._snr = snr
        self._num_labels = num_labels
        self._orig_molecule_id = orig_molecule_id
        self._scan_num = scan_num
        self._scan_dir = scan_dir
        self._chip_id = chip_id
        self._flowcell = flowcell
        self._run_id = run_id
        self._column = column
        self._start_fov = start_fov
        self._start_x = start_x
        self._start_y = start_y
        self._end_fov = end_fov
        self._end_x = end_x
        self._end_y = end_y

        self._labels = labels
        self._qual_scores = qual_scores
        pass



    def __repr__(self):
        ret = ''
        ret += "MoleculeId: %s\n" % (self._molecule_id)
        ret += "Length: %s\n" % (self._length)
        ret += "AvgIntensity: %s\n" % (self._ave_intesity)
        ret += "SNR: %s\n" % (self._snr)
        ret += "NumberofLabels: %s\n" % (self._num_labels)
        ret += "OriginalMoleculeId: %s\n" % (self._orig_molecule_id)
        ret += "ScanNumber: %s\n" % (self._scan_num)
        ret += "ScanDirection: %s\n" % (self._scan_dir)
        ret += "ChipId: %s\n" % (self._chip_id)
        ret += "Flowcell: %s\n" % (self._flowcell)
        ret += "RunId: %s\n" % (self._run_id)
        ret += "Column: %s\n" % (self._column)
        ret += "StartFOV: %s\n" % (self._start_fov)
        ret += "StartX: %s\n" % (self._start_x)
        ret += "StartY: %s\n" % (self._start_y)
        ret += "EndFOV: %s\n" % (self._end_fov)
        ret += "EndX: %s\n" % (self._end_x)
        ret += "EndY: %s\n" % (self._end_y)

        for label_num in self._labels:
            ret += "Label %s: %s\n" % (label_num, ','.join([str(x) for x in self._labels[label_num]]))
        for qs in self._qual_scores:
            ret += "Qual scores %s: %s\n" % (qs, ','.join([str(x) for x in self._qual_scores[qs]]))

        return ret


    def get_label_pos(self, motif_num):
        return self._labels[motif_num]


    def get_all_labels(self):
        ret = []
        for i in self._labels:
            for pos in self._labels[i]:
                ret.append((pos, i))

        return sorted(ret)


    def get_run(self):
        return self._header.get_run(self._chip_id, self._flowcell, self._run_id)


    def get_all_labels_pixels(self, bp_per_pixel):
        ret = []
        for i in self._labels:
            for pos in self._labels[i]:
                pixel = math.floor(pos/bp_per_pixel)
                ret.append((pixel, i))

        ret = sorted(ret)
        ret2 = []
        last = None
        for v in ret:
            if not last or last != v:
                ret2.append(v)
            last = v
        return ret2




    def get_qual_score(self, qs):
        return self._qual_scores[qs]


# LabelChannel
# MoleculeId
# Length
# AvgIntensity
# SNR
# NumberofLabels
# OriginalMoleculeId
# ScanNumber
# ScanDirection
# ChipId
# Flowcell
# RunId
# Column
# StartFOV
# StartX
# StartY
# EndFOV
# EndX
# EndY 

# label_channel
# molecule_id
# length
# ave_intesity
# snr
# num_labels
# orig_molecule_id
# scan_num
# scan_dir
# chip_id
# flowcell
# run_id
# column
# start_fov
# start_x
# start_y
# end_fov
# end_x
# end_y

    def parse_lines(header: BNXHeader, lines: list[str]):
        ann: dict[str, str] = {}
        labels: dict[int, list[int]] = {}
        qual_scores: dict[str, list[float]] = {}
        for line in lines:
            cols = line.strip().split('\t')
            if cols[0] == '0':
                for i, v in enumerate(cols[1:]):
                    ann[header.get_molecule_header_byidx(i+1)] = v
            elif cols[0][0] == 'Q':
                qual_scores[cols[0]] = [float(x) for x in cols[1:]]
            else:
                label_num = int(cols[0])
                labels[label_num] = [int(x) for x in cols[1:]]
            
        molecule_id = int(ann['MoleculeId']) if 'MoleculeId' in ann else None
        length = float(ann['Length']) if 'Length' in ann else None
        ave_intesity = float(ann['AvgIntensity']) if 'AvgIntensity' in ann else None
        snr = float(ann['SNR']) if 'SNR' in ann else None
        num_labels = int(ann['NumberofLabels']) if 'NumberofLabels' in ann else None
        orig_molecule_id = int(ann['OriginalMoleculeId']) if 'OriginalMoleculeId' in ann else None
        scan_num = int(ann['ScanNumber']) if 'ScanNumber' in ann else None
        scan_dir = int(ann['ScanDirection']) if 'ScanDirection' in ann else None
        chip_id = ann['ChipId'] if 'ChipId' in ann else None
        flowcell = int(ann['Flowcell']) if 'Flowcell' in ann else None
        run_id = int(ann['RunId']) if 'RunId' in ann else None
        column = int(ann['Column']) if 'Column' in ann else None
        start_fov = int(ann['StartFOV']) if 'StartFOV' in ann else None
        start_x = int(ann['StartX']) if 'StartX' in ann else None
        start_y = int(ann['StartY']) if 'StartY' in ann else None
        end_fov = int(ann['EndFOV']) if 'EndFOV' in ann else None
        end_x = int(ann['EndX']) if 'EndX' in ann else None
        end_y = int(ann['EndY']) if 'EndY' in ann else None

        return Molecule(header, molecule_id, length, ave_intesity, snr, num_labels, 
                       orig_molecule_id, scan_num, scan_dir, chip_id, flowcell, run_id, column, start_fov, 
                       start_x, start_y, end_fov, end_x, end_y, 
                       labels, qual_scores)



class BNXFile(object):
    def __init__(self, fname):
        self._fname = fname
        self._header = BNXHeader(fname)

    def molecules(self):
        f = gzip.open(self._fname, 'rt')
        buf = ''
        lines = []
        for line in f:
            if not line.strip() or line[0] == '#':
                continue

            if line[0] == '0':
                if lines:
                    yield Molecule.parse_lines(self._header, lines)
                    lines = []
                    buf = ''
                lines.append(line.strip())
            else:
                lines.append(line.strip())

        if lines:
            yield Molecule.parse_lines(self._header, lines)

        f.close()


        
class RefIndex(object):
    def __init__(self, fname):
        self._motifs = {}
        self._refs = {}

        cur_ref = None
        f = gzip.open(fname, 'rt')

        for line in f:
            if not line.strip():
                continue

            if line.startswith("#motif "):
                num = int(line.strip().split(':')[0].split(' ')[1])
                motif = line.strip().split(':')[1].split('/')[0]
                self._motifs[num] = motif
                continue

            if line[0] == '>':
                spl = line[1:].strip().split('\t')
                ref = spl[0]
                self._refs[ref] = []
                cur_ref = ref
            elif cur_ref:
                spl = line.split('\t')
                self._refs[ref].append((int(spl[1]), int(spl[2])))

        f.close()

    
    def get_refs(self):
        return self._refs.keys()


    def get_motifs(self):
        return self._motifs

    
    def get_ref_labels(self, ref):
        return self._refs[ref]


    def get_ref_labels_pixels(self, ref, bp_per_pixel):
        ret = []
        for pos, motif in self._refs[ref]:
            pixel = math.floor(pos / bp_per_pixel)
            ret.append((pixel, motif))

        ret = sorted(ret)
        ret2 = []
        last = None
        for v in ret:
            if not last or last != v:
                ret2.append(v)
            last = v
        return ret2


class Logger(object):
    def __init__(self, fobj=sys.stderr):
        self._last = ''
        self._fobj = fobj

    def write(self, msg, reset=False):
        if reset and self._last:
            self._fobj.write('\r')
            for i in range(len(self._last)+2):
                self._fobj.write(' ')
            self._fobj.write('\r')

        self._fobj.write(msg)
        self._fobj.flush()

        if msg.find('\n') == -1:
            self._last = msg
        else:
            self._last = ''
    
    def close(self):
        self._fobj.write('\n')
            
