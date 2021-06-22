from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from Bio.Seq import Seq

class Primer_range():
    def __init__(self, BED_FILE, REF_FILE):
        self.ranges_fwd = []
        self.ranges_rev = []
        self.amplicon = []
        self.primer_to_seq = {}
        ref_to_seq = {header:seq for header,seq in sfp(open(REF_FILE))}
        with open(BED_FILE, 'r') as bed_h:
            for bed_line in bed_h:
                ref,start,end,name,_,strand = bed_line.rstrip().split('\t')
                self.primer_to_seq[name] = ref_to_seq[ref][int(start):int(end)]
                if strand == "+":
                    self.ranges_fwd.append((int(start), int(end), strand, name))
                if strand == "-":
                    self.ranges_rev.append((int(start), int(end), strand, name))
        for index,fwd in enumerate(self.ranges_fwd):
            rev = self.ranges_rev[index]
            self.amplicon.append([fwd[0],rev[1],fwd[3].split("_LEFT")[0]])

    def is_contained(self, end_pos, end_type):
        '''
        This function return range of primer region containing a given position.
        Returns None if there is no containing interval.
        '''
        if end_type == "left":
            prev = self.ranges_fwd[0]
            for range in self.ranges_fwd:
                if range[0] > end_pos:
                    return [None,prev[3]]
                if range[0] <= end_pos and end_pos < range[1]:
                    return [range, range[3]]
                prev=range
        if end_type == "right":
            for range in self.ranges_rev:
                if range[0] > end_pos:
                    return [None,range[3]]
                if range[0] < end_pos and end_pos <= range[1]:
                    return [range,range[3]]
        return [None,None]

    def exact_match(self, fragment, name_l,name_r):
        seq_l = self.primer_to_seq[name_l]
        match_l = fragment.alignments[0].read_seq.startswith(seq_l)
        
        seq_r = self.primer_to_seq[name_r]
        match_r = fragment.alignments[1].read_seq.endswith(seq_r)

        return match_l&match_r



