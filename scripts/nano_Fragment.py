import Alignment

class Fragment():
    '''
    This object represents a molecule inserted between the adapters.
    '''
    def __init__(self, alignment):

        self.alignment = alignment
        self.read_name = alignment.read_name
        self.ref_start = alignment.ref_start
        self.ref_end = alignment.ref_end
        self.size = self.ref_end - self.ref_start
        self.read_seq = alignment.read_seq
        self.read_q = alignment.read_q

    def get_fastqlines(self):
        '''
        Returns a list of sequences and
        qualities of both reads.
        '''
        return ["@%s"%self.read_name,self.read_seq,"+",self.read_q]

    def slice(self, amp_def):
        '''
        Slices the fragment at the pos.
        '''
        def get_range(amp_def,alignment):
            start,end = amp_def
            range_slice = [len(alignment.read2refcorr),0]
            seq_index = [1e6,0]
            primer_LEFT, primer_RIGTH, ratio = 0,0,0
            for index,pos in enumerate(alignment.read2refcorr):
                if pos=="S":
                    continue
                if (pos<start):
                    primer_LEFT+=1
                    continue
                if (pos>end):
                    primer_RIGTH+=1
                    continue
                if pos<=seq_index[0]:
                    seq_index[0] = pos
                    range_slice[0] = index
                if pos>=seq_index[1]:
                    seq_index[1] = pos
                    range_slice[1] = index+1 # so that last pos in the range is included, not excluded as in range(0,10) would not output 10.
            if primer_RIGTH*primer_LEFT:
                ratio = (range_slice[1]-range_slice[0])/float(end-start)
            return range_slice,seq_index,ratio

        range_slice,_,ratio = get_range(amp_def,self.alignment)
        start,end = range_slice

        self.ratio = ratio
        self.alignment.read2refcorr = self.alignment.read2refcorr[start:end]
        self.alignment.read_q = self.alignment.read_q[start:end]
        self.alignment.read_seq = self.alignment.read_seq[start:end]

