import sys

import numpy as np
from Bio import SeqIO

import RNA


def read_alignment(fpath, ali_format='fasta'):
    with open(fpath) as f:
        for record in SeqIO.parse(f, ali_format):
            seq = str(record.seq)
            seq = seq.upper().replace("T", "U")
            yield seq


class Alignment:
    def __init__(self, alignment, md=RNA.md()):
        """Create alignment object from a given aligned file
        """
        self.md = md
        self.aligned = [seq for seq in alignment]
        self.ungapped = [seq.replace('-','') for seq in self.aligned]
        self.consensus = ""
        self.fc = None
        self.ps_energy = []      # 1-indexed
        self.ps_energy_up = []   # 1-indexed
        self.dominant_bps = []   # 1-indexed
        
        self._fold()

    @classmethod
    def from_file(cls, filepath, md=RNA.md(), ali_format='fasta'):
        return cls(read_alignment(fpath, ali_format), md)



    def _fold(self):
        """Get consensus structure and bpp
        """
        fc = RNA.fold_compound(self.aligned, self.md)
        ss, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        _, self.fee = fc.pf()

        self.consensus = ss
        self.fc = fc

        self.ps_energy, self.ps_energy_up = self._pseudo_energy()

        self.dominant_bps = [(ai, aj) for ai, aj in zip(*np.where(self.ps_energy<0))]


    def _pseudo_energy(self):
        """Create pseudo energy matrix
        """
        RT = RNA.GASCONST*(self.md.temperature + RNA.K0)/1000
        bpp = np.array(self.fc.bpp())
        bpp_pp = bpp.copy()
        np.set_printoptions(threshold=np.inf)
        paired = bpp>=0.99
        unpaired = bpp==0
        bpp[paired] = np.NaN
        bpp[unpaired] = np.NaN
        tmp = np.minimum(0, -RT*np.log(bpp/(1-bpp)))

        pos_paired = np.argwhere(paired==True)
        
        fc = RNA.fold_compound(self.aligned)
        ss, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)

        for bp in pos_paired:
            fc.hc_init()
            #this really enforces the base pair
            fc.hc_add_bp(int(bp[0]), int(bp[1]), RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE)
            _, fee_w = fc.pf()

            fc.hc_init()
            #this one forbids the base pair but accepts all other base pairs that include i and j
            fc.hc_add_bp(int(bp[0]), int(bp[1]), RNA.CONSTRAINT_CONTEXT_NONE)
            _, fee_wo = fc.pf()

            tmp[bp[0]][bp[1]] = fee_w - fee_wo

        tmp[unpaired] = 0

        bpp_pp += np.transpose(bpp_pp)

        paired_prob = np.sum(bpp_pp, axis=0)

        paired_pp = paired_prob >= 0.99
        unpaired_pp = paired_prob == 0
        paired_prob[paired_pp] = np.NaN
        paired_prob[unpaired_pp] = np.NaN

        pos_paired_pp = np.argwhere(paired_pp==True)

        tmp_pp = np.maximum(0, RT*np.log(paired_prob/(1-paired_prob)))

        for pos in pos_paired_pp:
            fc.hc_init()
            fc.hc_add_bp_nonspecific(int(pos), 0, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE)
            _, fee_p = fc.pf()

            fc.hc_init()
            fc.hc_add_up(int(pos), RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE)
            _, fee_up = fc.pf()

            tmp_pp[pos] = -(fee_p - fee_up)

        tmp_pp[unpaired_pp] = 0
        return tmp, tmp_pp


    def get_energy_of_aligned_seq(self, aligned, up_e):
        """Return pseudo energy matrix for aligned sequence
        """

        ungap = np.array(['X']+list(aligned)) != '-'
        l = np.sum(ungap)
        ungap_bps = np.outer(ungap, ungap)
        if up_e:
            res = self.ps_energy_up.copy()
            return res[ungap].reshape(l)
        else:
            res = self.ps_energy.copy()
            return res[ungap_bps].reshape((l,l))


    def add_sc(self, fc, aligned, up_e):
        """Add soft constraint to fc
        """
        to_add = self.get_energy_of_aligned_seq(aligned, up_e)
        if up_e:
            fc.sc_add_up(to_add.tolist())
        else:
            fc.sc_add_bp(to_add.tolist())


    def fc_of_aligned_seq(self, aligned, up_e=False):
        """Create fold compound for given aligned sequence with soft constraint
        """
        seq = aligned.replace('-', '')
        fc = RNA.fold_compound(seq, self.md)
        self.add_sc(fc, aligned, up_e)
        return fc
