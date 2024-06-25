import sys
from .data.base_complement import base_complement as basecom
from .data.biocodon import CODON_AA as condon_dic
from .data.molecular_weight import amino_acids_mw

class Seq:
    def __init__(self, seq, seqtype="nucl"):
        self.seq = seq
        self.seqtype = seqtype


    def __str__(self) -> str:
        return self.seq


    def complement(self):
        if self.seqtype == "nucl":
            return Seq("".join(basecom[base] for base in self.seq))
        else:
            print("Warning, the method only available for nucltides", file=sys.stderr)


    def reverse(self):
        seq = list(self.seq)
        seq.reverse()
        return Seq("".join(seq))


    def reverse_comp(self):
        if self.seqtype == "nucl":
            seq = self.seq
            seq = [basecom[base] for base in seq]
            seq.reverse()
            return Seq("".join(seq))
        else:
            print("Warning, the method only available for nucltides", file=sys.stderr)


    def translate(self):
        if self.seqtype == "nucl":
            seq = self.seq
            seq = [seq[idx: idx + 3] for idx in range(len(seq))[::3]]
            prot = ""
            for codon in seq:
                if len(codon) != 3:
                    break
                p = condon_dic[codon.upper()]
                prot += p
                if p == "*":
                    break
            return Seq(prot, seqtype="prot")
        else:
            print("Warning, the method only available for nucltides", file=sys.stderr)
        

    def molecular_weight(self):
        if self.seqtype == "nucl":
            return 0
        elif self.seqtype == "prot":
            weight = 0
            for p in self.seq:
                weight += amino_acids_mw[p.upper()]
            weight -= (len(self.seq) - 1) * 18
            return weight
        else:
            return 0
