import re


class GTF_BASE:
    def __init__(self):
        self.__seqname = None
        self.__source = None
        self.__feature = None
        self.__left = None
        self.__right = None
        self.__ori = None
        self.__score = None
        self.__frame = None
        self.__attr = {}
        self.__range = []

    def get_seqname(self):
        return self.__seqname
    
    def get_source(self):
        return self.__source
    
    def get_feature(self):
        return self.__feature
    
    def get_left(self):
        return self.__left
    
    def get_right(self):
        return self.__right
    
    def get_ori(self):
        return self.__ori
    
    def get_score(self):
        return self.__score
    
    def get_frame(self):
        return self.__frame
    
    def get_attr(self, key):
        return self.__attr.get(key)

    def get_attr_dic(self):
        return self.__attr

    def get_range(self):
        return self.__range


class GTF_START_STOP_CODON(GTF_BASE):
    def __init__(self, data_in): 
        super().__init__()

        seqname, source, feature, left, right, score, ori, frame, attrib \
                = data_in[0]
        self._GTF_BASE__seqname = seqname
        self._GTF_BASE__source = source
        self._GTF_BASE__feature = feature
        self._GTF_BASE__ori = ori
        self._GTF_BASE__score = score
        self._GTF_BASE__attr = attrib
        if "exon_number" in self._GTF_BASE__attr:
            self._GTF_BASE__attr.pop("exon_number")
        left = []
        right = []
        frame = []
        for line in data_in:
            left.append(int(line[3]))
            right.append(int(line[4]))
            frame.append(line[7])
        self._GTF_BASE__range = ([ele for ele in zip(left, right)])
        self._GTF_BASE__range.sort(key=lambda x:x[0])
        self._GTF_BASE__left = self._GTF_BASE__range[0][0]
        self._GTF_BASE__right = self._GTF_BASE__range[-1][1]
        self._GTF_BASE__frame = frame


class GTF_PROTEIN(GTF_BASE):
    def __init__(self, data_in):
        super().__init__()
        
        self.__start_stop_codon_ids = []
        self.__child = []

        cds_lines = []
        start_codon_lines = []
        stop_codon_lines = []
        for line in data_in:
            if line[2] == "CDS":
                cds_lines.append(line)
            elif line[2] == "start_codon":
                start_codon_lines.append(line)
            elif line[2] == "stop_codon":
                stop_codon_lines.append(line)
            else:
                print(f"feature is not known, {line}", file=sys.stderr)
        
        if cds_lines:
            seqname, source, feature, left, right, score, ori, frame, attrib \
                = cds_lines[0]
            self._GTF_BASE__seqname = seqname
            self._GTF_BASE__source = source
            self._GTF_BASE__feature = "protein"
            self._GTF_BASE__score = score
            self._GTF_BASE__ori = ori
            self._GTF_BASE__attr = attrib
            if "exon_number" in self._GTF_BASE__attr:
                self._GTF_BASE__attr.pop("exon_number")

            left = []
            right = []
            frame = []
            for line in cds_lines:
                left.append(int(line[3]))
                right.append(int(line[4]))
                frame.append(line[7])
            self._GTF_BASE__range = [ele for ele in zip(left, right)]
            self._GTF_BASE__range.sort(key=lambda x:x[0])
            self._GTF_BASE__left = self._GTF_BASE__range[0][0]
            self._GTF_BASE__right = self._GTF_BASE__range[-1][1]
            self._GTF_BASE__frame = frame

        elif start_codon_lines:
            seqname, source, feature, left, right, score, ori, frame, attrib \
                = start_codon_lines[0]
            self._GTF_BASE__seqname = seqname
            self._GTF_BASE__source = source
            self._GTF_BASE__feature = "protein"
            self._GTF_BASE__score = score
            self._GTF_BASE__ori = ori
            self._GTF_BASE__attr = attrib
            if "exon_number" in self._GTF_BASE__attr:
                self._GTF_BASE__attr.pop("exon_number")
        elif stop_codon_lines:
            seqname, source, feature, left, right, score, ori, frame, attrib \
                = stop_codon_lines[0]
            self._GTF_BASE__seqname = seqname
            self._GTF_BASE__source = source
            self._GTF_BASE__feature = "protein"
            self._GTF_BASE__score = score
            self._GTF_BASE__ori = ori
            self._GTF_BASE__attr = attrib
            if "exon_number" in self._GTF_BASE__attr:
                self._GTF_BASE__attr.pop("exon_number")

        if len(start_codon_lines) > 0:
            self.__start_stop_codon_ids.append("start_codon")
            self.__child.append(GTF_START_STOP_CODON(start_codon_lines))
        
        if len(stop_codon_lines) > 0:
            self.__start_stop_codon_ids.append("stop_codon")
            self.__child.append(GTF_START_STOP_CODON(stop_codon_lines))


    def get_protein_id(self):
        return self._GTF_BASE__attr["protein_id"]


    def get_start_stop_ids(self):
        return self.__start_stop_codon_ids


    def get_start_stop(self, start_stop_name):
        return self.__child[self.__start_stop_codon_ids.index(start_stop_name)]


class GTF_TRANSCRIPT(GTF_BASE):
    def __init__(self, data_in):
        super().__init__()

        self.__protein_ids = []
        self.__child = []

        transcript_line = []
        exon_lines = []
        cds_lines = []
        start_codon_lines = []
        stop_codon_lines = []
        for line in data_in:
            if line[2] == "exon":
                exon_lines.append(line)
            elif line[2] == "CDS":
                cds_lines.append(line)
            elif line[2] == "start_codon":
                start_codon_lines.append(line)
            elif line[2] == "stop_codon":
                stop_codon_lines.append(line)
            elif line[2] == "transcript":
                transcript_line.append(line)
            else:
                print(f"feature is not known, {line}", file=sys.stderr)

        if transcript_line:
            seqname, source, feature, left, right, score, ori, frame, attrib \
                = transcript_line[0]
            self._GTF_BASE__seqname = seqname
            self._GTF_BASE__source = source
            self._GTF_BASE__feature = feature
            self._GTF_BASE__left = left
            self._GTF_BASE__right = right
            self._GTF_BASE__score = score
            self._GTF_BASE__ori = ori
            self._GTF_BASE__frame = frame
            self._GTF_BASE__attr = attrib
            if exon_lines:
                left = []
                right = []
                for line in exon_lines:
                    left.append(int(line[3]))
                    right.append(int(line[4]))
                    self._GTF_BASE__range = ([ele for ele in zip(left, right)])
                    self._GTF_BASE__range.sort(key=lambda x:x[0])
            elif cds_lines:
                left = []
                right = []
                for line in cds_lines:
                    left.append(int(line[3]))
                    right.append(int(line[4]))
                    self._GTF_BASE__range([ele for ele in zip(left, right)])
                    self._GTF_BASE__range.sort(key=lambda x:x[0])
        elif exon_lines:
            seqname, source, feature, left, right, score, ori, frame, attrib \
                = exon_lines[0]
            self._GTF_BASE__seqname = seqname
            self._GTF_BASE__source = source
            self._GTF_BASE__feature = "transcript"
            self._GTF_BASE__score = score
            self._GTF_BASE__ori = ori
            self._GTF_BASE__frame = frame
            self._GTF_BASE__attr = attrib
            if "exon_number" in self._GTF_BASE__attr:
                self._GTF_BASE__attr.pop("exon_number")
            left = []
            right = []
            for line in exon_lines:
                left.append(int(line[3]))
                right.append(int(line[4]))
            self._GTF_BASE__range = [ele for ele in zip(left, right)]
            self._GTF_BASE__range.sort(key=lambda x:x[0])
            self._GTF_BASE__left = self._GTF_BASE__range[0][0]
            self._GTF_BASE__right = self._GTF_BASE__range[-1][1]

        elif cds_lines:
            seqname, source, feature, left, right, score, ori, frame, attrib \
                = cds_lines[0]
            if "gbkey" in attrib:
                attrib["gbkey"] = "mRNA"
            self._GTF_BASE__seqname = seqname
            self._GTF_BASE__source = source
            self._GTF_BASE__feature = "transcript"
            self._GTF_BASE__score = score
            self._GTF_BASE__ori = ori
            self._GTF_BASE__frame = frame            
            self._GTF_BASE__attr = attrib
            if "exon_number" in self._GTF_BASE__attr:
                self._GTF_BASE__attr.pop("exon_number")
            left = []
            right = []
            for line in cds_lines:
                left.append(int(line[3]))
                right.append(int(line[4]))
            self._GTF_BASE__range = [ele for ele in zip(left, right)]
            self._GTF_BASE__range.sort(key=lambda x:x[0])
            self._GTF_BASE__left = self._GTF_BASE__range[0][0]
            self._GTF_BASE__right = self._GTF_BASE__range[-1][1]
        else:
            print("error", file=sys.stderr)

        cds_lines = cds_lines + start_codon_lines + stop_codon_lines
        protein_data = {}
        for line in cds_lines:
            protein_id = line[-1].get("protein_id")
            if protein_id:
                if protein_id in protein_data:
                    protein_data[protein_id].append(line)
                else:
                    protein_data[protein_id] = [line]
                    self.__protein_ids.append(protein_id)
        assert len(self.__protein_ids) == len(set(self.__protein_ids))
        for protein_id in protein_data:
            self.__child.append(GTF_PROTEIN(protein_data[protein_id]))



    def get_trancript_id(self):
        return self._GTF_BASE__attr["transcript_id"]


    def get_protein_ids(self):
        return self.__protein_ids

    
    def get_protein(self, protein_id):
        if protein_id in self.__protein_ids:
            return self.__child[self.__protein_ids.index(protein_id)]
        else:
            return None


class GTF_GENE(GTF_BASE):
    def __init__(self, gene_data):
        super().__init__()
        self.__transcript_ids = []
        self.__child = []

        gene_line = gene_data[0]
        seqname, source, feature, left, right, score, ori, frame, attrib_dic \
            = gene_line
        assert feature == "gene"
        
        self._GTF_BASE__seqname = seqname
        self._GTF_BASE__source = source
        self._GTF_BASE__feature = feature
        self._GTF_BASE__left = int(left)
        self._GTF_BASE__right = int(right)
        self._GTF_BASE__score = score
        self._GTF_BASE__ori = ori
        self._GTF_BASE__frame = frame
        self._GTF_BASE__attr = attrib_dic
        self._GTF_BASE__range = [(self._GTF_BASE__left, self._GTF_BASE__right)]
        assert "gene_id" in self._GTF_BASE__attr

        transcript_lines = {}
        for line in gene_data[1:]:
            transcript_id = line[-1]["transcript_id"]
            if transcript_id in transcript_lines:
                transcript_lines[transcript_id].append(line)
            else:
                transcript_lines[transcript_id] = [line]
                self.__transcript_ids.append(transcript_id)
        assert len(self.__transcript_ids) == len(set(self.__transcript_ids))

        for transcript_id in transcript_lines:
            self.__child.append(GTF_TRANSCRIPT(transcript_lines[transcript_id]))


    def get_gene_id(self):
        return self._GTF_BASE__attr["gene_id"]


    def get_transcript_ids(self):
        return self.__transcript_ids


    def get_transcript(self, transcript_id):
        if transcript_id in self.__transcript_ids:
            return self.__child[self.__transcript_ids.index(transcript_id)]
        else:
            return None


class GTF:
    def __init__(self, gtf_file):
        self.__meta = None
        self.__gene_ids = []
        self.__child = []

        data = []
        meta = []
        fin = open(gtf_file, "r")
        for line in fin:
            line = line.rstrip()
            if line[0] != "#":
                data.append(line)
            else:
                meta.append(line)
        fin.close()

        self.__meta = meta

        retmp = re.compile(r'\s?(.+?)\s"(.+?)";')
        gtf_gene_data = {}
        gene_id = None
        for line in data:
            line = line.split("\t")
            feature = line[2]
            attrib = retmp.findall(line[-1])
            attrib_dic = {}
            unique_attr = ["gene", "gene_id", "transcript_id", "protein_id", "gbkey"]
            for ele in attrib:
                if ele[0] in unique_attr:
                    attrib_dic[ele[0]] = ele[1]
                elif ele[0] not in attrib_dic:
                    attrib_dic[ele[0]] = [ele[1]]
                else:
                    attrib_dic[ele[0]].append(ele[1])
            attrib = attrib_dic
            line[-1] = attrib
            if feature == "gene":
                assert "gene_id" in attrib
                gene_id = attrib["gene_id"]
                gtf_gene_data[gene_id] = [line]
                self.__gene_ids.append(gene_id)
            else:
                gtf_gene_data[gene_id].append(line)
        assert len(self.__gene_ids) == len(set(self.__gene_ids))
        for gene_id in gtf_gene_data:
            self.__child.append(GTF_GENE(gtf_gene_data[gene_id]))


    def __iter__(self):
        self.__start = 0
        self.__stop = len(self.__child)
        return self


    def __next__(self):
        if self.__start == self.__stop:
            raise StopIteration
        self.__start += 1
        return self.__child[self.__start - 1]


    def get_gene_ids(self):
        return self.__gene_ids


    def get_gene(self, gene_id):
        if gene_id in self.__gene_ids:
            return self.__child[self.__gene_ids.index(gene_id)]
        else:
            print(f"{gene_id} not found in GTF file.", file=sys.stderr)
            return None


    def get_meta(self):
        return self.__meta


if __name__ == "__main__":
    import sys
    import pdb
    breakpoint()
    gtf = GTF(sys.argv[1])
    for gene in gtf:
        print(gene.get_gene_id())
        for tran_id in gene.get_transcript_ids():
            print(tran_id)
            trans = gene.get_transcript(tran_id)
            for prot_id in trans.get_protein_ids():
                print(prot_id)
                prot = trans.get_protein(prot_id)
                print(prot.get_start_stop_ids())

