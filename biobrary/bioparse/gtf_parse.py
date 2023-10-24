class GTF_stop_codon:
    def __init__(self, stop_codon_data): 
        self.seqname = None
        self.source = None
        self.pos = None
        self.ori = None
        self.frame = None
        self.attr = None
        self.child = []
        seqname, source, feature, left, right, score, ori, frame, attrib_dic = \
            stop_codon_data[0]
        lefts = [int(left)]
        rights = [int(right)]
        for line in stop_codon_data[1:]:
            lefts.append(int(line[3]))
            rights.append(int(line[4]))
        self.seqname = seqname
        self.source = source
        self.pos = [[ele[0], ele[1]] for ele in zip(lefts, rights)]
        self.ori = ori
        self.frame = frame
        self.attr = attrib_dic
        


class GTF_start_codon:
    def __init__(self, start_codon_data):
        self.seqname = None
        self.source = None
        self.pos = None
        self.ori = None
        self.frame = None
        self.attr = None
        self.child = []
        seqname, source, feature, left, right, score, ori, frame, attrib_dic = \
            start_codon_data[0]
        lefts = [int(left)]
        rights = [int(right)]
        for line in start_codon_data[1:]:
            lefts.append(int(line[3]))
            rights.append(int(line[4]))
        self.seqname = seqname
        self.source = source
        self.pos = [[ele[0], ele[1]] for ele in zip(lefts, rights)]
        self.ori = ori
        self.frame = frame
        self.attr = attrib_dic


class GTF_CDS:
    def __init__(self, cds_data):
        self.seqname = None
        self.source = None
        self.pos = None
        self.ori = None
        self.frame = None
        self.attr = None

        seqname, source, feature, left, right, score, ori, frame, attrib_dic = \
            cds_data[0]

        lefts = [int(left)]
        rights = [int(right)]
        for line in cds_data[1:]:
            lefts.append(int(line[3]))
            rights.append(int(line[4]))
        self.seqname = seqname
        self.source = source
        self.pos = [[ele[0], ele[1]] for ele in zip(lefts, rights)]
        self.ori = ori
        self.frame = frame
        self.attr = attrib_dic


class GTF_exon:
    def __init__(self, exon_data):
        self.seqname = None
        self.source = None
        self.pos = None
        self.ori = None
        self.frame = None
        self.attr = None

        seqname, source, feature, left, right, score, ori, frame, attrib_dic = \
            exon_data[0]
        lefts = [int(left)]
        rights = [int(right)]
        for line in exon_data[1:]:
            lefts.append(int(line[3]))
            rights.append(int(line[4]))
        self.seqname = seqname
        self.source = source
        self.pos = [[ele[0], ele[1]] for ele in zip(lefts, rights)]
        self.ori = ori
        self.frame = frame
        self.attr = attrib_dic


class GTF_transcript:
    def __init__(self, data_in):
        self.seqname = None
        self.source = None
        self.pos = None
        self.ori = None
        self.frame = None
        self.attr = None
        self.child = []

        data_splited = {}
        for line in data_in:
            feature = line[2]
            if feature not in data_splited:
                data_splited[feature] = [line]
            else:
                data_splited[feature].append(line)


        child = []
        if "exon" in data_splited:
            child.append(GTF_exon(data_splited["exon"]))
        else:
            child.append(None)
        if "CDS" in data_splited:
            child.append(GTF_CDS(data_splited["CDS"]))
        else:
            child.append(None)
        if "start_codon" in data_splited:
            child.append(GTF_start_codon(data_splited["start_codon"]))
        else:
            child.append(None)
        if "stop_codon" in data_splited:
            child.append(GTF_stop_codon(data_splited["stop_codon"]))
        else:
            child.append(None)
        self.child = child

        if "transcript" in data_splited:
            assert len(data_splited["transcipt"]) == 1
            seqname, source, feature, left, right, score, ori, frame, attrib_dic\
                = data_splited["transcipt"][0]
            self.seqname = seqname
            self.source = source
            self.pos = [[int(left), int(right)]]
            self.ori = ori
            self.frame = frame
            self.attr = attrib_dic
            data_splited.pop("transcipt")
        else:
            assert "exon" in data_splited
            seqname, source, feature, left, right, score, ori, frame, attrib_dic\
                = data_splited["exon"][0]
            lefts = [int(left)]
            rights = [int(right)]
            for line in data_splited["exon"][1:]:
                lefts.append(int(line[3]))
                rights.append(int(line[4]))
            lefts.sort()
            rights.sort()
            self.seqname = seqname
            self.source = source
            self.pos = [[lefts[0], rights[-1]]]
            self.ori = ori
            self.frame = frame
            self.attr = attrib_dic
        


class GTF_gene:
    def __init__(self, gene_raw_dt):
        """
        Construct Gtf_gene object

        Parameters
        -------------
        gene_gtf_liens : lines related to one gene

        Returns
        ------------
        None

        self.data
        ---------------
        {transid: {exon: {exon_number}, CDS: [proteinid, {exon_number: {}}]}, ...}
        """

        self.seqname = None
        self.source = None
        self.pos = None
        self.ori = None
        self.frame = None
        self.attr = None
        self.child = []

        gene_line = gene_raw_dt[0]
        seqname, source, feature, left, right, score, ori, frame, attrib_dic = gene_line
        assert feature == "gene"
        self.seqname = seqname
        self.source = source
        self.pos = [[int(left), int(right)]]
        self.ori = ori
        self.attr = attrib_dic
        assert "gene_id" in self.attr


        gene_raw_dt = gene_raw_dt[1:]
        data = {}
        for line in gene_raw_dt:
            seqname, source, feature, left, right, score, ori, frame, attrib_dic = line
            transid = attrib_dic["transcript_id"][0]
            if transid not in data:
                data[transid] = [line]
            else:
                data[transid].append(line)
        
        for transid in data:
            self.child.append(GTF_transcript(data[transid]))
        
    
    def get_geneid(self):
        return self.attr["gene_id"]
    
    def get_gene_info(self):
        return [self.seqname, self.source, self.pos[0], self.ori]
    
    def get_gene_gbkey(self):
        if self.attr.get("gbkey"):
            return self.attr["gbkey"][0]
        else:
            return None

    def get_gene_biotype(self):
        if self.attr.get("gene_biotype"):
            return self.attr["gene_biotype"]
        else:
            return None


class GTF:
    def __init__(self, gtf_file):
        self.data = []
        self.meta = None
        self.geneids = None
        self.geneid_index = {}

        data = []
        fin = open(gtf_file, "r")
        for line in fin:
            data.append(line.rstrip())
        fin.close()

        meta = []
        meta_line_count = 0
        for line in data:
            if line[0] == "#":
                meta_line_count += 1
                meta.append(line)
            else:
                break
        self.meta = meta
        data = data[meta_line_count:]
        if data[-1][0] == "#":
            data = data[:-1]
        
        gtf_gene_data = {}
        index = 0
        key = None
        for line in data:
            line = line.split("\t")
            feature = line[2]
            attrib = line[-1]
            attrib = [ele.strip() for ele in attrib.split('";')[:-1]]
            attrib_dic = {}
            for ele in attrib:
                idx = ele.index(" ")
                key2 = ele[:idx]
                value = ele[idx+1:].strip('"')
                if len(value) > 0:
                    if value in attrib_dic:
                        attrib_dic[key2].append(value)
                    else:
                        attrib_dic[key2] = [value]
            line[-1] = attrib_dic

            if feature == "gene":
                index += 1
                key = index
                gtf_gene_data[key] = [line]
                continue
            gtf_gene_data[key].append(line)
        
        for key in gtf_gene_data:
            self.data.append(GTF_gene(gtf_gene_data[key]))
        
        index = 0
        for gene_data in self.data:
            geneid = gene_data.get_geneid()
            self.geneids.append(geneid)
            self.geneid_index[geneid] = index
            index += 1
        
        
    def __iter__(self):
        self.start = 0
        self.stop = len(self.data)
        return self

    def __next__(self):
        if self.start == self.stop:
            raise StopIteration
        self.start += 1
        return self.data[self.start - 1]


    def get_geneids(self):
        return self.geneids


    def get_gene(self, geneid):
        return self.data[self.geneid_index[geneid]]


if __name__ == "__main__":
    import sys
    gtf = GTF(sys.argv[1])
