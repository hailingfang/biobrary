import re
import gzip


class GTF_BASE:
    def __init__(self):
        self.__seqname = None
        self.__source = None
        self.__feature = None
        self.__ori = None
        self.__score = None
        self.__frame = None
        self.__attr = None
        self.__range = None

    def get_seqname(self):
        return self.__seqname
    
    def get_source(self):
        return self.__source
    
    def get_feature(self):
        return self.__feature
    
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


class GTF_STOP_CODON(GTF_BASE):
    def __init__(self):
        super().__init__()
        self._GTF_BASE__feature = 'start_codon'

    def init_by_gtf_line(self, data_line):
        seqid, source, frange, score, ori, frame, attr_dic = data_line
        self._GTF_BASE__seqname = seqid
        self._GTF_BASE__source = source
        self._GTF_BASE__score = score
        self._GTF_BASE__ori = ori
        self._GTF_BASE__frame = frame
        self._GTF_BASE__attr = attr_dic
        self._GTF_BASE__range = frange
        assert "gene_id" in self._GTF_BASE__attr
        assert "transcript_id" in self._GTF_BASE__attr
        assert "protein_id" in self._GTF_BASE__attr
        self.__gene_id = attr_dic["gene_id"]
        self.__transcript_id = attr_dic["transcript_id"]
        self.__protein_id = attr_dic["protein_id"]

    def get_gene_id(self):
        return self.__gene_id

    def get_transcript_id(self):
        return self.__transcript_id
    
    def get_protein_id(self):
        return self.__protein_id


class GTF_START_CODON(GTF_BASE):
    def __init__(self):
        super().__init__()
        self._GTF_BASE__feature = 'stop_codon'

    def init_by_gtf_line(self, data_line):
        seqid, source, frange, score, ori, frame, attr_dic = data_line
        self._GTF_BASE__seqname = seqid
        self._GTF_BASE__source = source
        self._GTF_BASE__score = score
        self._GTF_BASE__ori = ori
        self._GTF_BASE__frame = frame
        self._GTF_BASE__attr = attr_dic
        self._GTF_BASE__range = frange
        assert "gene_id" in self._GTF_BASE__attr
        assert "transcript_id" in self._GTF_BASE__attr
        assert "protein_id" in self._GTF_BASE__attr
        self.__gene_id = attr_dic["gene_id"]
        self.__transcript_id = attr_dic["transcript_id"]
        self.__protein_id = attr_dic["protein_id"]

    def get_gene_id(self):
        return self.__gene_id

    def get_transcript_id(self):
        return self.__transcript_id
    
    def get_protein_id(self):
        return self.__protein_id


class GTF_CDS(GTF_BASE):
    def __init__(self):
        super().__init__()
        self._GTF_BASE__feature = 'CDS'

    def init_by_gtf_line(self, data_line):
        seqid, source, frange, score, ori, frame, attr_dic = data_line
        self._GTF_BASE__seqname = seqid
        self._GTF_BASE__source = source
        self._GTF_BASE__score = score
        self._GTF_BASE__ori = ori
        self._GTF_BASE__frame = frame
        self._GTF_BASE__attr = attr_dic
        self._GTF_BASE__range = frange
        assert "gene_id" in self._GTF_BASE__attr
        assert "transcript_id" in self._GTF_BASE__attr
        assert "protein_id" in self._GTF_BASE__attr
        self.__gene_id = attr_dic["gene_id"]
        self.__transcript_id = attr_dic["transcript_id"]
        self.__protein_id = attr_dic["protein_id"]

    def get_gene_id(self):
        return self.__gene_id

    def get_transcript_id(self):
        return self.__transcript_id
    
    def get_protein_id(self):
        return self.__protein_id


class GTF_EXON(GTF_BASE):
    def __init__(self):
        super().__init__()
        self._GTF_BASE__feature = 'exon'

    def init_by_gtf_line(self, data_line):
        seqid, source, frange, score, ori, frame, attr_dic = data_line
        self._GTF_BASE__seqname = seqid
        self._GTF_BASE__source = source
        self._GTF_BASE__score = score
        self._GTF_BASE__ori = ori
        self._GTF_BASE__frame = frame
        self._GTF_BASE__attr = attr_dic
        self._GTF_BASE__range = frange
        assert "gene_id" in self._GTF_BASE__attr
        assert "transcript_id" in self._GTF_BASE__attr
        self.__gene_id = attr_dic["gene_id"]
        self.__transcript_id = attr_dic["transcript_id"]
    
    def get_gene_id(self):
        return self.__gene_id

    def get_transcript_id(self):
        return self.__transcript_id


class GTF_TRANSCRIPT(GTF_BASE):
    def __init__(self):
        super().__init__()
        self._GTF_BASE__feature = 'transcript'

    def init_by_gtf_line(self, data_line):
        seqid, source, left, right, score, ori, frame, attr_dic = data_line
        self._GTF_BASE__seqname = seqid
        self._GTF_BASE__source = source
        self._GTF_BASE__score = score
        self._GTF_BASE__ori = ori
        self._GTF_BASE__frame = [frame]
        self._GTF_BASE__attr = attr_dic
        self._GTF_BASE__range = [(left, right)]
        assert "gene_id" in self._GTF_BASE__attr
        assert "transcript_id" in self._GTF_BASE__attr
        self.__gene_id = attr_dic["gene_id"]
        self.__transcript_id = attr_dic["transcript_id"]
        self.__biotype = None
        self.__exon = None
        self.__CDS = None
        self.__start_codon = None
        self.__stop_codon = None

    def get_gene_id(self):
        return self.__gene_id

    def get_transcript_id(self):
        return self.__transcript_id

    def get_biotype(self):
        return self.__biotype

    def get_exon(self):
        return self.__exon

    def get_CDS(self):
        return self.__CDS

    def get_start_codon(self):
        return self.__start_codon

    def get_stop_codon(self):
        return self.__stop_codon


class GTF_GENE(GTF_BASE):
    def __init__(self):
        super().__init__()
        self._GTF_BASE__feature = 'gene'

    def init_by_gtf_line(self, data_line):
        seqid, source, left, right, score, ori, frame, attr_dic = data_line
        self._GTF_BASE__seqname = seqid
        self._GTF_BASE__source = source
        self._GTF_BASE__score = score
        self._GTF_BASE__ori = ori
        self._GTF_BASE__frame = [frame]
        self._GTF_BASE__attr = attr_dic
        self._GTF_BASE__range = [(left, right)]
        assert "gene" in attr_dic
        assert "gene_id" in attr_dic
        self.__gene = attr_dic["gene"]
        self.__gene_id = attr_dic["gene_id"]
        self.__transcript_id_transcript_dic = {}

    def get_gene_name(self):
        return self.__gene

    def get_gene_id(self):
        return self.__gene_id

    def get_transcript_id_s(self):
        return list(self.__transcript_id_transcript_dic.keys())

    def get_transcript_s(self):
        return list(self.__transcript_id_transcript_dic.values())

    def get_transcript(self, transcript_id):
        return self.__transcript_id_transcript_dic.get(transcript_id)


class GTF:
    def __init__(self):
        self.__meta = None
        self.__gene_id_gene_dic = {}
        self.__gene_name_gene_id_dic = {}

    def __iter__(self):
        self.__start = 0
        self.__stop = len(self.__gene_s)
        self.__gene_id_s = list(self.__gene_id_gene_dic.keys())
        return self

    def __next__(self):
        if self.__start == self.__stop:
            raise StopIteration
        self.__start += 1
        return self.__gene_id_gene_dic[self.__gene_id_s[self.__start - 1]]

    def get_gene_ids(self):
        return list(self.__gene_id_gene_dic.keys())

    def get_gene_s(self):
        return list(self.__gene_id_gene_dic.values())

    def get_gene(self, gene_id):
        return self.__gene_id_gene_dic.get(gene_id)

    def get_meta(self):
        return self.__meta


def __read_gtf(gtf_file):
    meta = []
    data = {}
    retmp = re.compile(r'\s?(.+?)\s"(.+?)";')
    
    #read file
    if gtf_file.endswith(".gz"):
        fin = gzip.open(gtf_file, "r")
        for line in fin:
            line = line.decode().rstrip()
            if line[0] == "#":
                meta.append(line)
                continue
            
            line = line.split("\t")
            seqid, source, feature, left, right, score, ori, frame, attr = line
            attr = retmp.findall(attr)
            attr_dic = {}
            for ele in attr:
                if ele[0] not in attr_dic:
                    attr_dic[ele[0]] = ele[1]
                else:
                    pre_v = attr_dic[ele[0]]
                    attr_dic[ele[0]] = [pre_v, ele[1]]
            if feature in data:
                data[feature].append([seqid, source, int(left), int(right), score, ori, frame, attr_dic])
            else:
                data[feature] = [[seqid, source, int(left), int(right), score, ori, frame, attr_dic]] 
        fin.close()
    else:
        fin = open(gtf_file, "r")
        for line in fin:
            line = line.rstrip()
            if line[0] == "#":
                meta.append(line)
                continue

            line = line.split("\t")
            seqid, source, feature, left, right, score, ori, frame, attr = line
            attr = retmp.findall(attr)
            attr_dic = {}
            for ele in attr:
                if ele[0] not in attr_dic:
                    attr_dic[ele[0]] = ele[1]
                else:
                    pre_v = attr_dic[ele[0]]
                    attr_dic[ele[0]] = [pre_v, ele[1]]
            if feature in data:
                data[feature].append([seqid, source, int(left), int(right), score, ori, frame, attr_dic])
            else:
                data[feature] = [[seqid, source, int(left), int(right), score, ori, frame, attr_dic]] 
        fin.close()

    return meta, data


def __meger_feature(data, feature, attr_key):
    if feature not in data:
        return

    data_fea = data[feature]
    group_dic = {}
    ignored_line_num = 0
    for line in data_fea:
        key = line[-1].get(attr_key)
        if key:
            if key not in group_dic:
                group_dic[key] = [line]
            else:
                group_dic[key].append(line)
        else:
            ignored_line_num += 1
            #print(f"{attr_key} not found in {line}")
    
    if ignored_line_num:
        print(f"{ignored_line_num} entries in {feature} were ommited, because less {attr_key} attribute.")

    grouped_data = []
    for key in group_dic:
        seqid_g = []
        source_g = []
        left_g = []
        right_g = []
        score_g = []
        ori_g = []
        frame_g = []
        attr_g = []
        for line in group_dic[key]:
            seqid, source, left, right, score, ori, frame, attr = line
            seqid_g.append(seqid)
            source_g.append(source)
            left_g.append(left)
            right_g.append(right)
            score_g.append(score)
            ori_g.append(ori)
            frame_g.append(frame)
            attr_g.append(attr)
        seqid = seqid_g[0]
        source = source_g[0]
        frange = [ele for ele in zip(left_g, right_g)]
        score = score_g[0]
        ori = ori_g[0]
        frame = frame_g
        attr = attr_g[0]
        grouped_data.append([seqid, source, frange, score, ori, frame, attr])
    data[feature] = grouped_data


def parse_gtf(gtf_file):
    meta, gtf_data = __read_gtf(gtf_file)
    #print(list(gtf_data.keys()))

    __meger_feature(gtf_data, "exon", "transcript_id")
    __meger_feature(gtf_data, "CDS", "protein_id")
    __meger_feature(gtf_data, "start_codon", "protein_id")
    __meger_feature(gtf_data, "stop_codon", "protein_id")

    gtf = GTF()
    gtf._GTF__meta = meta
    print(gtf.get_meta())

    assert 'gene' in gtf_data
    gene_lines = gtf_data['gene']
    gene_id_gene_dic = {}
    gene_name_gene_id_dic = {}
    for line in gene_lines:
        gene = GTF_GENE()
        gene.init_by_gtf_line(line)
        gene_name = gene.get_gene_name()
        gene_id = gene.get_gene_id()
        if gene_name not in gene_name_gene_id_dic:
            gene_name_gene_id_dic[gene_name] = [gene_id]
        else:
            gene_name_gene_id_dic[gene_name].append(gene_id)
        gene_id_gene_dic[gene_id] = gene
    gtf._GTF__gene_id_gene_dic = gene_id_gene_dic
    gtf._GTF__gene_name_gene_id_dic = gene_name_gene_id_dic

    if 'transcript' in gtf_data:
        transcript_lines = gtf_data['transcript']
        for line in transcript_lines:
            transcript = GTF_TRANSCRIPT()
            transcript.init_by_gtf_line(line)
            gene_id = transcript.get_gene_id()
            transcript_id = transcript.get_transcript_id()
            gene = gene_id_gene_dic[gene_id]
            gene._GTF_GENE__transcript_id_transcript_dic[transcript_id] = transcript

    if 'exon' in gtf_data:
        exon_lines = gtf_data['exon']
        for line in exon_lines:
            exon = GTF_EXON()
            exon.init_by_gtf_line(line)
            gene_id = exon.get_gene_id()
            transcript_id = exon.get_transcript_id()
            gene = gene_id_gene_dic[gene_id]
            if transcript_id in gene.get_transcript_id_s():
                transcript = gene.get_transcript(transcript_id)
                transcript._GTF_TRNASCRIPT__exon = exon
            else:
                transcript = GTF_TRANSCRIPT()
                transcript._GTF_TRANSCRIPT__gene_id = gene_id
                transcript._GTF_TRANSCRIPT__transcript_id = transcript_id
                transcript._GTF_TRANSCRIPT__exon = exon
                gene._GTF_GENE__transcript_id_transcript_dic[transcript_id] = transcript
            
    if 'CDS' in gtf_data:
        cds_lines = gtf_data['CDS']
        for line in cds_lines:
            cds = GTF_CDS()
            cds.init_by_gtf_line(line)
            gene_id = cds.get_gene_id()
            transcript_id = cds.get_transcript_id()
            gene = gene_id_gene_dic[gene_id]
            if transcript_id in gene.get_transcript_id_s():
                transcript = gene.get_transcript(transcript_id)
                transcript._GTF_TRANSCRIPT_biotype = 'mRNA'
                transcript._GTF_TRANSCRIPT__CDS = cds
            else:
                transcript = GTF_TRANSCRIPT()
                transcript._GTF_TRANSCRIPT__gene_id = gene_id
                transcript._GTF_TRANSCRIPT__transcript_id = transcript_id
                transcript._GTF_TRANSCRIPT__CDS = cds
                gene._GTF_GENE__transcript_id_transcript_dic[transcript_id] = transcript

    if 'start_codon' in gtf_data:
        start_codon_lines = gtf_data['start_codon']
        for line in start_codon_lines:
            start_codon = GTF_START_CODON()
            start_codon.init_by_gtf_line(line)
            gene_id = start_codon.get_gene_id()
            transcript_id = start_codon.get_transcript_id()
            gene = gene_id_gene_dic[gene_id]
            if transcript_id in gene.get_transcript_id_s():
                transcript = gene.get_transcript(transcript_id)
                transcript._GTF_TRANSCRIPT__start_codon = start_codon
            else:
                transcript = GTF_TRANSCRIPT()
                transcript._GTF_TRANSCRIPT__gene_id = gene_id
                transcript._GTF_TRANSCRIPT__transcript_id = transcript_id
                transcript._GTF_TRANSCRIPT__start_codon = start_codon
                gene._GTF_GENE__transcript_id_transcript_dic[transcript_id] = transcript

    if 'stop_codon' in gtf_data:
        stop_codon_lines = gtf_data['stop_codon']
        for line in stop_codon_lines:
            stop_codon = GTF_START_CODON()
            stop_codon.init_by_gtf_line(line)
            gene_id = stop_codon.get_gene_id()
            transcript_id = stop_codon.get_transcript_id()
            gene = gene_id_gene_dic[gene_id]
            if transcript_id in gene.get_transcript_id_s():
                transcript = gene.get_transcript(transcript_id)
                transcript._GTF_TRANSCRIPT__stop_codon = stop_codon
            else:
                transcript = GTF_TRANSCRIPT()
                transcript._GTF_TRANSCRIPT__gene_id = gene_id
                transcript._GTF_TRANSCRIPT__transcript_id = transcript_id
                transcript._GTF_TRANSCRIPT__stop_codon = stop_codon
                gene._GTF_GENE__transcript_id_transcript_dic[transcript_id] = transcript
    return gtf


def test_gtf(gtf_file):
    parse_gtf(gtf_file)


if __name__ == "__main__":
    import sys
    test_gtf(sys.argv[1])
