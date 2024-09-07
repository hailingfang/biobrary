"""
Founctions and Class to parse GTF file
"""

import re
import gzip
import sys


class GTF_BASE:
    def __init__(self):
        self._seq_name = None
        self._source = None
        self._feature = None
        self._range = None
        self._ori = None
        self._score = None
        self._frame = None
        self._attr = None

    def get_seq_name(self):
        return self._seqname
    
    def get_source(self):
        return self._source
    
    def get_feature(self):
        return self._feature
    
    def get_range(self):
        return self._range

    def get_ori(self):
        return self._ori
    
    def get_score(self):
        return self._score
    
    def get_frame(self):
        return self._frame
    
    def get_attr(self, key):
        return self._attr.get(key)

    def get_attr_dic(self):
        return self._attr


class GET_GTF_ATTR:
    def get_gene_name(self):
        return self._gene_name

    def get_gene_id(self):
        return self._gene_id

    def get_biotype(self):
        return self._biotype

    def get_transcript_id(self):
        return self._transcript_id

    def get_protein_id(self):
        return self._protein_id


def _init_by_gtf_line_helper(self, seqid, source, frange, score, ori, frame, attr_dic):
        self._seqname = seqid
        self._source = source
        self._range = frange
        self._score = score
        self._ori = ori
        self._frame = frame
        self._attr = attr_dic
        self._gene_name = attr_dic.get('gene') if attr_dic.get('gene') else attr_dic.get('gene_name')
        assert "gene_id" in attr_dic
        self._gene_id = attr_dic["gene_id"]


class GTF_STOP_CODON(GTF_BASE):
    def __init__(self):
        super().__init__()
        self._feature = 'stop_codon'
        self._gene_name = None
        self._gene_id = None
        self._transcript_id = None
        self._protein_id = None
        self._biotype = None

    def _init_by_gtf_line(self, seqid, source, frange, score, ori, frame, attr_dic):
        _init_by_gtf_line_helper(self, seqid, source, frange, score, ori, frame, attr_dic)
        assert "transcript_id" in attr_dic
        assert "protein_id" in attr_dic
        self._transcript_id = attr_dic["transcript_id"]
        self._protein_id = attr_dic["protein_id"]


class GTF_START_CODON(GTF_BASE, GET_GTF_ATTR):
    def __init__(self):
        super().__init__()
        self._feature = 'start_codon'
        self._gene_name = None
        self._gene_id = None
        self._transcript_id = None
        self._protein_id = None
        self._biotype = None

    def _init_by_gtf_line(self, seqid, source, frange, score, ori, frame, attr_dic):
        _init_by_gtf_line_helper(self, seqid, source, frange, score, ori, frame, attr_dic)
        assert "transcript_id" in attr_dic
        assert "protein_id" in attr_dic
        self._transcript_id = attr_dic["transcript_id"]
        self._protein_id = attr_dic["protein_id"]


class GTF_CDS(GTF_BASE, GET_GTF_ATTR):
    def __init__(self):
        super().__init__()
        self._feature = 'CDS'
        self._gene_name = None
        self._gene_id = None
        self._transcript_id = None
        self._protein_id = None
        self._biotype = None

    def _init_by_gtf_line(self, seqid, source, frange, score, ori, frame, attr_dic):
        _init_by_gtf_line_helper(self, seqid, source, frange, score, ori, frame, attr_dic)
        assert "transcript_id" in attr_dic
        assert "protein_id" in attr_dic
        self._transcript_id = attr_dic["transcript_id"]
        self._protein_id = attr_dic["protein_id"]


class GTF_EXON(GTF_BASE, GET_GTF_ATTR):
    def __init__(self):
        super().__init__()
        self._feature = 'exon'
        self._gene_name = None
        self._gene_id = None
        self._transcript_id = None
        self._biotype = None

    def _init_by_gtf_line(self, seqid, source, frange, score, ori, frame, attr_dic):
        _init_by_gtf_line_helper(self, seqid, source, frange, score, ori, frame, attr_dic)
        assert "transcript_id" in attr_dic
        self._transcript_id = attr_dic["transcript_id"]
        self._biotype = attr_dic.get('transcript_biotype')


class GTF_TRANSCRIPT(GTF_BASE, GET_GTF_ATTR):
    def __init__(self):
        super().__init__()
        self._feature = 'transcript'
        self._gene_name = None
        self._gene_id = None
        self._transcript_id = None
        self._biotype = None
        self._exon = None
        self._CDS = None
        self._start_codon = None
        self._stop_codon = None

    def _init_by_gtf_line(self, seqid, source, left, right, score, ori, frame, attr_dic):
        frange = [(left, right)]
        _init_by_gtf_line_helper(self, seqid, source, frange, score, ori, frame, attr_dic)
        assert "transcript_id" in attr_dic
        self._transcript_id = attr_dic["transcript_id"]
        self._biotype = attr_dic.get('transcript_biotype')

    def get_exon(self):
        return self._exon

    def get_CDS(self):
        return self._CDS

    def get_start_codon(self):
        return self._start_codon

    def get_stop_codon(self):
        return self._stop_codon


class GTF_GENE(GTF_BASE, GET_GTF_ATTR):
    def __init__(self):
        super().__init__()
        self._feature = 'gene'
        self._gene_name = None
        self._gene_id = None
        self._biotype = None
        self._transcript_id_transcript_dic = {}

    def _init_by_gtf_line(self, seqid, source, left, right, score, ori, frame, attr_dic):
        frange = [(left, right)]
        _init_by_gtf_line_helper(self, seqid, source, frange, score, ori, frame, attr_dic)
        self._biotype = attr_dic.get("gene_biotype")

    def get_transcript_id_s(self):
        return list(self._transcript_id_transcript_dic.keys())

    def get_transcript_s(self):
        return list(self._transcript_id_transcript_dic.values())

    def get_transcript(self, transcript_id):
        return self._transcript_id_transcript_dic.get(transcript_id)


class GTF:
    def __init__(self):
        self._meta = None
        self._gene_id_gene_dic = {}
        self._gene_name_gene_id_dic = {}

    def __iter__(self):
        self.__start = 0
        self.__stop = len(self._gene_id_gene_dic)
        self.__gene_id_s = list(self._gene_id_gene_dic.keys())
        return self

    def __next__(self):
        if self.__start == self.__stop:
            raise StopIteration
        self.__start += 1
        return self._gene_id_gene_dic[self.__gene_id_s[self.__start - 1]]

    def get_gene_id_s(self):
        return list(self._gene_id_gene_dic.keys())

    def get_gene_s(self):
        return list(self._gene_id_gene_dic.values())

    def get_gene(self, gene_id):
        return self._gene_id_gene_dic.get(gene_id)

    def get_meta(self):
        return self._meta

    def get_gene_id_by_name(self, gene_name):
        return self._gene_name_gene_id_dic.get(gene_name)

def _read_gtf(gtf_file):
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
                    if type(pre_v) == list:
                        attr_dic[ele[0]] = pre_v + [ele[1]]
                    else:
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


def _merge_feature(data, feature, attr_key):
    """
    merge lines of a feature by a attributes.
    """
    if feature not in data:
        return

    data_fea = data[feature]
    group_dic = {}
    ignored_line_num = 0
    for line in data_fea:
        seqid = line[0]
        if seqid in group_dic:
            pass
        else:
            group_dic[seqid] = {}
        key = line[-1].get(attr_key)
        if key:
            if key not in group_dic[seqid]:
                group_dic[seqid][key] = [line]
            else:
                group_dic[seqid][key].append(line)
        else:
            ignored_line_num += 1
            #print(f"{attr_key} not found in {line}")
    
    if ignored_line_num:
        print(f"{ignored_line_num} entries in {feature} were ommited, because less {attr_key} attribute.", file=sys.stderr)

    grouped_data = []
    for seqid in group_dic:
        for key in group_dic[seqid]:
            seqid_g = []
            source_g = []
            left_g = []
            right_g = []
            score_g = []
            ori_g = []
            frame_g = []
            attr_g = []
            for line in group_dic[seqid][key]:
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
    """
    The funtion to read and parse gtf file. It create GTF,
    GTF_GENE, GTF_TRANSCRIPT, GTF_CDS, GTF_START_CODON,
    GTF_STOP_CODON object and link those object together.
    """
    meta, gtf_data = _read_gtf(gtf_file)
    #print(list(gtf_data.keys()))

    _merge_feature(gtf_data, "exon", "transcript_id")
    _merge_feature(gtf_data, "CDS", "protein_id")
    _merge_feature(gtf_data, "start_codon", "protein_id")
    _merge_feature(gtf_data, "stop_codon", "protein_id")

    gtf = GTF()
    gtf._meta = meta

    assert 'gene' in gtf_data
    gene_lines = gtf_data['gene']
    gene_id_gene_dic = {}
    gene_name_gene_id_dic = {}
    for line in gene_lines:
        gene = GTF_GENE()
        seqid, source, left, right, score, ori, frame, attr_dic = line
        gene._init_by_gtf_line(seqid, source, left, right, score, ori, frame, attr_dic)
        gene_name = gene.get_gene_name()
        gene_id = gene.get_gene_id()
        assert gene_id != None
        if not gene_name:
            gene_name = gene_id
        if gene_name not in gene_name_gene_id_dic:
            gene_name_gene_id_dic[gene_name] = [gene_id]
        else:
            gene_name_gene_id_dic[gene_name].append(gene_id)
        gene_id_gene_dic[gene_id] = gene
    gtf._gene_id_gene_dic = gene_id_gene_dic
    gtf._gene_name_gene_id_dic = gene_name_gene_id_dic

    if 'transcript' in gtf_data:
        transcript_lines = gtf_data['transcript']
        for line in transcript_lines:
            transcript = GTF_TRANSCRIPT()
            seqid, source, left, right, score, ori, frame, attr_dic = line
            transcript._init_by_gtf_line(seqid, source, left, right, score, ori, frame, attr_dic)
            gene_id = transcript.get_gene_id()
            transcript_id = transcript.get_transcript_id()
            gene = gene_id_gene_dic.get(gene_id)
            if gene:
                gene._transcript_id_transcript_dic[transcript_id] = transcript
            else:
                print(f"{gene_id} not contain a gene feature line.", file=sys.stderr)

    if 'exon' in gtf_data:
        exon_lines = gtf_data['exon']
        for line in exon_lines:
            exon = GTF_EXON()
            seqid, source, frange, score, ori, frame, attr_dic = line
            exon._init_by_gtf_line(seqid, source, frange, score, ori, frame, attr_dic)
            gene_id = exon.get_gene_id()
            transcript_id = exon.get_transcript_id()
            gene = gene_id_gene_dic.get(gene_id)
            if gene:
                if transcript_id in gene.get_transcript_id_s():
                    transcript = gene.get_transcript(transcript_id)
                    if transcript._exon:
                        print("warning, multiple exons", transcript_id)
                    transcript._exon = exon
                else:
                    transcript = GTF_TRANSCRIPT()
                    transcript._gene_id = gene_id
                    transcript._transcript_id = transcript_id
                    transcript._exon = exon
                    gene._transcript_id_transcript_dic[transcript_id] = transcript
            else:
                    print(f"{gene_id} not contain a gene feature line.", file=sys.stderr)

    if 'CDS' in gtf_data:
        cds_lines = gtf_data['CDS']
        for line in cds_lines:
            cds = GTF_CDS()
            seqid, source, frange, score, ori, frame, attr_dic = line
            cds._init_by_gtf_line(seqid, source, frange, score, ori, frame, attr_dic)
            gene_id = cds.get_gene_id()
            transcript_id = cds.get_transcript_id()
            gene = gene_id_gene_dic.get(gene_id)
            if gene:
                if transcript_id in gene.get_transcript_id_s():
                    transcript = gene.get_transcript(transcript_id)
                    if transcript._CDS:
                        print("warning, multiple cds", transcript_id)
                    transcript._CDS = cds
                else:
                    transcript = GTF_TRANSCRIPT()
                    transcript._gene_id = gene_id
                    transcript._transcript_id = transcript_id
                    transcript._CDS = cds
                    gene._transcript_id_transcript_dic[transcript_id] = transcript
            else:
                print(f"{gene_id} not contain a gene feature line.", file=sys.stderr)

    if 'start_codon' in gtf_data:
        start_codon_lines = gtf_data['start_codon']
        for line in start_codon_lines:
            start_codon = GTF_START_CODON()
            seqid, source, frange, score, ori, frame, attr_dic = line
            start_codon._init_by_gtf_line(seqid, source, frange, score, ori, frame, attr_dic)
            gene_id = start_codon.get_gene_id()
            transcript_id = start_codon.get_transcript_id()
            gene = gene_id_gene_dic.get(gene_id)
            if gene:
                if transcript_id in gene.get_transcript_id_s():
                    transcript = gene.get_transcript(transcript_id)
                    transcript._start_codon = start_codon
                else:
                    transcript = GTF_TRANSCRIPT()
                    transcript._gene_id = gene_id
                    transcript._transcript_id = transcript_id
                    transcript._start_codon = start_codon
                    gene._transcript_id_transcript_dic[transcript_id] = transcript
            else:
                print(f"{gene_id} not contain a gene feature line.", file=sys.stderr)

    if 'stop_codon' in gtf_data:
        stop_codon_lines = gtf_data['stop_codon']
        for line in stop_codon_lines:
            stop_codon = GTF_START_CODON()
            seqid, source, frange, score, ori, frame, attr_dic = line
            stop_codon._init_by_gtf_line(seqid, source, frange, score, ori, frame, attr_dic)
            gene_id = stop_codon.get_gene_id()
            transcript_id = stop_codon.get_transcript_id()
            gene = gene_id_gene_dic.get(gene_id)
            if gene:
                if transcript_id in gene.get_transcript_id_s():
                    transcript = gene.get_transcript(transcript_id)
                    transcript._stop_codon = stop_codon
                else:
                    transcript = GTF_TRANSCRIPT()
                    transcript._gene_id = gene_id
                    transcript._transcript_id = transcript_id
                    transcript._stop_codon = stop_codon
                    gene._transcript_id_transcript_dic[transcript_id] = transcript
            else:
                print(f"{gene_id} not contain a gene feature line.", file=sys.stderr)
    return gtf


def test_gtf(gtf_file):
    """
    test gtf parse.
    """
    parse_gtf(gtf_file)


if __name__ == "__main__":
    import sys
    test_gtf(sys.argv[1])
