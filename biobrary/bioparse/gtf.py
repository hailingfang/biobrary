"""
Founctions and Class to parse GTF file
"""

import re
import gzip
import sys
from biobrary.misc import merge_islands


re_temp = re.compile(r'\s?(.+?)\s"(.+?)";')


class GTF_ENTRY:
    def __init__(self):
        self._seq = None
        self._source = None
        self._feature = None
        self._left = None
        self._right = None
        self._score = None
        self._ori = None
        self._frame = None
        self._attr_dic = {}

    def create_entry_by_gtf_line(self, line):
        line = line.split('\t')
        self._seqname = line[0]
        self._source = line[1]
        self._feature = line[2]
        self._left = int(line[3])
        self._right = int(line[4])
        self._score = line[5]
        self._ori = line[6]
        self._frame = line[7]
        attr_dic = {}
        attr = re_temp.findall(line[8])
        for ele in attr:
            if ele[0] not in attr_dic:
                attr_dic[ele[0]] = [ele[1]]
            else:
                attr_dic[ele[0]].append(ele[1])
        self._attr_dic = attr_dic

    def get_seqname(self):
        return self._seqname
    
    def get_source(self):
        return self._source
    
    def get_feature(self):
        return self._feature

    def get_left(self):
        return self._left
    
    def get_right(self):
        return self._right

    def get_score(self):
        return self._score
    
    def get_ori(self):
        return self._ori
    
    def get_frame(self):
        return self._frame
    
    def get_attr_dic(self):
        return self._attr_dic
    
    def get_attr(self, attr):
        return self._attr_dic.get(attr)
    
    def get_position(self):
        return (self._seqname, self._left, self._right, self._ori)


class GTF_TREE:
    def __init__(self):
        self.parent = None
        self.child = []


class GTF:
    def __init__(self):
        self._meta = None
        self._entries = []
        self._seqname_s = None
        self._source_s = None
        self._feature_type_s = None

    def get_meta(self):
        return self._meta
    
    def get_entries(self):
        return self._entries

    def get_entries_by_seqname(self, seqname):
        date_out = []
        for ent in self._entries:
            if ent.get_seqname() == seqname:
                date_out.append(ent)
        return date_out

    def get_entries_by_soure(self, source):
        data_out = []
        for ent in self._entries:
            if ent.get_source() == source:
                data_out.append(ent)
        return data_out
    
    def get_entries_by_feature(self, feature):
        data_out = []
        for ent in self._entries:
            if ent.get_feature() == feature:
                data_out.append(ent)
        return data_out

    def get_entries_by_left(self, border, side='right'):
        data_out = []
        if side == 'right':
            for ent in self._entries:
                if ent.get_left() >= border:
                    data_out.append(ent)
        elif side == 'left':
            for ent in self._entries:
                if ent.get_left() <= border:
                    data_out.append(ent)
        else:
            pass
        return data_out

    def get_entries_by_right(self, border, side='left'):
        data_out = []
        if side == 'left':
            for ent in self._entries:
                if ent.get_right() <= border:
                    data_out.append(ent)
        elif side == 'right':
            for ent in self._entries:
                if ent.get_right() >= border:
                    data_out.append(ent)
        else:
            pass
        return data_out

    def get_entries_by_score(self, score):
        data_out = []
        for ent in self._entries:
            if ent.get_score() == score:
                data_out.append(ent)
        return data_out

    def get_entries_by_ori(self, ori):
        data_out = []
        for ent in self._entries:
            if ent.get_ori() == ori:
                data_out.append(ent)
        return data_out

    def get_entries_by_frame(self, frame):
        data_out = []
        for ent in self._entries:
            if ent.get_frame() == frame:
                data_out.append(ent)
        return data_out

    def get_entries_by_attr(self, attr_name, attr_value):
        data_out = []
        for ent in self._entries:
            if ent.get_attr(attr_name) == attr_value:
                data_out.append(ent)
        return data_out

    def construct_gene_entries_by_transcript(self):
        transcript_entries = self.get_entries_by_feature('transcript')
        group_dic = {}
        for trans in transcript_entries:
            gene_id = trans.get_attr('gene_id')
            if gene_id in group_dic:
                group_dic[gene_id].apppend(trans)
            else:
                group_dic[gene_id] = [trans]

        for gene_id in group_dic:
            ent = GTF_ENTRY()
            gene_name = []
            seqname_s = set()
            ori_s = set()
            ent_range_s = []
            for trans in group_dic[gene_id]:
                gene_name.append(trans.get_attr('gene_name'))
                position = trans.get_position()
                seqname_s.add(position[0])
                ori_s.add(position[3])
                ent_range_s.append(position[1: 3])
            assert len(seqname_s) == 1
            assert len(ori_s) == 1
            merged_islands = merge_islands(ent_range_s)
            merged_islands.sort(key=lambda x:x[0])
            ent._seq = list(seqname_s)[0]
            ent._source = 'GTF_PARSER'
            ent._feature = 'gene'
            ent._left = merged_islands[0][0]
            ent._right = merged_islands[-1][-1]
            ent._score = '.'
            ent._ori = list(ori_s)[0]
            ent._frame = '.'
            ent._attr_dic = {'gene_id': gene_id}
            self._entries.append(ent)
        self._source_s.add('GTF_PARSER')
        self._feature_type_s.add('gene')



    def struct_fea_by_position(self, fea_parent, fea_child):
        pass


    def struc_fea_by_attribute(self, fea_parent, fea_child):
        pass


    def struct_feature_s(self, feature_s, method_s):
        pass



def parse_gtf(gtf_file):
    """
    """
    gtf = GTF()
    meta = []
    entry = []
    seqname_s = set()
    source_s = set()
    feature_type_s = set()

    if gtf_file.endswith('.gz'):
        for line in gzip.open(gtf_file, 'r'):
            line = line.decode().rstrip()
            if line[0] != '#':
                line_split = line.split('\t')
                seqname_s.add(line_split[0])
                source_s.add(line_split[1])
                feature_type_s.add(line_split[2])
                ent = GTF_ENTRY()
                ent.create_entry_by_gtf_line(line)
                entry.append(ent)
            else:
                meta.append(line)
    else:
        for line in open(gtf_file, 'r'):
            line = line.rstrip()
            if line[0] != '#':
                line_split = line.split('\t')
                seqname_s.add(line_split[0])
                source_s.add(line_split[1])
                feature_type_s.add(line_split[2])
                ent = GTF_ENTRY()
                ent.create_entry_by_gtf_line(line)
                entry.append(ent)
            else:
                meta.append(line)
    gtf._meta = meta
    gtf._entries = entry
    gtf._seqname_s = seqname_s
    gtf._source_s = source_s
    gtf._feature_type_s = feature_type_s

    return gtf


def test_gtf(gtf_file):
    """
    test gtf parse.
    """
    gtf = parse_gtf(gtf_file)
    entries = gtf.get_entries_by_feature('gene')
    for ent in entries:
        print(ent.get_position())



if __name__ == "__main__":
    import sys
    test_gtf(sys.argv[1])
