import sys
import re

class GBFF_locus:
    def __init__(self, locus_lines):
        
        if len(locus_lines) != 1:
            print("Warning, The LOCUS line should be one line.", file=sys.stderr)
        locus_lines = locus_lines[0]
        lin_fields = locus_lines.split()
        self.locus = lin_fields[1]
        self.seq_length = int(lin_fields[2])
    
    def __str__(self):
        return "\t".join([self.locus, str(self.seq_length)])


class GBFF_definition:
    def __init__(self, defin_lines):
        if len(defin_lines) != 1:
            print("Warning, The DEFINITION should be one line.", file=sys.stderr)
        defin_lines = defin_lines[0]
        self.definition = defin_lines


class GBFF_accession:
    def __init__(self, access_lines):
        if len(access_lines) != 1:
            print("Warning, The ACCESSION should contain only one line.", file=sys.stderr)
        access_lines = access_lines[0]
        self.accession = access_lines.split()[1]


class GBFF_version:
    def __init__(self, version_lines):
        if len(version_lines) != 1:
            print("Warning, The VERSION should contain only one line.", file=sys.stderr)
        version_lines = version_lines[0]
        self.version = version_lines.split()[1]


class GBFF_dlink:
    def __init__(self, dlink_lines):
        self.dlink = {}
        re_tmpt = re.compile(r".+\s(\w+):\s(.+)$")
        for l in dline_lines:
            re_res = re_tmpt.match(l)
            if re_res:
                key, value = re_res.groups()
                self.dlink[key] = value


class GBFF_keywords:
    def __init__(self, keywords_lines):
        if len(keywords_lines) != 1:
            print("Warning, The keywords lines should only have one line.", file=sys.stderr)
        keywords_lines = keywords_lines[0]
        self.keywords = keywords_lines.split()[1]


class GBFF_source:
    def __init__(self, source_lines):
        self.organism = None
        re_tmpt = re.compile(r"\s{2}ORGANISM\s(.+)$")
        for l in source_lines:
            re_res = re_tmpt.match(l)
            if re_res:
                self.organism = re_res.groups()[0]


class GBFF_reference:
    def __init__(self, reference_lines)
        self.reference = None


class GBFF_features:
    def _split_features_lines_(self, features_lines):
        blocks = []
        tmp = []
        for l in features_lines[1:]:
            if l[5] != " ":
                blocks.append(tmp)
                tmp = []
            tmp.append(l)
        blocks.append(tmp)
        blocks = blocks[1:]
        
        new_blocks = []
        for block in blocks:
            tmp = []
            ll = ""
            for l in block:
                if l[5] != " " or l.strip().startswith("/"):
                    tmp.append(ll)
                    ll = l.strip()
                else:
                    ll += l.strip()
            tmp.append(ll)
            tmp = tmp[1:]
            new_blocks.append(tmp)
        
        blocks = new_blocks
        new_blocks = []
        re_temp_f_range = re.compile(r'^(?:complement)?\(?(?:join)?\(?(.+?)\)?\)?$')
        re_temp = re.compile(r'^/(.+)="(.+)"$')
        for block in blocks:
            dt = [None, None, None, {}]
            first_l = block[0]
            feature, strand_range = first_l.split()
            dt[0] = feature
            if strand_range.startswith("compl"):
                dt[1] = "-"
            else:
                dt[1] = "+"

            dt[2] = re_temp_f_range.match(strand_range).groups()[0]
            
            for l in block:
                re_res = re_temp.match(l)
                key, value = re_res.groups()
                if key not dt[3]:
                    dt[3][key] = [value]
                else:
                    dt[3][key].append(value)
            new_blocks.append(dt)

        blocks = new_blocks
        new_blocks = []
        group = []
        group_range = (0, 0)
        for block in blocks:
            feature_name = block[0]
            f_range = block[2].split("..")
            if f_range[0][0].isdigit():
                f_start = int(f_range[0])
            else:
                f_start = int(f_range[0][1:])
            if f_range[1][0].isdigit():
                f_end = int(f_range[1])
            else:
                f_end = int(f_range[1][1:])
            

            if f_start >= group_range[0] and f_end <= group_range[1]:
                group.append(block)
            else:
                new_blocks.append(group)
                group = [block]
                group_range = [f_start, f_end]
        blocks = new_blocks
        
        return blocks

    def __init__(self, features_lines):
        self.features = GBFF_features._split_features_lines_(self, features_lines)
       

class GBFF_parse:
    def __init__(self, gbff_fname):
        self.gbff_fname = gbff_fname
        self.locus = []
        fin = open(gbff_fname, "r")
        for l in fin:
            if l[:5] == "LOCUS":
                self.locus.append(l.rstrip().split()[1])
    
#    def load_locus(self, locus):
#        fin = open(self.gbff_fname)
#        dt_block = []
#        swith = 0
#        
#        for lin in fin:
#            if lin[:5] == "LOCUS":
#                lin_fields = lin.rstrip().split()
#                if lin_fields[1] == locus:
#                    swith = 1
#                elif swith:
#                    break
#            if swith:
#                dt_block.append(lin.rstrip())
#                
#            
