import sys
import re

class FEATURES:
    def __init__(self, features_lines):
        if (len(features_lines)) != 1:
            print("Error, one LUCUS should only have one block of features data", file=sys.stderr)
            exit()
        features_lines = features_lines[0]
        self.data = FEATURES._split_features_lines(self, features_lines)
    
    def _split_features_lines(self, features_lines):
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
        re_temp = re.compile(r'^/(\w+)=?"?(.+?)?"?$')
        for block in blocks:
            #[feature, orientation, range, {key: value}]
            dt = [None, None, None, {}]
            first_l = block[0]
            feature, strand_range = first_l.split()
            dt[0] = feature
            if strand_range.startswith("compl"):
                dt[1] = "-"
            else:
                dt[1] = "+"

            dt[2] = re_temp_f_range.match(strand_range).groups()[0]
            
            for l in block[1:]:
                re_res = re_temp.match(l)
                key, value = re_res.groups()
                if key not in dt[3]:
                    dt[3][key] = [value]
                else:
                    dt[3][key].append(value)
            new_blocks.append(dt)
        blocks = new_blocks
        
        new_blocks = []
        gene_group = []
        gene_range = (0, 0)
        for block in blocks:
            feature_name = block[0]
            f_range = block[2].split(",")
            f1 = f_range[0].split("..")[0]
            f2 = f_range[-1].split("..")[-1]
            if f1[0].isdigit():
                f_start = int(f1)
            else:
                f_start = int(f1[1:])

            if f2[0].isdigit():
                f_end = int(f2)
            else:
                f_end = int(f2[1:])
            
            if f_start >= gene_range[0] and f_end <= gene_range[1]:
                gene_group.append(block)
            else:
                new_blocks.append(gene_group)
                gene_group = [block]
                if feature_name == "gene":
                    gene_range = [f_start, f_end]
        blocks = new_blocks[1:]
        return blocks


class LOCUS:
    def __init__(self, data_lines):
        self.data = data_lines
        self.locus_info = None
        self.defination = None
        self.accession = None
        self.version = None
        self.dblink = None
        self.keywords = None
        self.source = None
        self.reference = None
        self.comment = None
        self.features = None
        self.contig = None
        self.origin = None

    def _split_locus_anno_info(self, data_lines):
        data_blocks = {}
        tmp = []
        block = []
        for l in data_lines:
            if len(l) == 0:
                continue
            if l[0] != " ":
                tmp.append(block)
                block = []
            block.append(l)
        tmp.append(block)
        tmp = tmp[1:]

        for block in tmp:
            block_type = block[0].split()[0]
            if block_type not in data_blocks:
                data_blocks[block_type] = [block]
            else:
                data_blocks[block_type].append(block)
        return data_blocks
   

    def _parse_locus_info(self, dt_lines):
        self.locus_info = dt_lines[0][0].split()[1]

    def _parse_definition(self, dt_lines):
        self.definition = None

    def _parse_accession(self, dt_lines):
        self.accession = None

    def _parse_version(self, dt_lines):
        self.version = None

    def _parse_dblink(self, dt_lines):
        self.dbline = None

    def _parse_keywords(self, dt_lines):
        self.keywords = None

    def _parse_source(self, dt_lines):
        self.source = None

    def _parse_reference(self, dt_lines):
        self.reference = None

    def _parse_comment(self, dt_lines):
        self.comment = None

    def _parse_features(self, dt_lines):
        self.features = FEATURES(dt_lines)
    
    def _parse_contig(self, dt_lines):
        self.contig = None

    def _parse_origin(self, dt_lines):
        self.origin = ""
        ll = ""
        for line in dt_lines[0][1:]:
            ll += "".join(line.split()[1:])
        self.origin = ll

    def parse_locus(self):
        data_blocks  = LOCUS._split_locus_anno_info(self, self.data) 

        for key in data_blocks:
            dt_lines = data_blocks[key]
            if key == "LOCUS":
                self._parse_locus_info(dt_lines)
            elif key == "DEFINITION":
                self._parse_definition(dt_lines)
            elif key == "ACCESSION":
                self._parse_accession(dt_lines)
            elif key == "VERSION":
                self._parse_version(dt_lines)
            elif key == "DBLINK":
                self._parse_dblink(dt_lines)
            elif key == "KEYWORDS":
                self._parse_keywords(dt_lines)
            elif key == "SOURCE":
                self._parse_source(dt_lines)
            elif key == "REFERENCE":
                self._parse_reference(dt_lines)
            elif key == "COMMENT":
                self._parse_comment(dt_lines)
            elif key == "FEATURES":
                self._parse_features(dt_lines)
            elif key == "CONTIG":
                self._parse_contig(dt_lines)
            elif key == "ORIGIN":
                self._parse_origin(dt_lines)
            else:
                print("The INFOR type not known.", key, file=sys.stderr)
    

    def print_selected_seq(self, start, end):
        print(self.origin[start - 1, end])


    def print_gene_data(self):
        data = self.features.data
        seq = self.origin
        locus = self.locus_info
        comp_dic = {"a":"t", "g": "c", "c": "g", "t": "a", "n": "n"}

        for dt in data:
            dt_dic = {}
            if dt[0][0] == "gene":
                for ele in dt:
                    if ele[0] not in dt_dic:
                        dt_dic[ele[0]] = [ele]
                    else:
                        dt_dic[ele[0]].append(ele)
                gene_dt = dt_dic.pop("gene")[0]
                if "pseudo" in gene_dt[3]:
                    gene_ty = "pseudo"
                else:
                    gene_ty = "truegene"
                
                #print data
                gene_id = gene_dt[3]["gene"][0]
                print("\t".join([">gene", gene_id, gene_dt[1], gene_dt[2], gene_ty, locus]))
                gene_range = gene_dt[2].split("..")
                if gene_range[0][0].isdigit():
                    g_s = int(gene_range[0])
                else:
                    g_s = int(gene_range[0][1:])
                if gene_range[1][0].isdigit():
                    g_e = int(gene_range[1])
                else:
                    g_e = int(gene_range[1][1:])
                gene_seq = seq[g_s - 1: g_e]
                if gene_dt[1] == "-":
                    gene_seq = gene_seq[::-1]
                    gene_seq = "".join([comp_dic[e] for e in gene_seq])
                print(gene_seq)

                if "mRNA" in dt_dic:
                    mRNA_dt = dt_dic.pop("mRNA")
                    for mRNA_ele in mRNA_dt:
                        trans_id = mRNA_ele[3].get("transcript_id")
                        if trans_id:
                            trans_id = trans_id[0] + "|" + gene_id
                        else:
                            trans_id = "-" + "|" + gene_id
                        print("\t".join(["@mRNA", trans_id, mRNA_ele[1], mRNA_ele[2]]))
                        seq_range = mRNA_ele[2].split(",")
                        seq_range = [e.split("..") for e in seq_range]
                        tmp = []
                        for ele in seq_range:
                            if len(ele) == 2:
                                start, end = ele
                                if start[0].isdigit():
                                    start = int(start)
                                else:
                                    start = int(start[1:])
                                if end[0].isdigit():
                                    end = int(end)
                                else:
                                    end = int(end[1:])
                                tmp.append([start, end])
                            else:
                                start = ele[0]
                                if start[0].isdigit():
                                    start = int(start)
                                else:
                                    start = int(start[1:])
                                tmp.append([start])

                        seq_range = tmp
                        mRNA_seq = ""
                        for ele in seq_range:
                            if len(ele) == 2:
                                mRNA_seq += seq[ele[0] - 1: ele[1]]
                            else:
                                mRNA_seq += seq[ele[0] - 1]
                        if mRNA_ele[1] == "-":
                            mRNA_seq = mRNA_seq[::-1]
                            mRNA_seq = "".join([comp_dic[e] for e in mRNA_seq])
                        print(mRNA_seq)
                 
                if "CDS" in dt_dic:
                    cds_dt = dt_dic.pop("CDS")
                    for cds_ele in cds_dt:
                        prot_id = cds_ele[3].get("protein_id")
                        if prot_id:
                            prot_id = prot_id[0] + "|" + gene_id
                        else:
                            prot_id = "-" + "|" + gene_id
                        print("\t".join(["@CDS", prot_id, cds_ele[1], cds_ele[2], cds_ele[3]["codon_start"][0]]))
                        seq_range = cds_ele[2].split(",")
                        seq_range = [e.split("..") for e in seq_range]
                        cds_seq = ""
                        for ele in seq_range:
                            if len(ele) == 2:
                                start, end = ele
                                if start[0].isdigit():
                                    start = int(start)
                                else:
                                    start = int(start[1:])
                                if end[0].isdigit():
                                    end = int(end)
                                else:
                                    end = int(end[1:])
                                cds_seq += seq[start - 1: end]
                            if len(ele) == 1:
                                start = ele[0]
                                if start[0].isdigit():
                                    start = int(start)
                                else:
                                    start = int(start[1:])

                                cds_seq += seq[start - 1]
                        if cds_ele[1] == "-":
                            cds_seq = cds_seq[::-1]
                            cds_seq = "".join([comp_dic[e] for e in cds_seq])
                        print(cds_seq)
                        prot_seq = cds_ele[3].get("translation")
                        if prot_seq:
                            print(prot_seq[0])
                        else:
                            print("-")
                                    



class GBFF_parser:
    """
    Parse gbff file.
    gbff_data = GBFF_parse(gbff_file_name)
    """
    def __init__(self, gbff_fname):
        self.gbff_fname = gbff_fname
        self.locus = set()
        self.loaded_locus_id = []
        self.locus_data = []
    
    def get_locus_ids(self):
        if self.locus:
            return self.locus

        fin = open(gbff_fname, "r")
        for l in fin:
            if l[:5] == "LOCUS":
                self.locus.add(l.rstrip().split()[1])
        fin.close()
        return self.locus

    def parse_data(self, item_locus=[]):
        self.loaded_locus = item_locus

        data_split_by_locus = []
        fin = open(self.gbff_fname, "r")
        dt_block = []
        for line in fin:
            if line[:2] == "//":
                locus_id = dt_block[0].split()[1]
                if locus_id in item_locus or len(item_locus) == 0:
                    data_split_by_locus.append(dt_block)
                dt_block = []
                continue
            dt_block.append(line.rstrip())
        fin.close()
        
        for dt_block in data_split_by_locus:
            self.locus_data.append(LOCUS(dt_block))

        return self.locus_data


if __name__ == "__main__":
    gbkf = sys.argv[1]
    gbk_dt = GBFF_parse(gbkf)
    #locuses = gbk_dt.parse_data(["NC_007112"])
    locuses = gbk_dt.parse_data()

    for locus in locuses:
        locus.parse_locus() 
        locus.print_gene_data()

