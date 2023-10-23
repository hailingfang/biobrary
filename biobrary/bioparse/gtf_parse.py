import sys
import re


class Gtf_gene:
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
        self.geneid = None
        self.data = None
        self.seqname = None
        self.source = None
        self.left = None
        self.right = None
        self.ori = None
        self.gene_biotype = None
        self.attr = None

        data = {}
        for line in gene_raw_dt:
            seqname, source, feature, left, right, score, ori, frame, attrib = line
            attrib = [ele.strip() for ele in attrib.split(";")[:-1]]
            attrib_dic = {}
            for ele in attrib:
                idx = ele.index(" ")
                key = ele[:idx]
                value = ele[idx+1:].strip('"')
                if len(value) > 0:
                    if value in attrib_dic:
                        attrib_dic[key].append(value)
                    else:
                        attrib_dic[key] = [value]

            if feature == "gene":
                assert len(attrib_dic["gene_id"]) == 1
                assert frame == "."
                self.geneid = attrib_dic["gene_id"][0]
                self.seqname = seqname
                self.source = source
                self.left = int(left)
                self.right = int(right)
                self.ori = ori
                self.gene_biotype = attrib_dic["gene_biotype"][0]
                self.attr = attrib_dic
            elif feature == "transcript":
                assert frame == "."
                geneid = attrib_dic["gene_id"][0]
                transid = attrib_dic["transcript_id"][0]
                assert geneid == self.geneid
                assert transid not in data
                data[transid] = {}
                data[transid]["seqname"] = seqname
                data[transid]["source"] = source
                data[transid]["left"] = int(left)
                data[transid]["right"] = int(right)
                data[transid]["ori"] = ori
                data[transid]["transcript_biotype"] = \
                    attrib_dic["transcript_biotype"][0]
                data[transid]["attr"] = attrib_dic
            elif feature == "exon":
                assert frame == "."
                geneid = attrib_dic["gene_id"][0]
                transid = attrib_dic["transcript_id"][0]
                exon_number = attrib_dic["exon_number"][0]
                assert geneid == self.geneid
                if transid in data and ("exon" not in data[transid]):
                    data[transid]["exon"] = {}
                else:
                    data[transid] = {}
                    data[transid]["exon"] = {}
                data[transid]["exon"][exon_number] = {}
                data[transid]["exon"][exon_number]["left"] = int(left)
                data[transid]["exon"][exon_number]["right"] = int(right)
                data[transid]["exon"][exon_number]["ori"] = ori
            elif feature == "CDS":
                geneid = attrib_dic["gene_id"][0]
                print(geneid)
                transid = attrib_dic["transcript_id"][0]
                proteinid = attrib_dic.get("protein_id")
                if proteinid:
                    proteinid = proteinid[0]
                else:
                    proteinid = None
                exon_number = attrib_dic["exon_number"][0]
                assert geneid == self.geneid
                assert transid in data
                if "CDS" not in data[transid]:
                    data[transid]["CDS"] = [proteinid, {}]
                data[transid]["CDS"][1][exon_number] = {}
                data[transid]["CDS"][1][exon_number]["left"] = int(left)
                data[transid]["CDS"][1][exon_number]["right"] = int(right)
                data[transid]["CDS"][1][exon_number]["ori"] = ori
                data[transid]["CDS"][1][exon_number]["frame"] = frame
            elif feature == "start_codon":
                geneid = attrib_dic["gene_id"][0]
                transid = attrib_dic["transcript_id"][0]
                proteinid = attrib_dic.get("protein_id")
                if proteinid:
                    proteinid = proteinid[0]
                else:
                    proteinid = None
                exon_number = attrib_dic["exon_number"][0]

                assert geneid == self.geneid
                assert transid in data
                assert proteinid == data[transid]["CDS"][0]
                if "start_codon" not in data[transid]:
                    data[transid]["start_codon"] = [proteinid, {}]
                data[transid]["start_codon"][1][exon_number] = {}    
                data[transid]["start_codon"][1][exon_number]["left"] = int(left)
                data[transid]["start_codon"][1][exon_number]["right"] = int(right)
                data[transid]["start_codon"][1][exon_number]["ori"] = ori
                data[transid]["start_codon"][1][exon_number]["frame"] = frame
            elif feature == "stop_codon":
                geneid = attrib_dic["gene_id"][0]
                transid = attrib_dic["transcript_id"][0]
                proteinid = attrib_dic.get("protein_id")
                if proteinid:
                    proteinid = proteinid[0]
                else:
                    proteinid = None
                exon_number = attrib_dic["exon_number"][0]
                assert geneid == self.geneid
                assert transid in data
                assert proteinid == data[transid]["CDS"][0]
                if "stop_codon" not in data[transid]:
                    data[transid]["stop_codon"] = [proteinid, {}]
                data[transid]["stop_codon"][1][exon_number] = {}
                data[transid]["stop_codon"][1][exon_number]["left"] = int(left)
                data[transid]["stop_codon"][1][exon_number]["right"] = int(right)
                data[transid]["stop_codon"][1][exon_number]["ori"] = ori
                data[transid]["stop_codon"][1][exon_number]["frame"] = frame
            else:
                print(f"Warning, {feature} not recognized", file=sys.stderr)
        self.data = data
        self._check_transcript_data()

        

    def _check_transcript_data(self):
        """
        some gtf file do not have transcript line, make up info for such files
        """
        for transid in self.data:
            if "left" not in self.data[transid]:
                lefts = []
                rights = []
                oris = []
                bio_type = []
                if "exon" in self.data[transid]:
                    for exon_num in self.data[transid]["exon"]:
                        lefts.append(self.data[transid]["exon"][exon_num]["left"])
                        rights.append(self.data[transid]["exon"][exon_num]["right"])
                        oris.append(self.data[transid]["exon"][exon_num]["ori"])
                        bio_type.append(self.data[transid]["exon"][exon_num]\
                                        .get("gbkey"))
                    lefts.sort()
                    rights.sort()
                    assert len(set(oris)) == 1
                    assert len(set(bio_type)) == 1
                    self.data[transid]["left"] = lefts[0]
                    self.data[transid]["right"] = rights[-1]
                    self.data[transid]["ori"] = oris[0]
                    self.data[transid]["transcript_biotype"] = bio_type[0]
                    self.data[transid]["seqname"] = self.seqname
                    self.data[transid]["source"] = self.source
                else:
                    print("Error,transcript have no exon", file=sys.stderr)
        
        
    def get_geneid(self):
        return self.geneid
    
    def get_geneid_range(self):
        return (self.left, self.right, self.ori)
    
    def get_gene_type(self):
        return self.gene_biotype

    def if_protein_coding(self):
        if self.gene_biotype == "protein_coding":
            return True
        return False

    def get_transcript(self):
        trans_data = []
        for transid in self.data:
            data_entry = []
            data_entry.append(transid)
            data_entry.appen([self.data[transid]["left"],
                              self.data[transid]["right"]])
            data_entry.append(self.data[transid]["ori"])
            data_entry.append(self.data[transid]["transcript_biotype"])
            trans_data.append(data_entry)
        
        return trans_data


    def get_exon(self):
        exon_data = []
        for transid in self.data:
            if "exon" in self.data[transid]:
                data_entry = []
                data_entry.append(transid)
                data_entry.append([])
                for exon_num in self.data[transid]["exon"]:
                    data_entry[1].append([exon_num,
                                self.data[transid]["exon"][exon_num]["left"],
                                self.data[transid]["exon"][exon_num]["right"]])
                data_entry.append(self.data[transid]["exon"][exon_num]["ori"])
                exon_data.append(data_entry)
            
        return exon_data
    

    def get_CDS(self):
        CDS_data = []
        for transid in self.data:
            if "CDS" in self.data[transid]:
                data_entry = []
                data_entry.append(transid)
                data_entry.append(self.data[transid]["CDS"][0])
                data_entry.append([])
                for exon_num in self.data[transid]["CDS"][1]:
                    data_entry[2].append([exon_num,
                                self.data[transid]["CDS"][exon_num]["left"],
                                self.data[transid]["CDS"][exon_num]["left"]])
                data_entry.append(self.data[transid]["CDS"][exon_num]["ori"])
                CDS_data.append(data_entry)
        return CDS_data

        
    def get_start_codon(self):
        start_codon = []
        for transid in self.data:
            if "start_codon" in self.data[transid]:
                start_codon.append(self.data[transid]["start_codon"][0])
                start_codon.append([[self.data[transid]["start_condn"][1]["left"]],\
                                    self.data[transid]["start_codon"][1]["right"]])
                start_codon.append(self.data[transid]["start_codon"][1]["ori"])
                start_codon.append(self.data[transid]["start_codon"][1]["frame"])
        return start_codon
                

    def get_stop_codon(self):
        stop_codon = []
        for transid in self.data:
            if "stop_codon" in self.data[transid]:
                stop_codon.append(self.data[transid]["stop_codon"][0])
                stop_codon.append([[self.data[transid]["start_condn"][1]["left"]],\
                                    self.data[transid]["stop_codon"][1]["right"]])
                stop_codon.append(self.data[transid]["stop_codon"][1]["ori"])
                stop_codon.append(self.data[transid]["stop_codon"][1]["frame"])
        return stop_codon
    



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

            if feature == "gene":
                index += 1
                key = index
                gtf_gene_data[key] = [line]
                continue
            gtf_gene_data[key].append(line)
        
        for key in gtf_gene_data:
            self.data.append(Gtf_gene(gtf_gene_data[key]))
        

    def __iter__(self):
        pass


    def __next__(self):
        pass


    def get_geneids(self):
        geneids = []
        for gene in self.data:
            geneids.append(gene.get_gene())
        return geneids


    def get_gene(self, geneid):
        return self.data[self.geneid_index[geneid]]


if __name__ == "__main__":
    import sys
    gtf = GTF(sys.argv[1])
