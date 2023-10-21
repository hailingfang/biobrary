import sys
import re


class Gtf_gene:
    def __init__(self, gene_raw_dt):
        self.gene_id = None
        self.gene = None
        self.transcript = None
        self.exon = None
        self.CSD = None
        self.start_stop = None
        
        gene_dt = {}
        transcript_dt = {}
        exon_dt = {}
        cds_dt = {}
        start_stop_dt = {}


        for line in gene_raw_dt:
            seqname, source, feature, left, right, score, ori, frame, attrib = line
            attrib = [ele.strip() for ele in attrib.split(";")][:-1]
            
            attrib_dic = {"db_xref": [], "gene_synonym": [], "transl_except": [], "inference": []}
            multi_attr = ["db_xref", "gene_synonym", "inference", "transl_except"]
            for ele in attrib:
                idx = ele.index(" ")
                key = ele[:idx]
                value = ele[idx+1:].strip('"')
                if len(value) > 0:
                    if key in attrib_dic and key in multi_attr:
                        attrib_dic[key].append(value)
                    elif key not in attrib_dic:
                        attrib_dic[key] = value
                    else:
                        attrib_dic[key] = value
                        print(f"To developer: {key} is a new attributate", file=sys.stderr)


            if feature == "gene":
                gene_dt["gene_id"] = attrib_dic["gene_id"]
                gene_dt["refseq"] = seqname
                gene_dt["source"] = source
                gene_dt["left"] = int(left)
                gene_dt["right"] = int(right)
                gene_dt["score"] = score
                gene_dt["ori"] = ori
                gene_dt["attr"] = {}

                gene_dt["attr"]["gene"] = attrib_dic["gene"]
                gene_dt["attr"]["gene_biotype"] = attrib_dic["gene_biotype"]
                gene_dt["attr"]["gbkey"] = attrib_dic["gbkey"]
                gene_dt["attr"]["db_xref"] = attrib_dic["db_xref"]
                gene_dt["attr"]["description"] = attrib_dic.get("description", None)
                gene_dt["attr"]["pseudo"] = attrib_dic.get("pseudo", "false")
                gene_dt["attr"]["gene_synonym"] = attrib_dic.get("gene_synonym", None)

                self.gene_id = attrib_dic["gene_id"]
                
            
            elif feature == "transcript":
                trans_dt_ele = {}
                transid = attrib_dic["transcript_id"]
                gbkey = attrib_dic["gbkey"]
                trans_biotype = attrib_dic["transcript_biotype"]
                product = attrib_dic.get("product", None)
                model_evidence = attrib_dic.get("model_evidence", None)

                trans_dt_ele["trans_id"] = transid
                trans_dt_ele["left"] = int(left)
                trans_dt_ele["right"] = int(right)
                trans_dt_ele["attr"] = {}

                trans_dt_ele["attr"]["gbkey"] = gbkey
                trans_dt_ele["attr"]["transript_biotype"] = trans_biotype
                trans_dt_ele["attr"]["product"] = product
                trans_dt_ele["attr"]["model_evidence"] = model_evidence
                
                transcript_dt[transid] = trans_dt_ele


            elif feature == "exon":
                transid = attrib_dic["transcript_id"]
                exon_number = attrib_dic["exon_number"]
                if transid not in exon_dt:
                    exon_dt[transid] = {"exon_range": [(int(left), int(right))], "exon_number": [int(exon_number)]}
                else:
                    exon_dt[transid]["exon_range"].append((int(left), int(right)))
                    exon_dt[transid]["exon_number"].append(int(exon_number))


            elif feature == "CDS":
                transid = attrib_dic["transcript_id"]
                proteinid = attrib_dic.get("protein_id", None)
                exon_number = attrib_dic["exon_number"]
                if transid not in cds_dt:
                    cds_dt[transid] = {"CDS_range": [(int(left), int(right))], "protein_id": proteinid, "frame": [int(frame)], "exon_number": [int(exon_number)]}
                else:
                    cds_dt[transid]["CDS_range"].append((int(left), int(right)))
                    cds_dt[transid]["frame"].append(int(frame))
                    cds_dt[transid]["exon_number"].append(int(exon_number))
                start_stop_dt[transid] = {"start_codon": None, "stop_codon": None, "protein_id": proteinid}


            elif feature == "start_codon":
                transid = attrib_dic["transcript_id"]
                proteinid = attrib_dic.get("protein_id", None)
                start_stop_dt[transid]["start_codon"] = (int(left), int(right))


            elif feature == "stop_codon":
                transid = attrib_dic["transcript_id"]
                proteinid = attrib_dic.get("protein_id", None)
                start_stop_dt[transid]["stop_codon"] = (int(left), int(right))

            else:
                print(f"{feature} not expected by the utility", file=sys.stderr)


            self.gene = gene_dt
            self.transcript = transcript_dt
            self.exon = exon_dt
            self.CDS = cds_dt
            self.start_stop = start_stop_dt
    

    def get_gene(self):
        return self.gene_id, self.gene


    def get_transcript(self):
        return self.gene_id, self.transcript


    def get_exon(self):
        return self.gene_id, self.exon


    def get_CSD(self):
        return self.gene_id, self.CDS


    def get_start_stop(self):
        return self.gene_id, self.start_stop


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
        for line in dt:
            if line[0] == "#":
                meta_line_count += 1
                meta.append(line)
            else:
                break
        self.meta = meta
        dt = dt[meta_line_count:]
        if dt[-1][0] == "#":
            dt = dt[:-1]
        
        gtf_gene_data = {}
        index = 0
        key = None
        for line in dt:
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
    pass