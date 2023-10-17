class GFF_parser:
    """
    This class was used to parse for gff file version3.
    Thist class just offer a data structure for GFF file.
    """
    def __init__(self, file_name):

        def struc_data(file_name):
            data_out = {'meta_dt': [], 'information': [], 'fasta_seq': []}
            seq_tmp_re = re.compile('##FASTA')
            i = 0
            for line in open(file_name):
                line = line.rstrip()
                if line[0] == '#':
                    data_out['meta_dt'].append(line)
                    if seq_tmp_re.match(line):
                        i = 1
                elif i == 0:
                    data_out['information'].append(line)
                elif i == 1:
                    data_out['fasta_seq'].append(line)
                else:
                    raise Exception('opps')
            tmp = {}
            for line in data_out['information']:
                line_ele = line.split('\t')
                contig = line_ele[0]
                source = line_ele[1]
                seq_type = line_ele[2]
                position = (int(line_ele[3]), int(line_ele[4]))
                filde_6 = line_ele[5]
                strand = line_ele[6]
                phase = line_ele[7]
                attributes = line_ele[8].split(';')
                attributes_tmp = {}
                for ele in attributes:
                    ele = ele.split('=')
                    attributes_tmp[ele[0]] = ele[1]
                attributes = attributes_tmp
                if contig not in tmp:
                    tmp[contig] = {}
                if source not in tmp[contig]:
                    tmp[contig][source] = {}
                if seq_type not in tmp[contig][source]:
                    tmp[contig][source][seq_type] = {}
                tmp[contig][source][seq_type][position] = {'score':filde_6, \
                    'strand':strand, 'phase':phase, 'attributes':attributes}
            data_out['information'] = tmp
            return data_out

        self.data = struc_data(file_name)
        self.data_struc="{'mata_dt':[str], \
            'information':{contig:{source:{seq_type:{position:{'score':str, \
            'strand':str, 'phase':str, attributes:{...}}}}}}, \
            'seq':[]}"


