import re


class Fasta_parser:
    """
    This class was writed to parse fasta file.
    """
    def __init__(self, file_name):
        self.data = {}
        self.heads = []
        with open(file_name) as f_in:
            for line in f_in:
                line = line.strip()
                if line[0] == '>':
                    head = line[1:]
                    self.heads.append(head)
                    self.data[head] = ''
                else:
                    self.data[head] += line
            f_in.close()
        self.print_width = 80

    def __str__(self):
        print_str = ''
        for head in self.data:
            print_str += '>' + head + '\n'
            for i in range(len(self.data[head]))[::self.print_width]:
                print_str += self.data[head][i: i + self.print_width] + '\n'
            print_str = print_str.rstrip()
        return print_str

    def set_print_width(self, width):
        self.print_width = width

    def __iter__(self):
        self.all_keys = self.heads
        self.all_keys.reverse()
        return self

    def __next__(self):
        if len(self.all_keys) == 0:
            raise StopIteration
        key = self.all_keys.pop()
        out = self.data[key]
        return [key, out]

    def trim_head(self):
        data_head_trimed = {}
        for head in self.data:
            data_head_trimed[head.split()[0]] = self.data[head]
        self.data = data_head_trimed


class Gff_parser:
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




#-------------------------------------------------------------------------------
class Gbk_parse:
    pass




#-------------------------------------------------------------------------------
class Blast_parse:
    """
        This is a Class for Blast result. the out format of blast is 0.
    """
    def __init__(self,file_name):

        def split_matadata_query_res(file_lines):
            matadata=[]
            query_res=[]
            i=0
            begin=re.compile(r'Query=\s.+')
            end=re.compile(r'\s{2}Database:\s.+')
            for line in file_lines:
                line=line.rstrip()
                if len(line)<1: continue
                if begin.match(line): i=1
                if end.match(line): i=0

                if i==0:
                    matadata.append(line)
                if i==1:
                    query_res.append(line)
            return matadata,query_res

        def split_query_records(query_res):
            query_records=[]
            query_begin=re.compile(r'Query=\s.+')
            tmp=[]
            for line in query_res:
                if query_begin.match(line):
                    query_records.append(tmp)
                    tmp=[]
                tmp.append(line)
            query_records.append(tmp)
            query_records=query_records[1:]
            return query_records

        def split_per_query(query_records):
            data_out={}
            query_re=re.compile(r'Query=.*')
            result_table_re_1=re.compile(r'Sequences\sproducing\ssignificant.*')
            result_table_re_2=re.compile(r'\*\*\*\*\*.*')
            query_mata_re_1=re.compile(r'Lambda.*')
            query_mata_re_2=re.compile(r'Effective\ssearch\sspace\sused:.*')
            query_name_re=re.compile(r'Query=\s(.+?)Length.*')
            query_length_re=re.compile(r'.*Length=(\d+).*')
            for query in query_records:
                head=[]
                body=[]
                tail=[]
                i=j=k=0
                for line in query:
                    if query_re.match(line):i=1
                    if result_table_re_1.match(line) or result_table_re_2.match(line):j=1
                    if query_mata_re_1.match(line) or query_mata_re_2.match(line):k=1
                    if i==1 and j==0 and k==0:
                        head.append(line)
                    if i==1 and j==1 and k==0:
                        body.append(line)
                    if i==1 and j==1 and k==1:
                        tail.append(line)
                head=''.join(head)
                #print(head)
                query_name=query_name_re.match(head).groups()[0]
                query_len=query_length_re.match(head).groups()[0]
                data_out[query_name]={'query_len':query_len,'matadata':tail,'body':body}
            return data_out

        def split_hits(query_records):
            significent_alig_re=re.compile(r'Sequences\sproducing\ssignificant\salignments:.+')
            hit_re=re.compile(r'>\s(.+)')
            for query in query_records:
                alig_summary=[]
                hits=[]
                i=j=0
                for line in query_records[query]['body']:
                    if i==0 and significent_alig_re.match(line):i=1
                    if j==0 and hit_re.match(line): j=1
                    if i==1 and j==0:
                        alig_summary.append(line)
                    if i==1 and j==1:
                        hits.append(line)
                query_records[query]['alig_summary']=alig_summary
                query_records[query].pop('body')
                tmp=[]
                block=[]
                for line in hits:
                    if hit_re.match(line):
                        tmp.append(block)
                        block=[]
                    block.append(line)
                tmp.append(block)
                hits=tmp[1:]
                query_records[query]['hits']=hits
            return query_records

        def split_hit_segment(query_records):
            head_re=re.compile(r'>\s.*')
            score_re=re.compile(r'\sScore\s=\s(.+?),')
            hit_name_re=re.compile(r'>\s(.+?)Length.*')
            hit_len_re=re.compile(r'.*Length=(\d+).*')
            for query in query_records:
                hits={}
                for hit in query_records[query]['hits']:
                    head=[]
                    body=[]
                    i=j=0
                    for line in hit:
                        if i==0 and head_re.match(line): i=1
                        if j==0 and score_re.match(line): j=1
                        if i==1 and j==0:
                            head.append(line)
                        if i==1 and j==1:
                            body.append(line)

                    head=''.join(head)
                    #print(head)
                    hit_name=hit_name_re.match(head).groups()[0]
                    hit_len=hit_len_re.match(head).groups()[0]

                    hits[hit_name]={'hit_len':hit_len}
                    segment={}
                    i=1
                    for line in body:
                        if score_re.match(line):
                            key='segment_'+str(i)
                            if key not in segment: segment[key]=[]
                            i+=1
                        segment[key].append(line)
                    hits[hit_name]['segments']=segment
                query_records[query]['hits']=hits
            return query_records

        def struc_query_subjct(qs_line):
            q_range=[]
            s_range=[]
            q=''
            s=''
            for line in qs_line:
                line=line.split()
                if len(line)==4:
                    if line[0]=='Query':
                        q_range+=[line[1],line[3]]
                        q+=line[2]
                    elif line[0]=='Sbjct':
                        s_range+=[line[1],line[3]]
                        s+=line[2]
            q_range=[q_range[0],q_range[-1]]
            s_range=[s_range[0],s_range[-1]]
            return [q_range,s_range,q,s]

        def extrect_alig_info(alig_info):
            alig_info_dic={}
            score_re=re.compile(r'.*?(Score)\s=\s(.+?)(?:,|$)')
            expect_re=re.compile(r'.*?(Expect)\s=\s(.+?)(?:,|$)')
            identities_re=re.compile(r'.*?(Identities)\s=\s(.+?)(?:,|$)')
            gap_re=re.compile(r'.*?(Gaps)\s=\s(.+?)(?:,|$)')
            positives_re=re.compile(r'.*?(Positives)\s=\s(.+?)(?:,|$)')
            strand_re=re.compile(r'.*?(Strand)=(.+?)(?:,|$)')
            split_identities_gaps_positives_re=re.compile(r'.*?(\d+)/(\d+)\s\((.+)\).*')
            split_score_re=re.compile(r'(.+?)\sbits\s')
            patterns=[score_re,expect_re,identities_re,gap_re,positives_re,strand_re]
            alig_len=None
            gap_num=None
            match_num=None
            mismatch_num=None
            for line in alig_info:
                for pattern in patterns:
                    re_res=pattern.match(line)
                    if re_res:
                        if re_res.groups()[0]=='Identities':
                            re_split=split_identities_gaps_positives_re.match(re_res.groups()[1])
                            #print(re_split)
                            alig_info_dic['Identities']=re_split.groups()[2]
                            match_num=re_split.groups()[0]
                            #print('>>>>>>>',match_num)
                            alig_len=re_split.groups()[1]
                        elif re_res.groups()[0]=='Gaps':
                            re_split=split_identities_gaps_positives_re.match(re_res.groups()[1])
                            #print(re_res.groups()[1])
                            #print(re_split)
                            alig_info_dic['Gaps']=re_split.groups()[2]
                            gap_num=re_split.groups()[0]
                        elif re_res.groups()[0]=='Positives':
                            re_split=split_identities_gaps_positives_re.match(re_res.groups()[1])
                            alig_info_dic['Positives']=re_split.groups()[2]
                        elif re_res.groups()[0]=='Score':
                            #print(re_res.groups()[1])
                            alig_info_dic['Score']=split_score_re.match(re_res.groups()[1]).groups()[0]
                        else:
                            alig_info_dic[re_res.groups()[0]]=re_res.groups()[1]
            if alig_len and gap_num and match_num:
                mismatch_num=str(int(alig_len)-int(match_num)-int(gap_num))
            alig_info_dic['alig_len']=alig_len
            alig_info_dic['gap_num']=gap_num
            alig_info_dic['match_num']=match_num
            alig_info_dic['mismatch_num']=mismatch_num
            return alig_info_dic

        def struc_segment(query_records):
            score_re=re.compile(r'\sScore\s=\s(.+?),')
            query_re=re.compile(r'Query.+')
            for query in query_records:
                for hit in query_records[query]['hits']:
                    for segment in query_records[query]['hits'][hit]['segments']:
                        qs_lines=[]
                        alig_info=[]
                        i=j=0
                        for line in query_records[query]['hits'][hit]['segments'][segment]:
                            if score_re.match(line): i=1;j=0
                            if query_re.match(line): i=0;j=1
                            if i==1:
                                alig_info.append(line)
                            elif j==1:
                                qs_lines.append(line)
                        query_subjct=struc_query_subjct(qs_lines)
                        alig_info=extrect_alig_info(alig_info)
                        query_records[query]['hits'][hit]['segments'][segment]={'alignment':query_subjct,'alig_info':alig_info}
            return query_records

        def struc_data(f_name):
            data_out={}
            f_in=open(f_name,'r')
            a_lines=f_in.readlines()
            query_matadata,query_res=split_matadata_query_res(a_lines)
            del a_lines
            data_out['mata_data']=query_matadata
            query_records=split_query_records(query_res)
            query_records=split_per_query(query_records)
            query_records=split_hits(query_records)
            query_records=split_hit_segment(query_records)
            query_records=struc_segment(query_records)
            data_out['query_records']=query_records
            return data_out

        self.data=struc_data(file_name)
        self.data_struc="{'mata_data':[str,...],'query_records':{query_name:{'mata_data':[],'query_len':str,'aligment_summary':[],'hits':{hit_name:{'hit_len':str,'segments':{segment:{'alignment':[q_range,s_range,q,s],'alig_info':{'Score'?:str,'Expect'?:str,...}},...}}}},...}}"


    def get_summary(self):
        data_list=[]
        for query in self.data['query_records']:
            query_name=query
            query_length=self.data['query_records'][query]['query_len']
            for hit in self.data['query_records'][query]['hits']:
                hit_name=hit
                hit_len=self.data['query_records'][query]['hits'][hit]['hit_len']
                for segment in self.data['query_records'][query]['hits'][hit]['segments']:
                    segment_name=segment
                    alig_len=self.data['query_records'][query]['hits'][hit]['segments'][segment]['alig_info']['alig_len']
                    match_num=self.data['query_records'][query]['hits'][hit]['segments'][segment]['alig_info']['match_num']
                    mismatch_num=self.data['query_records'][query]['hits'][hit]['segments'][segment]['alig_info']['mismatch_num']
                    gap_num=self.data['query_records'][query]['hits'][hit]['segments'][segment]['alig_info']['gap_num']
                    identities=self.data['query_records'][query]['hits'][hit]['segments'][segment]['alig_info']['Identities']
                    score=self.data['query_records'][query]['hits'][hit]['segments'][segment]['alig_info']['Score']
                    query_start,query_end=self.data['query_records'][query]['hits'][hit]['segments'][segment]['alignment'][0]
                    subjct_start,subjct_end=self.data['query_records'][query]['hits'][hit]['segments'][segment]['alignment'][1]
                    data_list.append([query_name,query_length,hit_name,hit_len,segment_name,alig_len,match_num,mismatch_num,gap_num,
                    query_start,query_end,subjct_start,subjct_end,identities,score])
        return data_list


#        for query in self.data['query_records']:
#            query_len=self.data['query_records'][query]['query_len']
#            for sbjct in self.data['query_records'][query]['subjects']:
#                sbjct_len=self.data['query_records'][query]['subjects'][sbjct]['subject_len']
#                for segment in self.data['query_records'][query]['subjects'][sbjct]['match_segment']:
#                    match_len,alig_len=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['identities'][0].split('/')
#                    identities=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['identities'][1]
#                    gaps=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['gaps'][0].split('/')[0]
#                    score=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['score'][0].split('/')[0]
#                    expect=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['expect']
#                    strand='_'.join(self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['strand'])
#                    q_start,q_end=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['query_range']
#                    s_start,s_end=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['sbjct_range']
#                    mismatch=int(alig_len)-int(match_len)-int(gaps)
#
#                    data_list.append([query,query_len,sbjct,sbjct_len,segment,alig_len,match_len,mismatch,gaps,strand,q_start,q_end,s_start,s_end,identities,expect,score])
#                    #print(query,query_len,sbjct,sbjct_len,segment,alig_len,match_len,mismatch,gaps,strand,q_start,q_end,s_start,s_end,identities,expect,score)
#        #order result list.
#        order_dic={}
#        sbjct_score={}
#        for line in data_list:
#            if line[2] not in order_dic:
#                order_dic[line[2]]={}
#            order_dic[line[2]][line[4]]=line
#            if line[2] not in sbjct_score:
#                sbjct_score[line[2]]=[]
#            sbjct_score[line[2]].append(int(line[-1]))
#        tmp=[]
#        for sbjct in sbjct_score:
#            tmp.append([sbjct,max(sbjct_score[sbjct])])
#        tmp.sort(reverse=True,key=lambda x:x[1])
#        tmp=[ele[0] for ele in tmp]
#        #print(tmp)
#        data_out=[]
#        for key in tmp:
#            contain=[]
#            for segment in order_dic[key]:
#                contain.append(order_dic[key][segment])
#            contain.sort(reverse=True,key=lambda x:int(x[-1]))
#            for ele in contain:
#                data_out.append(ele)
#
#        return data_out
#-------------------------------------------------------------------------------

#===============================================================================
if __name__ == '__main__':
    import sys
    print(sys.argv)
    dt = Gff_parser(sys.argv[1])
    for contig in dt.data['information']:
        print('>', contig)
        for source in dt.data['information'][contig]:
            print('>>', source)
            for seq_type in dt.data['information'][contig][source]:
                print('>>>', seq_type)
                i = 0
                for position in dt.data['information'][contig][source][seq_type]:
                    #print(position)
                    i += 1
                    if i <10:
                        print(dt.data['information'][contig][source][seq_type][position])
                    else:
                        break
