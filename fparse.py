import copy
import re
import argparse
class Fasta_parse:
    """
    this class was writed to parse fasta file. 

    And there are three methods within class. fist is __inint__, to add data and seq_heads attibution to Instance. The second is fasta_print, to print self.data back to fasta file, and the length of fasta file line can be modifed by argrment 'line_len', dauflt length is 80 characters. The third is join_lines, which will modified self.data, and make value into a line

    And Instance contain tow attibution. one is data, it is a dictionary contain the information of fasta file, and the key is head line, and the value is the seq. second attribution is seq_heads this is a list, contain all head line of fasta file and by the fasta file order.
    """

    def __init__(self,file_name):
        def check_file_format(file_name):
            all_first_char=[line[0] for line in open(file_name)]
            all_char=[char for char in all_first_char if char != '>']
            #string=''.join(all_char)

            i=0
            if len(all_first_char)>0 and all_first_char[0]=='>': pass
            else: i=1
            #if string.isalpha: pass
            #else: i=1

            if i:
                return False
            else:
                return True


        check_res=check_file_format(file_name)
        #print(check_res)
        if check_res:
            self.fileformat='fasta'

            data={}
            seq_heads=[]
            for line in open(file_name,'r'):
                line=line.rstrip()
                if line[0]=='>':
                    key=line[1:]
                    seq_heads.append(key)
                    data[key]=[]
                else:
                    data[key].append(line)
            self.data=data
            self.seq_heads=seq_heads
            self.print_len=80
            self.data_struc='{head:[seq,...]}'

        else:
            print('The file may not be a fasta file.')
            exit()


    def join_lines(self):
        for key in self.data:
            self.data[key]=''.join(self.data[key])


    def fasta_print(self,line_len=80):
        for key in self.seq_heads:
            print('>'+key)
            whole_line=''.join(self.data[key])

            i=0
            part=whole_line[i:i+line_len]
            while part:
                print(part)
                i+=line_len
                part=whole_line[i:i+line_len]
   
        self.print_len=line_len
#-----------------------------------------------------------------


class Gff_parse:
    pass
#-----------------------------------------------------------------



class Gbk_parse:
    pass
#-----------------------------------------------------------------



class Blast_parse:
    """
        This is a Class for Blast result. the out format of blast is 0.
    """
    def __init__(self,file_name):
        def struc_data(f_name):
            f_in=open(f_name,'r')
            data_out={}

            matadata=[]
            query_res=[]
            i=0
            begin=re.compile(r'Query=\s.+')
            end=re.compile(r'\s{2}Database:\s.+')
            for line in f_in:
                line=line.rstrip()
                if len(line)<1: continue
                if begin.match(line): i=1
                if end.match(line): i=0

                if i==0:
                    matadata.append(line)
                if i==1:
                    query_res.append(line)
            f_in.close()
            data_out['mata_data']=matadata

            #split query records.
            query_records={}
            for line in query_res:
                if begin.match(line):
                    key=re.match(r'Query=\s(.+)',line).groups()[0]
                    if key not in query_records: query_records[key]=[]
                query_records[key].append(line)

            #split hits.
            length_re=re.compile(r'Length=(\d+)')
            hit_re=re.compile(r'>\s(.+)')
            mata_re=re.compile(r'Lambda.+')
            hit_end=re.compile(r'Effective\ssearch.+')
            for query in query_records:
                i=j=0
                length=None
                hit={}
                mata_dt=[]
                for line in query_records[query]:
                    if not length and length_re.match(line):
                        length=length_re.match(line).groups()[0]
                    if mata_re.match(line) or hit_end.match(line):i=1
                    if hit_re.match(line):
                        key=hit_re.match(line).groups()[0]
                        if not key in hit: hit[key]=[]
                        j=1
                    if i==1:
                        mata_dt.append(line)
                    elif j==1:
                        hit[key].append(line) 
                    else:
                        continue
                query_records[query]={}
                query_records[query]['mata_data']=mata_dt
                query_records[query]['query_len']=length
                query_records[query]['subjects']=hit

            #split hits match segments.
            score_re=re.compile(r'\sScore\s=\s(.+?),')
            for query in query_records:
                for hit in query_records[query]['subjects']:
                    hit_len=None
                    segment={}
                    i=1
                    j=0
                    for line in query_records[query]['subjects'][hit]:
                        #print(line)
                        if not hit_len and length_re.match(line):
                            hit_len=length_re.match(line).groups()[0]

                        if score_re.match(line):
                            key='segment_'+str(i)
                            if key not in segment: segment[key]=[]
                            i+=1
                            j=1
                        if j==1:
                            segment[key].append(line)

                    query_records[query]['subjects'][hit]={}
                    query_records[query]['subjects'][hit]['subject_len']=hit_len
                    query_records[query]['subjects'][hit]['match_segment']=segment


            def struc_segment(a_line):
                q_range=[]
                s_range=[]
                q=''
                s=''
                for line in a_line:
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
                return q_range,s_range,q,s

                    
            #structure hits match segments informations.
            query_ali=re.compile(r'Query\s{2}.+')
            expect_re=re.compile(r'.+Expect\s=\s(.+)')
            identities_re=re.compile(r'\sIdentities\s=\s(.+?),')
            gaps_re=re.compile(r'.+Gaps\s=\s(.+)')
            strand_re=re.compile(r'\sStrand=(.+)')
            for query in query_records:
                for hit in query_records[query]['subjects']:
                    for segment in query_records[query]['subjects'][hit]['match_segment']:
                        segment_info={}
                        query_sbjct=[]
                        i=0
                        for line in query_records[query]['subjects'][hit]['match_segment'][segment]:
                            if score_re.match(line):
                                score=score_re.match(line).groups()[0]
                                score=score.split()
                                score=[score[0],score[2].strip('()')]
                            if expect_re.match(line):
                                expect=expect_re.match(line).groups()[0]
                            if identities_re.match(line):
                                identities=identities_re.match(line).groups()[0]
                                identities=identities.split()
                                identities=[identities[0],identities[1].strip('()')]
                            if gaps_re.match(line):
                                gaps=gaps_re.match(line).groups()[0]
                                gaps=gaps.split()
                                gaps=[gaps[0],gaps[1].strip('()')]
                            if strand_re.match(line):
                                strand=strand_re.match(line).groups()[0]
                                strand=strand.split('/')
                            if query_ali.match(line):i=1
                            if i==1:query_sbjct.append(line)
                        query_range,sbjct_range,query_seq,sbjct_seq=struc_segment(query_sbjct)
                        segment_info['score']=score
                        segment_info['expect']=expect
                        segment_info['identities']=identities
                        segment_info['gaps']=gaps
                        segment_info['strand']=strand
                        segment_info['query_range']=query_range
                        segment_info['sbjct_range']=sbjct_range
                        segment_info['query_seq']=query_seq
                        segment_info['sbjct_seq']=sbjct_seq

                        query_records[query]['subjects'][hit]['match_segment'][segment]=segment_info

            data_out['query_records']=query_records

            return data_out


        self.data=struc_data(file_name)
        self.data_struc="{'mata_data':[str,...],'query_records':{query_name:{'mata_data':[],'query_len':str,'subjects':{subject_name:{'subject_len':str,''match_segment:{segment:{'score':[str,str],'expect':str,'identities':[str,str],'gaps':[str,str],'strand':[str,str],'query_range':[str,str],'sbjct_range':[str,str],'query_seq':str,'sbjct_seq':str},...}}}},...}}"


    def give_result_list(self):
        data_list=[]
        for query in self.data['query_records']:
            query_len=self.data['query_records'][query]['query_len']
            for sbjct in self.data['query_records'][query]['subjects']:
                sbjct_len=self.data['query_records'][query]['subjects'][sbjct]['subject_len']
                for segment in self.data['query_records'][query]['subjects'][sbjct]['match_segment']:
                    match_len,alig_len=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['identities'][0].split('/')
                    identities=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['identities'][1]
                    gaps=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['gaps'][0].split('/')[0]
                    score=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['score'][0].split('/')[0]
                    expect=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['expect']
                    strand='_'.join(self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['strand'])
                    q_start,q_end=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['query_range']
                    s_start,s_end=self.data['query_records'][query]['subjects'][sbjct]['match_segment'][segment]['sbjct_range']
                    mismatch=int(alig_len)-int(match_len)-int(gaps)

                    data_list.append([query,query_len,sbjct,sbjct_len,segment,alig_len,match_len,mismatch,gaps,strand,q_start,q_end,s_start,s_end,identities,expect,score])
                    #print(query,query_len,sbjct,sbjct_len,segment,alig_len,match_len,mismatch,gaps,strand,q_start,q_end,s_start,s_end,identities,expect,score)
        #order result list.
        order_dic={}
        sbjct_score={}
        for line in data_list:
            if line[2] not in order_dic:
                order_dic[line[2]]={}
            order_dic[line[2]][line[4]]=line
            if line[2] not in sbjct_score:
                sbjct_score[line[2]]=[]
            sbjct_score[line[2]].append(int(line[-1]))
        tmp=[]
        for sbjct in sbjct_score:
            tmp.append([sbjct,max(sbjct_score[sbjct])])
        tmp.sort(reverse=True,key=lambda x:x[1])
        tmp=[ele[0] for ele in tmp]
        #print(tmp)
        data_out=[]
        for key in tmp:
            contain=[]
            for segment in order_dic[key]:
                contain.append(order_dic[key][segment])
            contain.sort(reverse=True,key=lambda x:int(x[-1]))
            for ele in contain:
                data_out.append(ele)

        return data_out




#-----------------------------------------------------------------




#=================================================================
if __name__ == '__main__':
    import sys
    
    fdt=Blast_parse(sys.argv[1])

#    for query in fdt.data['query_records']:
#        print(query)
#        print(fdt.data['query_records'][query]['query_len'])
#        for hit in fdt.data['query_records'][query]['subjects']:
#            print(hit)
#            for segment in fdt.data['query_records'][query]['subjects'][hit]['match_segment']:
#                print(segment)
#                for key in fdt.data['query_records'][query]['subjects'][hit]['match_segment'][segment]: 
#                    print(key,fdt.data['query_records'][query]['subjects'][hit]['match_segment'][segment][key])

    #print(fdt.data)
    dlist=fdt.give_result_list()
    for ele in dlist:
        print(ele)





