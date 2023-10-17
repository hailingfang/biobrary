import os
import sys

class Fasta_block:
    """
        Data block class
    """
    def __init__(self, data):
        self.head = None
        self.head_line = None
        self.content = []
        data = data.rstrip().split("\n")
        self.head_line = data[0]
        self.data_lines = data[1:]
        self.head_id = self.head_line[1:].split()[0]
        
    def parse_headline(self):
        self.attributes = {}
    
    def join_data_lines(self):
        pass


class FASTA_parser:
    """
        This class was writed to parse fasta file.
    """
    def __init__(self, fasta_file):
        if not os.path.exists(fasta_file):
            print(f"error, file {fasta_file} not found", file=sys.stderr)
            exit()
        self.file_name = fasta_file
        self.head_seek = {}
        self.head_index = {}
        self.index_head = {}
        self.print_width = 80

        fin = open(fasta_file, "r" )
        head_index = 0
        while True:
            line = fin.readline()
            if not line:
                break
            if line[0] == ">":
                head = line[1:].rstrip().split()[0]
                self.head_seek[head] = fin.tell() - len(line) 
                self.head_index[head] = head_index
                self.index_head[head_index] = head
                head_index += 1
        fin.close()
        self.file_handle = open(fasta_file, "r")


    def read(self, head_id=None):
        print(self.head_seek[head_id])
        print(self.head_index[head_id])
        print(self.index_head[self.head_index[head_id]])
        block_data = None
        if head_id:
            seek_pos = self.head_seek.get(head_id, "NA")
            head_index = self.head_index.get(head_id, "NA")
            if (seek_pos == "NA") or (head_index == "NA"):
                print(f"error, head id {head_id} not found", file=sys.stderr)
                exit()
            fin = open(self.file_name, "r")
            fin.seek(seek_pos, os.SEEK_SET)

            head_index_next = head_index + 1
            read_len = None
            if head_index_next in self.index_head:
                head_next = self.index_head[head_index_next]
                head_next_seek_pos = self.head_seek.get(head_next, None)
                if head_next_seek_pos:
                    read_len = head_next_seek_pos - seek_pos
            if read_len:
                block_data = fin.read(read_len)
            else:
                block_data = fin.read()
            print(block_data) 
        else:
            head_line = self.file_handle.readline()
            if head_line[0] != ">":
                print("error, the first line should be head line.", file=sys.stderr)
                exit()
            block_data = None
        return  Fasta_block(block_data)


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


if __name__ == "__main__":
    '''
        For test purpose.
    '''
    import sys
    fasta_ins = Fasta_parser(sys.argv[1])
    #block_ins = fasta_ins.read("XR_943749.3")
    #block_ins = fasta_ins.read("NM_000014.6")
    block_ins = fasta_ins.read("XR_953308.2")
    print(block_ins.head_line)
    print(block_ins.data_lines)
