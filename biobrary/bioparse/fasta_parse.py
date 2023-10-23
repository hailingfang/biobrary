import os
import sys


class Fasta_seq:
    """
    Fasta data block class.
    """
    def __init__(self, seqid, seqdata, seqid_append=None, seqid_make_func=None):
        self.seqid = seqid
        self.seqdata = seqdata
        self.seqid_append = seqid_append
        self.seqlen = len(self.seqdata)
        self.seqid_make_func = seqid_make_func
    

    def __str__(self):
        return ">" + self.seqid + "\n" + self.seqdata


    def substr(self, left, right):
        """return substr of seqdata"""

        if left > right or left > self.seqlen or right > self.seqlen:
            print("Wraning, illeagl left or right border", file=sys.stderr)            
        return self.seqdata[left - 1: right]
    
    def get_seqid(self):
        return self.seqid

    def print(self, width=80, seqid_make_func=None, file=sys.stdout):
        """print fastq seqdata"""

        if not seqid_make_func and self.seqid_make_func:
            seqid_make_func = self.seqid_make_func
        if seqid_make_func:
            seqid = seqid_make_func(self.seqid, self.seqid_append)
            print(">" + seqid, file=file)
            for i in range(self.seqlen)[::width]:
                print(self.seqdata[i: i + width], file=file)
        else:
            print(">" + self.seqid, file=file)
            for i in range(self.seqlen)[::width]:
                print(self.seqdata[i: i + width], file=file)




class FASTA:
    """
    This class was writed to parse fasta file.
    """
    def __init__(self, fasta_file, seqid_format_func=None):
        if not os.path.exists(fasta_file):
            print(f"error, file {fasta_file} not found", file=sys.stderr)
            exit()

        self.file_name = fasta_file
        self.indexed = False
        self.seqid_order = []
        self.seqid_info = {}
        self.seqid_seek = {}
        self.seqid_index = {}
        self.index_seqid = {}
        self.seqid_format_func = seqid_format_func
        self.index_file()


    def __iter__(self):
        self.start = 0
        self.stop = len(self.seqid_order)
        return self


    def __next__(self):
        if self.start == self.stop:
            raise StopIteration
        self.start += 1
        return self.get_seq(self.seqid_order[self.start - 1])


    def index_file(self):
        """
        Index file seek position for every seqid.
        """
        if not self.indexed:
            seqid_format_func = self.seqid_format_func
            fin = open(self.file_name, "r" )
            index = 0
            while True:
                line = fin.readline()
                if not line:
                    break
                if line[0] == ">":
                    line = line.rstrip()
                    seqid_raw = line[1:].split()[0]
                    if seqid_format_func:
                        seqid = seqid_format_func(seqid_raw)
                    else:
                        seqid = seqid_raw

                    self.seqid_order.append(seqid)
                    if seqid in self.seqid_info:
                        print(f"Warning, seqid {seqid} duplicated", file=sys.stderr)
                    else:
                        self.seqid_info[seqid] = line[1:]

                    self.seqid_seek[seqid] = fin.tell() - len(line) - 1
                    self.seqid_index[seqid] = index
                    self.index_seqid[index] = seqid
                    index += 1
            fin.close()
            self.indexed = True


    def get_seqids(self):
        """reture all seqid by order along with the fasta file"""
        return self.seqid_order
    
    
    def get_seqid_info(self, seqid):
        """return seq information of seqid"""

        return self.seqid_info[seqid]
    

    def get_seq(self, seqid):
        """
        Make fasta seq instance

        Parameters
        --------------
        seqid : a string
            fastq sequence id.

        Returns
        ------------
        Fasta_data : A instance of Fasta data
        """
        dt_seek_pos = self.seqid_seek.get(seqid)
        if dt_seek_pos != 0 and not dt_seek_pos:
            print("sequence id not found ...", file=sys.stderr)
            exit()
        next_seqid = self.index_seqid.get(self.seqid_index[seqid] + 1)
        next_dt_seek_pos = self.seqid_seek.get(next_seqid)

        fin = open(self.file_name, "r")
        fin.seek(dt_seek_pos)
        if next_dt_seek_pos:
            seq_data = fin.read(next_dt_seek_pos - dt_seek_pos)
        else:
            seq_data = fin.read()
        fin.close()
        seq_data = seq_data.split("\n")
        seq_seqid_line = " ".join(seq_data[0].split()[1:])
        seq_data = "".join(seq_data[1:])

        return Fasta_seq(seqid, seq_data, seqid_append=seq_seqid_line)




if __name__ == "__main__":
    """For test purpose."""

    fasta = FASTA(sys.argv[1])
    
