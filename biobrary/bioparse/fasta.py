"""
Founctions and Classes to parse FASTA file.

The main entry to parse FASTA file is parse_fasta. The function read
a fasta file and reture a FASTA class. The FASTA class hold FASTA_ENTRY
objects. And the information of FASTA_ENTRY objects can be inquired
by method of FASTA_ENTRY class.
"""

import gzip


class FASTA_ENTRY:
    """
    The FASTA_ENTRY class, storing the information of entries and inquiry information
    from the entries.
    """
    def __init__(self, seq_id, head_line, entry_start, entry_size, fasta_file):
        self._seq_id = seq_id
        self._head_line = head_line
        self._entry_start = entry_start
        self._entry_size = entry_size
        self._file = fasta_file
        self._seq_loaded = False
        self._seq = None
        self._seq_len = None

    def __str__(self):
        return self.get_seq_id()

    def _load_seq(self):
        if not self._seq_loaded:
            if self._file.endswith(".gz"):
                fin = gzip.open(self._file, "r")
                fin.seek(self._entry_start, 0)
                self._seq = ''.join(fin.read(self._entry_size).decode().split('\n')[:-1])
                self._seq_loaded = True
                fin.close()
            else:
                fin = open(self._file, "r")
                fin.seek(self._entry_start, 0)
                self._seq = ''.join(fin.read(self._entry_size).split('\n')[:-1])
                self._seq_loaded = True
                fin.close()

    def get_seq_id(self):
        """
        Get the sequence id.

        Returns
        -----------
        seq_id: string, the seq_id of the entry.
        """
        return self._seq_id

    def get_head_line(self):
        """
        Get the head line of entry.

        Returns
        -----------------
        head_line: string, head line of the entry.
        """
        return self._head_line

    def get_seq(self):
        """
        Get the sequence of the entry.

        Returns
        ---------------
        seq: string, the sequence of the entry.
        """
        if self._seq_loaded:
            return self._seq
        else:
            self._load_seq()
            return self._seq

    def format_entry_str(self, width=80):
        """
        To format the entry head lines and sequence for printing.

        Returns
        --------------
        string: the formated string for printing.
        """
        lines = []
        block_idx = list(range(0, self._seq_len, width))
        block_idx.append(self._seq_len)
        lines.append(self._head_line)
        for idx in range(len(block_idx) - 1):
            lines.append(self._seq[block_idx[idx]: block_idx[idx + 1]])
        return "\n".join(lines)


class FASTA:
    """
    The FASTA class to parse FASTA file.
    """
    def __init__(self):
        self._seq_id_entry_dic = {}
        self._seq_id_list = []

    def __iter__(self):
        self._start = 0
        self._stop = len(self._seq_id_list)
        return self

    def __next__(self):
        if self._start < self._stop:
            self._start += 1
            return self._seq_id_entry_dic[self._seq_id_list[self._start - 1]]
        else:
            raise StopIteration

    def get_seq_id_s(self):
        """
        Get seq_id_s of all entries.

        Returns
        --------------
        seq_list: a list which contain all sequence ids of entries.
        """
        return self._seq_id_list

    def get_seq_entry(self, seq_id):
        """
        Get entry by seq_id.

        Parameters
        ----------------
        seq_id: the id of the entry.

        Returns
        ----------------
        FASTA entry: the FASTA ENTRY object inquired by the id.
        """
        return self._seq_id_entry_dic.get(seq_id)


def _read_fasta(fasta_file):
    """
    To read the FASTA file, and split the file contents to entries.
    """
    data = {}
    head_len = []
    head_line_s = []
    file_pos_s = []
    if fasta_file.endswith(".gz"):
        fin = gzip.open(fasta_file, "r")
        while True:
            line = fin.readline()
            if line:
                if line[0] == 62:
                    line = line.decode()
                    head_len.append(len(line))
                    head_line_s.append(line.rstrip())
                    file_pos_s.append(fin.tell())
            else:
                break
        file_pos_s.append(fin.tell())
        head_len.append(0)
        head_len = head_len[1:]
        fin.close()
    else:
        fin = open(fasta_file, "r")
        while True:
            line = fin.readline()
            if line:
                if line[0] == '>':
                    head_len.append(len(line))
                    head_line_s.append(line.rstrip())
                    file_pos = fin.tell()
                    file_pos_s.append(file_pos)
            else:
                break
        file_pos_s.append(fin.tell())
        head_len.append(0)
        head_len = head_len[1:]
        fin.close()
    for idx, head_line in enumerate(head_line_s):
        #sequence data start position and its offset.
        data[head_line] = [file_pos_s[idx], file_pos_s[idx + 1] - file_pos_s[idx] - head_len[idx]] 
    return data


def _default_head_parse_fun(head_line):
    """
    The default function to parse head lines of FASTA entry.
    """
    return head_line.split()[0][1:]


def parse_fasta(fasta_file, load_seq=False, head_parse_func=_default_head_parse_fun):
    """
    Read FASTA file, and return a FASTA object.

    Parameters
    -------------
    fasta_file: string, the name of FASTA file, plain text or gz file.
    load_seq: bool, default=False, Load the sequence into memory when create FASTA object.
    head_parse_func: fuction, default=None, The function used to parse head line,
        the function should return seq_id.

    Returns
    ------------
    FASTA_object: An instance of FASTA class.
    """
    fasta_data = _read_fasta(fasta_file)
    fasta = FASTA()

    for head_line in fasta_data:
        seq_id = head_parse_func(head_line)
        entry_start, entry_size = fasta_data[head_line]
        fasta_entry = FASTA_ENTRY(seq_id, head_line, entry_start, entry_size, fasta_file)
        fasta._seq_id_entry_dic[seq_id] = fasta_entry
        fasta._seq_id_list.append(seq_id)
    if load_seq:
        if fasta_file.endswith('gz'):
            dt = {}
            seq_id = None
            fin = gzip.open(fasta_file, 'r')
            for line in fin:
                line = line.decode().rstrip()
                if line[0] == '>':
                    seq_id = head_parse_func(line)
                    dt[seq_id] = []
                else:
                    dt[seq_id].append(line)
            fin.close()
            for key in dt:
                dt[key] = ''.join(dt[key])
            for seq_id in fasta.get_seq_id_s():
                entry = fasta.get_seq_entry(seq_id)
                seq = dt[seq_id]
                entry._seq = seq
                entry._seq_len = len(seq)
                entry._seq_loaded = True
        else:
            dt = {}
            seq_id = None
            fin = open(fasta_file, 'r')
            for line in fin:
                line = line.rstrip()
                if line[0] == '>':
                    seq_id = head_parse_func(line)
                    dt[seq_id] = []
                else:
                    dt[seq_id].append(line)
            fin.close()
            for key in dt:
                dt[key] = ''.join(dt[key])
            for seq_id in fasta.get_seq_id_s():
                entry = fasta.get_seq_entry(seq_id)
                seq = dt[seq_id]
                entry._seq = seq
                entry._seq_len = len(seq)
                entry._seq_loaded = True
    return fasta


def test_fasta(fasta_file):
    """
    The testing function. Test the read of the FASTA file.

    Parameters
    ---------------
    fasta_file: the file name of FASTA file.
    """
    fasta = parse_fasta(fasta_file, load_seq=True)
    for ent in fasta:
        print(ent.format_entry_str())
        break


if __name__ == "__main__":
    import sys
    test_fasta(sys.argv[1])
    
