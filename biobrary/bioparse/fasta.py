"""
Founctions and Class to parse FASTA file.
"""

import gzip


class FASTA_ENTRY:
    def __init__(self, seq_id, seq_info, entry_start, entry_size, fasta_file):
        self._seq_id = seq_id
        self._seq_info = seq_info
        self._entry_start = entry_start
        self._entry_size = entry_size
        self._file = fasta_file
        self._seq_loaded = False
        self._seq = None

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
        return self._seq_id

    def get_seq_info(self):
        return self._seq_info

    def get_seq(self):
        if self._seq_loaded:
            return self._seq
        else:
            self._load_seq()
            return self._seq


class FASTA:
    def __init__(self):
        self.seq_id_entry_dic = {}

    def get_seq_id_s(self):
        return list(self.seq_id_entry_dic.keys())

    def get_seq_entry(self, seq_id):
        return self.seq_id_entry_dic[seq_id]


def _read_fasta(fasta_file):
    data = {}
    file_pos_s = []
    head_len = []
    seq_id_s = []
    seq_info_s = []
    if fasta_file.endswith(".gz"):
        fin = gzip.open(fasta_file, "r")
        while True:
            line = fin.readline()
            if line:
                if line[0] == 62:
                    line = line.decode()
                    head_len.append(len(line))
                    line = line.rstrip().split()
                    seq_id = line[0][1:]
                    seq_info = ' '.join(line[1:])
                    file_pos = fin.tell()
                    seq_id_s.append(seq_id)
                    seq_info_s.append(seq_info)
                    file_pos_s.append(file_pos)
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
                    line = line.rstrip().split()
                    seq_id = line[0][1:]
                    seq_info = ' '.join(line[1:])
                    file_pos = fin.tell()
                    seq_id_s.append(seq_id)
                    seq_info_s.append(seq_info)
                    file_pos_s.append(file_pos)
            else:
                break
        file_pos_s.append(fin.tell())
        head_len.append(0)
        head_len = head_len[1:]
        fin.close()
    for idx, seq_id in enumerate(seq_id_s):
        data[seq_id] = [seq_info_s[idx], file_pos_s[idx], file_pos_s[idx + 1] - file_pos_s[idx] - head_len[idx]]
    return data


def parse_fasta(fasta_file, load_seq=False):
    fasta_data = _read_fasta(fasta_file)
    fasta = FASTA()
    for seq_id in fasta_data:
        seq_info, entry_start, entry_size = fasta_data[seq_id]
        fasta_entry = FASTA_ENTRY(seq_id, seq_info, entry_start, entry_size, fasta_file)
        fasta.seq_id_entry_dic[seq_id] = fasta_entry
    if load_seq:
        if fasta_file.endswith('gz'):
            dt = {}
            seq_id = None
            fin = gzip.open(fasta_file, 'r')
            for line in fin:
                line = line.decode().rstrip()
                if line[0] == '>':
                    seq_id = line[1:].split()[0]
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
                entry._seq_loaded = True
        else:
            dt = {}
            seq_id = None
            fin = open(fasta_file, 'r')
            for line in fin:
                line = line.rstrip()
                if line[0] == '>':
                    seq_id = line[1:].split()[0]
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
                entry._seq_loaded = True
    return fasta


def test_fasta(fasta_file):
    fasta = parse_fasta(fasta_file, load_seq=True)
    for seq_id in fasta.get_seq_id_s():
        entry = fasta.get_seq_entry(seq_id)
        seq = entry.get_seq()
        print(seq_id, len(seq))


if __name__ == "__main__":
    import sys
    test_fasta(sys.argv[1])
    
