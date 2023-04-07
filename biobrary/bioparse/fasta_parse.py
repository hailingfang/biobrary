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


