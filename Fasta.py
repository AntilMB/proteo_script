class Fasta():
    def __init__(self, file=None, only_first=False):
        self.data = {}        
        if file is None:
            return
        
        self.read_fasta(file, only_first) 
        
    def get_key(self, stirng):
        return stirng[1:].split('|')[0].upper()
        
    def read_fasta(self, files, only_first=False):
        x = open(files)
        key = None
        seq = ""
        for line in x:
            line = line.strip()
            if key is None and line.startswith(">"):
                key = line[1:] if not only_first else self.get_key(line)
            elif not(key is None) and line.startswith(">"):
                self.data[key] = seq
                seq = ""
                key = line[1:] if not only_first else self.get_key(line)
            else:
                seq += line.upper()
        self.data[key] = seq
