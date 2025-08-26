import logging
import gzip
import logging.handlers
from datetime import datetime
from pathlib import Path


###########################
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"16s_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('16S_Analysis')
logger.setLevel(logging.INFO)



def read_fasta(filepath):
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"the input file {filepath} doesn't exist!")
    opener = gzip.open if filepath.suffix == '.gz' else open
    with opener(filepath, 'rt') as f:
        header, seq = '', ''
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield (header, seq)
                header, seq = line, ''
            else:
                seq += line
        if header:
            yield (header, seq)


def read_fastq(filepath):
    filepath = str(filepath)
    opener = gzip.open if filepath.endswith('.gz') else open
    with opener(filepath, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            if not header.startswith('@'):
                raise ValueError(f"Invalid FASTQ header: '{header}' (must start with '@')")
            seq = f.readline().strip()
            _ = f.readline().strip()
            qual = f.readline().strip()
            yield (header, seq, qual)


def save_files(outfile, reads):
    with open(outfile, 'w') as f:
        for read in reads:
            # fastq format
            if len(read) == 3:
                header, seq, qual = read
                f.write(f"{header}\n{seq}\n+\n{qual}\n")
            # fasta format
            else:
                header, seq = read
                if header.startswith('@'):
                    header = '>' + header[1:]
                f.write(f"{header}\n{seq}\n")


# generate kmers
def yield_kmers(seq, k):
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if 'N' not in kmer:
            yield kmer
