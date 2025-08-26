import time
import gzip
import logging
import sqlite3
from pathlib import Path
from config import db_config

###############
# this script is used for constructing kmer database by using reference data
###############

start = time.time()

### ====== logs ====== ###
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.FileHandler("../kmer_db_builder.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


### ====== base functions ====== ###
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


### ====== main functions ====== ###
# generate kmers
def yield_kmers(seq, k):
    seq = seq.upper()
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i:i + k]
        if 'N' not in kmer:
            yield kmer


def process_sequence(seq, source, label, k_list, cursor):
    if seq and source and label:
        for k in k_list:
            kmers = set(yield_kmers(seq, k))
            if kmers:
                table_name = f"kmer_labels_k{k}"
                cursor.executemany(
                    f'INSERT OR IGNORE INTO {table_name} VALUES (?,?,?)',
                    [(kmer, label, source) for kmer in kmers])


def create_table(cursor, k_list):
    for k in k_list:
        table_name = f"kmer_labels_k{k}"
        cursor.execute(f"""
            CREATE TABLE IF NOT EXISTS {table_name}(
                kmer TEXT NOT NULL,
                label TEXT NOT NULL, 
                source TEXT,
                PRIMARY KEY (kmer, label))
        """)
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_kmer_{k} ON {table_name}(kmer)")


# construct kmer SQLite database, easy to query
def build_kmer_db(ref_fasta, kmer_db, k_list):
    with sqlite3.connect(kmer_db) as conn:
        cursor = conn.cursor()
        create_table(cursor, k_list)

        if ref_fasta:
            fasta_path = Path(ref_fasta)
            for header, seq in read_fasta(fasta_path):
                parts = header.split()
                if len(parts) >= 3:
                    source = parts[0]
                    label = ' '.join(parts[2:])
                    process_sequence(seq, source, label, k_list, cursor)
                else:
                    logger.error(f"Invalid header skipped: {header}")
        conn.commit()

        for k in k_list:
            cursor.execute(f"SELECT COUNT(DISTINCT kmer), COUNT(DISTINCT label) FROM kmer_labels_k{k}")
            kmer_count, label_count = cursor.fetchone()
            logger.info(f"k={k}: {kmer_count:,} kmers -> {label_count:,} labels")

    logger.info(f"Database constructed: {kmer_db}")


def main(ref_fasta, kmer_db, k_list):
    logger.info("Building k-mer database ...")
    build_kmer_db(ref_fasta, kmer_db, k_list)


if __name__ == '__main__':
    main(**db_config)
    logger.info(f"All kmer tables created in {time.time() - start:.2f}s")
