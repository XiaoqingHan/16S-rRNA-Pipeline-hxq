import time
from collections import Counter
from config import step3_config
from multiprocessing import Pool, cpu_count
from common_func import *
from db_utils import *
import json


### ====== base functions ====== ###
def build_kmer_cache(seqs, k):
    seq_kmers = {}
    kmer_counter = Counter()
    for header, seq in seqs:
        kmers = list(yield_kmers(seq, k))
        seq_kmers[(header, seq)] = kmers
        kmer_counter.update(kmers)
    return seq_kmers, kmer_counter


def dominant_label(ref_kmer, kmers, kmer_counter):
    label_counter = Counter()
    for km in kmers:
        labels = ref_kmer.get(km, ())
        # use square root to avoid weight be 0 because conserved kmers frequency can be larger than 10000
        weight = 1 if kmer_counter[km] <= 10 else 1 / kmer_counter[km] ** 0.5
        for l in labels:
            label_counter[l] += weight
    if not label_counter:
        return None, 0
    return label_counter.most_common(1)[0]


def chimera_check(kmers, ref_kmer, kmer_counter, min_match_thresh):
    mid = len(kmers) // 2
    # avoid too short
    if mid < 1 or len(kmers) - mid < 1:
        return False, None, None
    l_kmers, r_kmers = kmers[:mid], kmers[mid:]
    l_label, l_cnt = dominant_label(ref_kmer, l_kmers, kmer_counter)
    r_label, r_cnt = dominant_label(ref_kmer, r_kmers, kmer_counter)
    l_ratio = l_cnt / len(l_kmers)
    r_ratio = r_cnt / len(r_kmers)
    is_chimera = (l_label != r_label and l_ratio >= min_match_thresh and r_ratio >= min_match_thresh)
    return is_chimera, l_label, r_label


def process_sequence(args):
    (header, seq), kmers, kmer_counter, ref_kmer, min_match_thresh = args
    if not kmers:
        return (header, seq), False
    is_chimera, l_label, r_label = chimera_check(kmers, ref_kmer, kmer_counter, min_match_thresh)
    if is_chimera:
        logger.info(f"Chimera detected: {header.split(' ')[0]}, left label: {l_label}, right label: {r_label}")
    return (header, seq), is_chimera


def find_input_files(group_dir):
    input_files = {}
    for sample_dir in group_dir.iterdir():
        if sample_dir.is_dir():
            merged = sample_dir / "merged.fasta"
            filtered = sample_dir / "filtered.fasta"
            if merged.exists():
                input_files[sample_dir.name] = merged
            elif filtered.exists():
                input_files[sample_dir.name] = filtered
            else:
                logger.warning(f"No valid fasta file for {sample_dir.name}.")
    return input_files


### ====== main function ====== ###
def main(project_path, groups, kmer_db, min_match_thresh, num_proc):

    project_path = Path(project_path)

    summary_path = project_path / "kmer_results" / "kmer_summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"{summary_path} not found.")
    with open(summary_path) as f:
        best_k = json.load(f)["best_k"]

    if not groups:
        raise ValueError("No groups specified.")
    for group in groups:
        group_dir = project_path / group
        if not group_dir.exists() or not group_dir.is_dir():
            logger.warning(f"Group directory {group_dir} not found, skipping.")
            continue

        logger.info(f"{'=' * 50}")
        logger.info(f"Processing group: {group}")

        sample_files = find_input_files(group_dir)

        for sample_name, input_file in sample_files.items():
            logger.info(f"{'-' * 50}")
            logger.info(f"Processing sample {sample_name}")

            non_chimeras_file = group_dir / sample_name / "non_chimeras.fasta"

            #if non_chimeras_file.exists():
            #    logger.info(f"Output files already exist, skipping {sample_name}.")
            #    continue

            seqs = list(read_fasta(input_file))
            seq_kmers, kmer_counter = build_kmer_cache(seqs, best_k)
            logger.info(f"Loaded {len(seqs)} sequences.")

            with open_kmer_db(kmer_db) as conn:
                ref_kmer = fetch_all_kmer_labels(conn, best_k)

            # parallel running
            seq_items = sorted(seq_kmers.items(), key=lambda x: len(x[1]), reverse=True)
            num_processes = min(num_proc or cpu_count(), len(seq_items))

            args_list = [
                ((header, seq), kmers, kmer_counter, ref_kmer, min_match_thresh)
                for (header, seq), kmers in seq_items]

            chimera_count = []
            non_chimeras = []
            with Pool(processes=num_processes) as pool:
                results = pool.map(process_sequence, args_list)
                for (header, seq), is_chimera in results:
                    if is_chimera:
                        chimera_count.append((header, seq))
                    else:
                        non_chimeras.append((header, seq))

            logger.info(
                f"Non-chimeras rate: {len(non_chimeras) / len(seqs):.2%}, {len(non_chimeras)} sequences passed.")
            save_files(non_chimeras_file, non_chimeras)
            logger.info(f"Saved non-chimeric sequences successfully.")


if __name__ == '__main__':
    start = time.time()

    main(**step3_config)
    logger.info(f"Analysis completed in {time.time() - start:.2f}s.")
