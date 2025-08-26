import random
import time
import json
import pandas as pd
import numpy as np
from itertools import combinations
from collections import Counter
from config import step2_config
from db_utils import *
from common_func import *
from multiprocessing import Pool, cpu_count


###############
# this script is run for multiple samples in a group
# - inputs: ref.fasta, merged.fasta, unmerged.fasta
# construct kmer database - build_kmer_db()
# select samples randomly - select_samples()
# test kmer influence - kmer_impact()
# - outputs: top_n.csv, kmer_jaccard.csv
###############


### ====== base functions ====== ###
# select samples randomly
def select_samples(directory, n, random_seed):
    dir = Path(directory)
    group_samples = {}
    samples = [s for s in dir.iterdir() if s.is_dir()]
    if samples:
        group_samples[dir.name] = samples
    if random_seed is not None:
        random.seed(random_seed)
    selected = {}
    for group, samples in group_samples.items():
        selected[group] = random.sample(samples, n)
    return selected


### ====== main functions ====== ###
# connect kmer db and following analysis
def load_kmer_maps(kmer_db, k_list):
    kmer_maps = {}
    with open_kmer_db(kmer_db) as conn:
        for k in k_list:
            table_name = f"kmer_labels_k{k}"
            kmer_maps[k] = fetch_all_kmer_labels(conn, k, table_name)
    return kmer_maps


# calculate pairwise jaccard similarities between two top_n results
def jaccard(results):
    ks = sorted(results.keys())
    jacc_map = {k: {} for k in ks}
    for k1, k2 in combinations(ks, 2):
        set1 = set(results[k1]['top_taxon'])
        set2 = set(results[k2]['top_taxon'])
        inter = len(set1 & set2)
        union = len(set1 | set2)
        j = inter / union if union > 0 else 0.0
        jacc_map[k1][k2] = j
        jacc_map[k2][k1] = j
    for k in ks:
        jacc_map[k][k] = 1.0
    return jacc_map


def process_batch(args):
    batch, k, kmer_map = args
    taxon_counter = Counter()
    for header, seq in batch:
        votes = Counter()
        for km in yield_kmers(seq, k):
            labels = kmer_map.get(km, [])
            if labels:
                weight = 1.0 / len(labels)
                for label in labels:
                    votes[label] += weight
        if votes:
            dominant = votes.most_common(1)[0][0]
        else:
            dominant = "Unclassified"
        taxon_counter[dominant] += 1
    return taxon_counter


# compare the influence of different kmers
def kmer_impact(seqs, kmer_maps, k_list, batch_size, top_n):
    results = {}
    seq_list = list(seqs)
    logger.info(f"Total sequences to process: {len(seq_list)}")

    k_times = {} # record runtime at different k

    for k in k_list:
        logger.info(f"Starting analysis for k={k}")
        start_k = time.time()
        taxon_counter = Counter()
        batches = []
        for i in range(0, len(seq_list), batch_size):
            batch = seq_list[i:i + batch_size]
            batches.append((batch, k, kmer_maps[k]))

        num_processes = min(cpu_count(), len(batches))

        # process batches in parallel
        with Pool(processes=num_processes) as pool:
            batch_results = pool.map(process_batch, batches)

        # merge results from all batches
        for batch_counter in batch_results:
            taxon_counter.update(batch_counter)

        top = taxon_counter.most_common(top_n)
        results[k] = {
            'top_taxon': [g for g, _ in top],
            'top_counts': [c for _, c in top],
            'unique_taxon': len(taxon_counter)
        }

        end_k = time.time()
        k_times[k] = end_k - start_k

        logger.info(f"Found {results[k]['unique_taxon']} taxon, time used: {k_times[k]:.2f}s.")

    return results, k_times


def main(project_path, groups, kmer_db, k_list, n_sample, random_seed, seq_mode, batch_size, top_n):

    if seq_mode == "single":
        file = "filtered.fasta"
    elif seq_mode == "paired":
        file = "merged.fasta"

    res_dir = project_path / "kmer_results"
    res_dir.mkdir(exist_ok=True)

    if not groups or len(groups) == 0:
        raise ValueError("Must specify at least one group name.")

    groups = [project_path / g for g in groups]
    for g in groups:
        if not g.is_dir():
            raise FileNotFoundError(f"Group folder not found: {g}")

    logger.info(f"Selected groups: {[g.name for g in groups]}")

    all_topn = []
    all_jacc = {}

    kmer_maps = load_kmer_maps(kmer_db, k_list)

    for group in groups:
        samples = select_samples(group, n_sample, random_seed).get(group.name, [])
        logger.info(f"{'=' * 50}")
        logger.info(f"Selected samples in group {group.name}: {[s.name for s in samples]}")

        if not samples:
            continue

        for sample in samples:
            logger.info(f"{'-' * 50}")
            sample_name = f"{group.name}/{sample.name}"
            logger.info(f"Processing sample: {sample_name}")
            path = sample / file
            if not path.exists():
                logger.warning(f"Skipping {sample_name}: {file} not found")
                continue
            seqs = read_fasta(path)
            results, k_times = kmer_impact(seqs, kmer_maps, k_list, batch_size, top_n)

            # 存 Jaccard
            all_jacc[sample_name] = jaccard(results)

            # 存 top-n 和时间
            for k in k_list:
                all_topn.append({
                    "group": group.name,
                    "sample": sample.name,
                    "k": k,
                    "top_taxa": results[k]["top_taxon"],
                    "counts": results[k]["top_counts"],
                    "unique_taxa": results[k]["unique_taxon"],
                    "time_sec": k_times[k]
                })

    all_topn = pd.DataFrame(all_topn)
    all_topn.to_csv(res_dir / "top_n.csv", index=False)

    logger.info(f"{'=' * 50}")

    scores = {}  # scores[k] = [list of all jaccard values for that k vs others]
    for sample, s_mat in all_jacc.items():
        for k, k_mat in s_mat.items():
            for k_, j in k_mat.items():
                if k_ == k:
                    continue
                scores.setdefault(k, []).append(j)

    avg_scores = {int(k): np.mean(vals) for k, vals in scores.items()}
    std_scores = {int(k): np.std(vals) for k, vals in scores.items()}
    avg_time = all_topn.groupby('k')['time_sec'].mean().to_dict()
    # optimal k: largest jaccard similarity, most stable, minimal runtime
    best_k = max(avg_scores,key=lambda k: (avg_scores[k],-std_scores[k], -avg_time.get(k, float('inf'))))

    logger.info("K-mer comparison results:")
    for k in sorted(avg_scores, key=int):
        time_sec = avg_time.get(k, 0.0)
        logger.info(f"k={k:>2} | avg Jaccard={avg_scores[k]:.3f} | std Jaccard={std_scores[k]:.3f} | avg time={time_sec:.2f}s")
    logger.info(f"=> Optimal k (global across all groups) = {best_k}")

    summary = {
        "avg_scores": avg_scores,
        "std_scores": std_scores,
        "avg_time": avg_time,
        "best_k": best_k
    }
    with open(res_dir / "kmer_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

if __name__ == '__main__':
    start = time.time()

    main(**step2_config)
    logger.info(f"All samples tested with different kmers completed in {time.time() - start:.2f}s")

