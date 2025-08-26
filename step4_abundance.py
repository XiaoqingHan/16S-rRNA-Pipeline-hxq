import time
import pandas as pd
from collections import Counter
from config import step4_config
from common_func import *
from db_utils import *
from multiprocessing import Pool, cpu_count
import json


# get taxon of each seq by kmer hits
def classify(seq, k, kmer_map):
    votes = Counter()
    for kmer in yield_kmers(seq, k):
        if kmer in kmer_map:
            votes.update(kmer_map[kmer])
    if not votes:
        return "Unclassified"
    return votes.most_common(1)[0][0]


# calculate abundance
def compute_abundance(seqs, sample_id, k, kmer_map):
    counter = Counter()
    for seq in seqs:
        taxon = classify(seq, k, kmer_map)
        counter[taxon] += 1
    total = sum(counter.values())
    rows = []
    for taxon, cnt in counter.items():
        rows.append({
            'Sample_ID': sample_id,
            'Taxon': taxon,
            'Abundance': cnt,
            'Relative_Abundance': cnt / total * 100
        })
    df = pd.DataFrame(rows)
    df = df.sort_values('Abundance', ascending=False).reset_index(drop=True)
    return df


def process_sample(kmer_db, best_k, sample_dir, file):
    sample_id = sample_dir.name
    non_chi = sample_dir / file

    if not non_chi.exists():
        logger.info(f"Skipping {sample_id}, file not found.")
        return sample_id, None

    try:
        logger.info(f"{'=' * 50}")
        logger.info(f"Processing {sample_id}...")
        seqs = [seq for _, seq in read_fasta(non_chi)]

        with open_kmer_db(kmer_db) as conn:
            kmer_map = fetch_all_kmer_labels(conn, best_k)

        abundance_df = compute_abundance(seqs, sample_id, best_k, kmer_map)

        output_file = sample_dir / "taxon_abundance.csv"
        abundance_df.to_csv(output_file, index=False)
        logger.info(f"Saved {sample_id} taxon abundance file successfully.")

        return sample_id, len(seqs)

    except Exception as e:
        logger.error(f"Error in {sample_id}: {str(e)}")
        return sample_id, str(e)


### ====== main ====== ###
def main(kmer_db, project_path, groups, file):
    project_path = Path(project_path)
    group_dirs = [project_path / g for g in groups]

    for g in group_dirs:
        if not g.is_dir():
            raise FileNotFoundError(f"Group folder not found: {g}")

    # read best_k from previous result
    summary_path = project_path / "kmer_results" / "kmer_summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"kmer_summary.json not found: {summary_path}")

    with open(summary_path) as f:
        best_k = json.load(f).get("best_k")
        if best_k is None:
            raise ValueError("best_k not found in kmer_summary.json")

    sample_dirs = [d for g in group_dirs for d in g.iterdir() if d.is_dir()]
    if not sample_dirs:
        logger.warning("No sample directories found in any group.")
        return

    params = [(kmer_db, best_k, d, file) for d in sample_dirs]

    with Pool(processes=min(cpu_count(), len(sample_dirs))) as pool:
        pool.starmap(process_sample, params)

    logger.info(f"Processed {len(sample_dirs)} samples using k={best_k}.")


# edit here
if __name__ == '__main__':
    start = time.time()

    main(**step4_config)

    logger.info(f"Steps completed in {time.time() - start:.2f}s.")
