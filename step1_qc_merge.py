import time
from config import *
from multiprocessing import Pool
from common_func import *


###############
# this script is run for all samples in a group
# - inputs: fastq.gz
# remove primer - trim()
# quality control - qc()
# merge two sequences - merge_paired_reads()
# qc merged seqs - merge_qc()
# - outputs: merge.fasta, unmerged.fasta
###############


### ====== main functions ====== ###
# calculate average phred score per read, phred=ascii-33 (Illumina >= 1.8)
def phred_score(base_quality):
    if not base_quality:
        return 0
    return sum(ord(i) - 33 for i in base_quality) / len(base_quality)


# calculate the mismatch number of base
def hamming_distance(seq1, seq2):
    count = 0
    if len(seq1) == len(seq2):
        for base1, base2 in zip(seq1, seq2):
            if base1 != base2:
                count += 1
    return count


# remove primer, no need to consider 3'(seq[-len(primer):]) cuz primer usually happens in 5', allow 2 mismatch
# primer not always starts from first position, need slide
def detect_primer(seq, qual, primer, max_trim_rate, max_offset):
    seq = seq.upper()
    primer = primer.upper()
    trimmed_flag = False
    search_range = min(max_offset + 1, len(seq) - len(primer) + 1)
    for i in range(search_range):
        window = seq[i:i + len(primer)]
        if hamming_distance(window, primer) <= max_trim_rate:
            seq, qual = seq[i + len(primer):], qual[i + len(primer):]
            trimmed_flag = True
            break
    return seq, qual, trimmed_flag


# trim primer
def trim(reads, primer, max_trim_rate, max_offset):

    trimmed_reads = []
    trimmed_cnt = 0
    primer_len = len(primer)
    max_hamming_trim = max(1, int(primer_len * max_trim_rate))

    for header, seq, qual in reads:
        seq, qual, trimmed_flag = detect_primer(seq, qual, primer, max_hamming_trim, max_offset)
        # keep reads that detected primer
        if trimmed_flag:
            trimmed_cnt += 1
            trimmed_reads.append((header, seq, qual))

    return trimmed_cnt, trimmed_reads


# quality control, use Q30 strategy
def qc(reads, min_len, min_q, max_n):
    filtered_reads = []
    passed = 0
    filter_reasons = {"low_quality": 0, "short_length": 0, "high_N": 0}

    for header, seq, qual in reads:
        n_content = seq.count('N') / len(seq)
        avg_q = phred_score(qual)
        if len(seq) >= min_len and avg_q >= min_q and n_content <= max_n:
            filtered_reads.append((header, seq, qual))
            passed += 1
        else:
            if avg_q < min_q:
                filter_reasons["low_quality"] += 1
            elif len(seq) < min_len:
                filter_reasons["short_length"] += 1
            elif n_content > max_n:
                filter_reasons["high_N"] += 1

    logger.info(f"QC breakdown: {filter_reasons}")
    return passed, filtered_reads


# get reversed, complementary sequence of reads 2, including IUPAC
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
                  'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B'}
    return ''.join(complement.get(base.upper()) for base in reversed(seq))


# normalize headers
def normal_header(h):
    h = h.lstrip('@>').split()[0]
    return h


# merge reads1 and reads2, allow 1 mismatch, R1+reverse_complement(R2)
def merge_paired_reads(reads1, reads2, min_overlap, max_merge_rate):
    reads1 = list(reads1)
    reads2 = list(reads2)
    # record the original header
    header_map1 = {normal_header(h): h for h, _, _ in reads1}
    header_map2 = {normal_header(h): h for h, _, _ in reads2}

    # unique keys and values
    seq1 = {}
    for h, r, _ in reads1:
        norm_h = normal_header(h)
        if norm_h not in seq1:
            seq1[norm_h] = r.upper()
    seq2 = {}
    for h, r, _ in reads2:
        norm_h = normal_header(h)
        if norm_h not in seq2:
            seq2[norm_h] = r.upper()

    common_headers = []
    for norm_h in seq1:
        if norm_h in seq2:
            # verify headers
            orig_h1 = header_map1[norm_h]
            orig_h2 = header_map2[norm_h]
            if orig_h1.split()[0] == orig_h2.split()[0]:
                common_headers.append(norm_h)
    for norm_h in common_headers[:3]:
        logger.info("Normalized common headers: %r ↔ %r", header_map1[norm_h], header_map2[norm_h])

    merged_seqs = {}
    unmerged1 = {}
    unmerged2 = {}
    for h in common_headers:
        s1, s2 = seq1[h], seq2[h]
        merged = False
        for i in range(min(len(s1), len(s2)), min_overlap - 1, -1):
            if hamming_distance(s1[-i:], s2[:i]) / i <= max_merge_rate:
                merged_seqs[f">{h}"] = s1 + s2[i:]
                merged = True
                break
        if not merged:
            unmerged1[f">{h}/1"] = s1
            unmerged2[f">{h}/2"] = reverse_complement(s2)  # keep original

    merged_list = [(header, seq) for header, seq in merged_seqs.items()]
    unmerged1_list = [(header, seq) for header, seq in unmerged1.items()]
    unmerged2_list = [(header, seq) for header, seq in unmerged2.items()]

    return merged_list, unmerged1_list, unmerged2_list


# check the quality of merged seqs
def merge_qc(merged_seqs, min_len_m, max_len_m, max_n_m):
    filtered = []
    for header, seq in merged_seqs:
        seq_len = len(seq)
        n_content = seq.count('N') / seq_len
        if (min_len_m <= seq_len <= max_len_m) and (n_content <= max_n_m):
            filtered.append((header, seq))
    return filtered


# run for single sample
def process_sample(sample_dir, params, seq_mode):
    sample_name = sample_dir.name
    logger.info(f"{'=' * 50}")
    logger.info(f"Analyzing sample: {sample_name}...")

    # if (sample_dir / "merged.fasta").exists() or (sample_dir / "filtered.fasta").exists():
    #    logger.info(f"Output files already exist, skipping {sample_name}.")
    #    return

    if seq_mode == "paired":
        try:
            fastq_r1 = next(sample_dir.glob('*_1.fastq.gz'))
            fastq_r2 = next(sample_dir.glob('*_2.fastq.gz'))
        except StopIteration:
            logger.error(f"Cannot find paired-end fastq files.")
            return
        out_qc1 = sample_dir / "filtered_R1.fasta"
        out_qc2 = sample_dir / "filtered_R2.fasta"
        out_merged = sample_dir / "merged.fasta"
        out_unmerge1 = sample_dir / "unmerged_1.fasta"
        out_unmerge2 = sample_dir / "unmerged_2.fasta"

        reads1 = list(read_fastq(fastq_r1))
        reads2 = list(read_fastq(fastq_r2))
        count1, count2 = len(reads1), len(reads2)
        logger.info(f"Processed R1: {count1}, R2: {count2} read pairs.")

        # trim primer
        trim_cnt1, trim1 = trim(reads1,
                                primer=params["primer_1"],
                                max_trim_rate=params["max_trim_rate"],
                                max_offset=params["max_offset"])

        trim_cnt2, trim2 = trim(reads2,
                                primer=params["primer_2"],
                                max_trim_rate=params["max_trim_rate"],
                                max_offset=params["max_offset"])
        logger.info(f"R1 trimmed rate:  {trim_cnt1 / count1:.2%}, {trim_cnt1} reads trimmed primer.")
        logger.info(f"R2 trimmed rate:  {trim_cnt2 / count2:.2%}, {trim_cnt2} reads trimmed primer.")

        # qc
        passed1, filt1 = qc(
            trim1,
            min_len=params["min_len"],
            min_q=params["min_q"],
            max_n=params["max_n"]
        )
        logger.info(f"QC rate of R1: {passed1 / count1:.2%}, {passed1} reads passed.")

        passed2, filt2 = qc(
            trim2,
            min_len=params["min_len"],
            min_q=params["min_q"],
            max_n=params["max_n"]
        )
        logger.info(f"QC rate of R2: {passed2 / count2:.2%}, {passed2} reads passed.")

        # merge
        filt2_rev = [(h, reverse_complement(s), q) for h, s, q in filt2]
        merged, un1, un2 = merge_paired_reads(
            filt1, filt2_rev,
            min_overlap=params["min_overlap"],
            max_merge_rate=params["max_merge_rate"]
        )
        logger.info(f"Merge success rate: {len(merged) / min(passed1, passed2):.2%}, {len(merged)} read pairs merged.")

        # qc again
        qc_merged = merge_qc(
            merged,
            min_len_m=params["min_len_m"],
            max_len_m=params["max_len_m"],
            max_n_m=params["max_n_m"]
        )
        logger.info(f"Post-merge QC: {len(qc_merged)}/{len(merged)} merged sequences kept.")

        lengths = [len(seq) for _, seq in qc_merged]
        logger.info(
            f"Merged length stats: min={min(lengths)}, max={max(lengths)}, avg={sum(lengths) / len(lengths):.1f}bp.")

        save_files(out_qc1, [(r[0].replace('@', '>'), r[1]) for r in filt1])
        save_files(out_qc2, [(r[0].replace('@', '>'), r[1]) for r in filt2])
        save_files(out_merged, qc_merged)
        save_files(out_unmerge1, un1)
        save_files(out_unmerge2, un2)

    elif seq_mode == "single":
        try:
            fastq = next(sample_dir.glob('*.fastq.gz'))
        except StopIteration:
            logger.error(f"{sample_name}: Cannot find single-end fastq, skip.")
            return

        out_qc = sample_dir / "filtered.fasta"
        reads = list(read_fastq(fastq))
        count = len(reads)
        logger.info(f"Processed single-end: {count}.")

        trim_cnt, trim_ = trim(reads,
                               primer=params["primer_1"],
                               max_trim_rate=params["max_trim_rate"],
                               max_offset=params["max_offset"])

        passed, filt = qc(
            trim_,
            min_len=params["min_len"],
            min_q=params["min_q"],
            max_n=params["max_n"]
        )
        logger.info(f"QC rate of R1: {passed / count:.2%}, {passed} reads passed.")

        fasta_reads = [(h.replace('@', '>'), seq) for h, seq, _ in filt]
        save_files(out_qc, fasta_reads)

    logger.info(f"{sample_name} done.")


# paralleling
def process_wrapper(sample_dir, seq_mode):
    try:
        process_sample(Path(sample_dir), step1_config, seq_mode)
    except Exception as e:
        logger.error(f"{sample_dir.name} failed: {e}", exc_info=True)


### ====== main ====== ###
def main():
    start = time.time()

    seq_mode = step1_config["seq_mode"]
    n_procs = step1_config["n_procs"]

    sample_dirs = []
    for d in directory.iterdir():
        if d.is_dir() and not d.name.startswith('.'):
            sample_dirs.extend([s for s in d.iterdir()])
    logger.info(f"Found {len(sample_dirs)} samples across all groups")

    task_args = [(d, seq_mode) for d in sample_dirs]

    with Pool(processes=n_procs) as pool:
        pool.starmap(process_wrapper, task_args)

    logger.info(f"All samples completed in {time.time() - start:.2f}s.")


if __name__ == '__main__':
    main()
