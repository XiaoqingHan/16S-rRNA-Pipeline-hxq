from pathlib import Path

#############################
# paths, common params
#############################
root_dir = Path("/Users/hanxiaoqing/Downloads/fq/projects/")
project = "PRJNA555320"  # PRJNA555320
groups = ["CF", "SECC"]
kmer_db = "/Users/hanxiaoqing/Downloads/fq/ref/ref_v16.kmer"
k_list = [15, 21, 25, 31, 35]
seq_mode = "paired"  # paired, single

directory = root_dir / project



#############################
# kmer db
#############################
db_config = dict(
    k_list=k_list,
    kmer_db=kmer_db,
    ref_fasta="/Users/hanxiaoqing/Downloads/fq/ref/HOMD_16S_rRNA_RefSeq_V16.02_full.fasta"
)



#############################
# qc, merge
#############################
step1_config = dict(
    input_dir=directory,
    seq_mode=seq_mode,
    primer_1="GTGCCAGCMGCCGCGGTAA",  # v3-v4 CCTACGGRRBGCASCAGKVRVGAAT  # v4 GTGCCAGCMGCCGCGGTAA
    primer_2="GGACTACHVGGGTWTCTAAT",  # v3-v4 GGACTACNVGGGTWTCTAATCC  # v4 GGACTACHVGGGTWTCTAAT
    min_len=100,
    min_q=30,
    min_overlap=20,
    max_trim_rate=0.15,  # mismatch for primer
    max_offset=10,
    max_n=0.05,
    max_merge_rate=0.05,  # mismatch for overlap
    min_len_m=200,
    max_len_m=500,
    max_n_m=0.01,
    n_procs=6
)



#############################
# select k
#############################
step2_config = dict(
    project_path=directory,
    groups=groups,
    kmer_db=kmer_db,
    k_list=k_list,
    n_sample=10,
    random_seed=42,
    seq_mode=seq_mode,
    batch_size=1000,
    top_n=10
)



#############################
# non-chimera
#############################
step3_config = dict(
    project_path=directory,
    groups=groups,
    kmer_db=kmer_db,
    min_match_thresh=0.8,
    num_proc=None,
)



#############################
# abundance
#############################
step4_config = dict(
    kmer_db=kmer_db,
    project_path=directory,
    groups=groups,
    file="non_chimeras.fasta"
)



#############################
# top n
#############################
step5_config = dict(
    project_path=directory,
    groups=groups,
    top_n=10,
    genus_focus=["Streptococcus", "Veillonella"],
    min_fraction=0.05
)
