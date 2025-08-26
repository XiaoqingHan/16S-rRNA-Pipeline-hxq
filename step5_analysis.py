import re
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from config import step5_config
from common_func import logger
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu
import time
from statsmodels.stats.multitest import multipletests


# transform abundance table, combine all samples into one matrix
def abundance_matrix(dir_path):
    base = Path(dir_path)

    # relative and absolute abundance for species and genera
    abs_dict_spe, rel_dict_spe = {}, {}
    abs_dict_gen, rel_dict_gen = {}, {}

    for sample_dir in base.iterdir():
        if not sample_dir.is_dir() or sample_dir.name == "kmer_results":
            continue
        sid = sample_dir.name
        abund = sample_dir / "taxon_abundance.csv"
        if not abund.exists():
            continue

        df = pd.read_csv(abund)

        # species
        species_abs = df.groupby('Taxon')['Abundance'].sum()
        total_spe = species_abs.sum()
        for taxon, cnt in species_abs.items():
            abs_dict_spe.setdefault(taxon, {})[sid] = cnt
            rel_dict_spe.setdefault(taxon, {})[sid] = cnt / total_spe * 100

        # genus
        df['Genus'] = df['Taxon'].str.split().str[0]
        genus_abs = df.groupby('Genus')['Abundance'].sum()
        total_gen = genus_abs.sum()
        for taxon, cnt in genus_abs.items():
            abs_dict_gen.setdefault(taxon, {})[sid] = cnt
            rel_dict_gen.setdefault(taxon, {})[sid] = cnt / total_gen * 100

    abs_df_spe = pd.DataFrame.from_dict(abs_dict_spe, orient='index').fillna(0)
    rel_df_spe = pd.DataFrame.from_dict(rel_dict_spe, orient='index').fillna(0)
    abs_df_gen = pd.DataFrame.from_dict(abs_dict_gen, orient='index').fillna(0)
    rel_df_gen = pd.DataFrame.from_dict(rel_dict_gen, orient='index').fillna(0)

    # sort by sample id
    def _sort_key(col_name):
        m = re.search(r'([A-Za-z]+)(\d+)$', col_name)
        return int(m.group(2)) if m else col_name

    sorted_cols = sorted(abs_df_spe.columns, key=_sort_key)

    abs_df_spe = abs_df_spe[sorted_cols]
    rel_df_spe = rel_df_spe[sorted_cols]
    abs_df_gen = abs_df_gen[sorted_cols]
    rel_df_gen = rel_df_gen[sorted_cols]

    abs_df_spe.index.name = 'Species'
    rel_df_spe.index.name = 'Species'
    abs_df_gen.index.name = 'Genus'
    rel_df_gen.index.name = 'Genus'

    return abs_df_spe, rel_df_spe, abs_df_gen, rel_df_gen


# only taxa with non-zero counts in at least 5% of the samples are retained
def filter_min_samples(df, min_fraction):
    n_samples = df.shape[1]
    min_count = max(1, int(n_samples * min_fraction))
    keep = df.gt(0).sum(axis=1) >= min_count
    filtered_df = df.loc[keep]
    n_filtered = (~keep).sum()
    return df.loc[keep], n_filtered


# calculate top n genera, remark the rest as other
def top_genera(df, n):
    groups = {col: re.sub(r"\d+", "", col) for col in df.columns}
    # calculate mean abundance
    group_mean = df.T.groupby(groups).mean().T
    # calculate top n taxa for each group
    top_taxa = set()
    for g in group_mean.columns:
        top_n = group_mean[g].nlargest(n).index.tolist()
        top_taxa.update(top_n)

    # group low abundance samples as other and insert other into the last row
    all_groups = group_mean.copy()
    all_groups.loc['Other'] = all_groups.loc[~all_groups.index.isin(top_taxa)].sum()
    all_groups = all_groups.loc[list(top_taxa) + ['Other']]

    # calculate top n taxa for each sample
    top_taxa_ = set()
    for s in df.columns:
        top_n = df[s].nlargest(n).index
        top_taxa_.update(top_n)
    df_ = df.copy()
    all_samples = df_.loc[df_.index.isin(top_taxa_)].copy()
    all_samples.loc['Other'] = df_.loc[~df_.index.isin(top_taxa_)].sum()

    return all_groups, all_samples


# plot the composition of each sample
def barplot(df, outpath, fname="composition.png"):
    df_plot = df.T

    plt.figure(figsize=(10, 8))
    ax = df_plot.plot(
        kind='bar',
        stacked=True,
        colormap='tab20',
        edgecolor='white',
        width=0.8,
        alpha=0.9,
        ax=plt.gca()
    )
    ax.set_title("Top Genera Composition")
    ax.set_xlabel('Samples', fontsize=8)
    ax.set_ylabel('Relative Abundance (%)', fontsize=8)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha='right',
        fontsize=6,
        rotation_mode='anchor'
    )
    ax.tick_params(axis='y', labelsize=6)
    ax.legend(
        bbox_to_anchor=(1, 1),
        title='Genus',
        fontsize=5
    )
    plt.tight_layout()
    plt.savefig(Path(outpath) / fname, dpi=300, bbox_inches='tight')
    plt.close()


def boxplot(df, output, genus_list, test_results_df):
    df_filtered = df[df['Genus'].isin(genus_list)].copy()

    groups = sorted(df_filtered['Group'].unique())
    palette = sns.color_palette("husl", len(groups))
    palette = dict(zip(groups, palette))

    plt.figure(figsize=(12, 8))
    g = sns.catplot(
        data=df_filtered,
        x='Group',
        y='Abundance',
        hue='Group',
        col='Genus',
        col_order=genus_list,
        kind='box',
        palette=palette,
        height=6,
        aspect=1,
        sharey=False
    )

    for ax in g.axes.flat:
        genus_name = ax.get_title().split(' = ')[1]
        sig_info = test_results_df[test_results_df['Genus'] == genus_name]

        if not sig_info.empty:
            pval = sig_info['pvalue_adj'].values[0]
            sig = sig_info['significant'].values[0]

            if sig:
                ax.set_title(f"{genus_name} (p_adj={pval:.3e}, *)")
            else:
                ax.set_title(f"{genus_name} (p_adj={pval:.3e})")

    g.fig.subplots_adjust(top=0.85)
    g.fig.suptitle("Non-parametric Test (FDR corrected)", fontsize=16)
    g.set_axis_labels("", "Relative Abundance (%)")
    g.add_legend(title="Group")

    plt.savefig(Path(output) / "boxplot.png", dpi=300, bbox_inches='tight')
    plt.close()



def heatmap(df, output, top_n):
    mean_abundance = df.mean(axis=1).sort_values(ascending=False)

    top_genera = mean_abundance.head(top_n).index.tolist()
    other_genera = mean_abundance.index.difference(top_genera)

    plot_data = df.loc[top_genera].copy()
    plot_data.loc['Other'] = df.loc[other_genera].sum()

    groups = {col: re.sub(r"\d+", "", col) for col in plot_data.columns}
    group_mean = plot_data.T.groupby(groups).mean().T

    non_other = group_mean.drop(index=['Other'], errors='ignore')
    sorted_non_other = non_other.sum(axis=1).sort_values(ascending=False).index.tolist()
    final_order = sorted_non_other + (['Other'] if 'Other' in group_mean.index else [])

    group_mean = group_mean.loc[final_order]

    plt.figure(figsize=(8, 12))
    plt.title("Oral Microbiome Core Genera Abundance", fontsize=16)
    ax = sns.heatmap(
        group_mean,
        annot=True,
        fmt='.2f',
        cmap='YlOrRd',
        cbar_kws={'label': 'Relative abundance (%)'},
        linewidths=0.5,
        annot_kws={'size': 9}
    )

    ax.set_yticks([i + 0.5 for i in range(len(group_mean.index))])
    ax.set_yticklabels([f"$\\it{{{g}}}$" for g in group_mean.index], rotation=0)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(Path(output) / "heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()


# non-parametric tests
def nonpara_test(df, genera):
    df_filtered = df[df['Genus'].isin(genera)].copy()

    df_melted = df_filtered.melt(
        id_vars='Genus',
        var_name='Sample',
        value_name='Abundance'
    )

    df_melted['Group'] = df_melted['Sample'].str.extract(r'([A-Z]+)')[0]
    groups = sorted(df_melted['Group'].unique())

    test_results = {}
    pvals = []
    for genus in genera:
        genus_data = df_melted[df_melted['Genus'] == genus]
        group_data = [genus_data[genus_data['Group'] == g]['Abundance'] for g in groups]
        # Mann-Whitney U test for two groups
        if len(groups) == 2:
            stat, pval = mannwhitneyu(*group_data)
            test_name = "Mann-Whitney U"
        # Kruskal-Wallis test for more than 2 groups
        elif len(groups) > 2:
            stat, pval = kruskal(*group_data)
            test_name = "Kruskal-Wallis"
        else:
            raise ValueError("At least two groups.")

        test_results[genus] = {
            'test': test_name,
            'statistic': stat,
            'pvalue': pval
        }
        pvals.append(pval)

    reject, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')
    for i, genus in enumerate(genera):
        test_results[genus]['pvalue_adj'] = pvals_corrected[i]
        test_results[genus]['significant'] = reject[i]

    return test_results


def main(project_path, groups, top_n, genus_focus, min_fraction):
    project_path = Path(project_path)
    out_path = project_path / "abund_results"
    out_path.mkdir(exist_ok=True)

    group_top_rel = {}

    for group_name in groups:
        g_dir = project_path / group_name
        if not g_dir.is_dir():
            logger.warning(f"Group folder not found: {g_dir}")
            continue
        logger.info(f"Processing group: {group_name}")

        sample_dirs = [d for d in g_dir.iterdir() if d.is_dir()]
        if not sample_dirs:
            logger.warning(f"No samples found in group {group_name}")
            continue

        abs_gen_dict, rel_gen_dict = {}, {}

        for sample_dir in sample_dirs:
            abund_file = sample_dir / "taxon_abundance.csv"
            if not abund_file.exists():
                logger.warning(f"No abundance file for sample {sample_dir.name}")
                continue
            df = pd.read_csv(abund_file)
            df['Genus'] = df['Taxon'].str.split().str[0]
            genus_abs = df.groupby('Genus')['Abundance'].sum()
            total_gen = genus_abs.sum()
            for taxon, cnt in genus_abs.items():
                abs_gen_dict.setdefault(taxon, {})[sample_dir.name] = cnt
                rel_gen_dict.setdefault(taxon, {})[sample_dir.name] = cnt / total_gen * 100

        rel_gen = pd.DataFrame.from_dict(rel_gen_dict, orient='index').fillna(0)
        abs_gen = pd.DataFrame.from_dict(abs_gen_dict, orient='index').fillna(0)

        # filter genera
        rel_gen, n_filtered_gen = filter_min_samples(rel_gen, min_fraction)
        logger.info(f"{group_name}: filtered out {n_filtered_gen} genera")

        abs_gen.to_csv(out_path / f"{group_name}_abs_abund_gen.csv")
        rel_gen.to_csv(out_path / f"{group_name}_rel_abund_gen.csv")

        # top n, barplot
        top_gen, compo_gen = top_genera(rel_gen, top_n)
        top_gen.to_csv(out_path / f"{group_name}_top_n_gen.csv")
        barplot(compo_gen, out_path, fname=f"{group_name}_composition.png")

        group_top_rel[group_name] = rel_gen

    if not group_top_rel:
        logger.error("No group data found. Exiting.")
        return

    shared_genera = list(set.intersection(*(set(df.index) for df in group_top_rel.values())))
    all_rel = pd.concat([df.loc[shared_genera] for df in group_top_rel.values()], axis=1)

    mean_abundance = all_rel.mean(axis=1).sort_values(ascending=False)
    top_genera_for_stats = mean_abundance.head(top_n).index.tolist()
    logger.info(f"Top {top_n} genera: {top_genera_for_stats}")
    # heatmap
    heatmap(all_rel, out_path, top_n)

    all_samples = all_rel.reset_index().rename(columns={'index': 'Genus'})

    # top n + focused genera
    all_genera_to_test = list(set(top_genera_for_stats + genus_focus))
    logger.info(f"Genera for statistical testing: {all_genera_to_test}")

    # non-parametric test for
    test_results = nonpara_test(all_samples, all_genera_to_test)

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    pvals = [test_results[genus]['pvalue'] for genus in all_genera_to_test]
    reject, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')

    for i, genus in enumerate(all_genera_to_test):
        test_results[genus]['pvalue_adj'] = pvals_corrected[i]
        test_results[genus]['significant'] = reject[i]
        test_results[genus]['type'] = 'top' if genus in top_genera_for_stats else 'focus'

    df_results = pd.DataFrame.from_dict(test_results, orient='index').reset_index().rename(columns={'index': 'Genus'})
    df_results.to_csv(out_path / "non_parametric.csv", index=False)

    if genus_focus:
        logger.info(f"Creating boxplots for focus genera: {genus_focus}")

        all_samples_long = all_samples.melt(
            id_vars='Genus',
            var_name='Sample',
            value_name='Abundance'
        )
        all_samples_long['Group'] = all_samples_long['Sample'].str.extract(r'([A-Z]+)')[0]
        boxplot(all_samples_long, out_path, genus_focus, df_results)

    logger.info("All analysis completed.")


if __name__ == "__main__":
    start = time.time()

    main(**step5_config)
    logger.info(f"All steps completed in {time.time() - start:.2f}s.")

'''
not necessary analysis
# calculate alpha diversity, the diversity within one sample, can test data quality and sample balance
def alpha_diversity(matrix):

    results = {'Sample': [], 'shannon': [], 'simpson': []}
    for sample in matrix.columns:
        counts = np.array(matrix[sample])
        nonzero = counts[counts>0]
        # Shannon index: H=-sum(p_i*log p_i)
        p = nonzero/nonzero.sum()
        H = -np.sum(p*np.log(p))
        # Simpson index: D=1-sum(p_i^2)
        D = 1-np.sum(p**2)
        results['Sample'].append(sample)
        results['Shannon'].append(H)
        results['Simpson'].append(D)
    alpha_df = pd.DataFrame(results).set_index('Sample')

    return alpha_df

# beta diversity, measuring differences in community composition between samples
# correlation analysis, the interactions between genus

'''
