#!/usr/bin/env python
"""
filter_kegg_decoder.py

Takes the full TSV produced by KEGG-decoder (all ~150 pathways) and
generates:
  1. A filtered TSV with only the informative / selected pathways
     (main figure for publication)
  2. An SVG/PNG heatmap of the filtered version
  3. (optional) The complete original stays untouched for the supplementary

Usage:
    python filter_kegg_decoder.py -i KEGG_decoder_output.tsv -o filtered_output.tsv --fig heatmap.png --format png --dpi 300 --mode variance_category

Main options:
    --mode variance          removes zero-variance columns (default)
    --mode threshold         keeps only pathways above a completeness threshold
    --mode select            uses a hand-curated list (see SELECTED_PATHWAYS below)
    --mode category          aggregates pathways into functional macro-categories
    --mode category_select   keeps individual pathways from chosen categories,
                             in the order given, plotted as separated blocks
    --mode variance_category filters by variance and groups remaining pathways
                             into macro-categories (plus 'Other' for leftovers)
    --format {png,svg}       export format (default: png)
    --dpi INT                image resolution in DPI (default: 300)

Dependencies: pandas, seaborn, matplotlib, scipy
"""

import argparse
import os
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial import distance


# ---------------------------------------------------------------------------
# Edit this list for the hand-curated editorial selection (--mode select)
# ---------------------------------------------------------------------------
SELECTED_PATHWAYS = [
    "glycolysis", "gluconeogenesis", "TCA Cycle",
    "dissim nitrate reduction", "DNRA", "nitrite reduction",
    "nitric oxide reduction", "nitrous-oxide reduction", "nitrogen fixation",
    "sulfur assimilation", "dissimilatory sulfate < > APS",
    "dissimilatory sulfite < > sulfide", "thiosulfate oxidation",
    "sulfide oxidation",
    "NiFe hydrogenase", "NAD-reducing hydrogenase",
    "Photosystem II", "Photosystem I",
    "Wood-Ljungdahl", "rTCA Cycle", "CBB Cycle",
    "F-type ATPase", "V-type ATPase",
]

# ---------------------------------------------------------------------------
# Functional macro-categories for aggregation (--mode category / category_select / variance_category)
# Extend/edit based on what you need to show
# ---------------------------------------------------------------------------
CATEGORIES = {
    "Nitrogen cycle": [
        "dissim nitrate reduction", "nitrite oxidation", "DNRA",
        "nitrite reduction", "nitric oxide reduction",
        "nitrous-oxide reduction", "nitrogen fixation",
        "hydroxylamine oxidation", "ammonia oxidation (amo/pmmo)",
    ],
    "Sulfur cycle": [
        "sulfur assimilation", "dissimilatory sulfate < > APS",
        "dissimilatory sulfite < > APS", "dissimilatory sulfite < > sulfide",
        "thiosulfate oxidation", "alt thiosulfate oxidation doxAD",
        "alt thiosulfate oxidation tsdA", "thiosulfate disproportionation",
        "sulfur reductase sreABC", "thiosulfate/polysulfide reductase",
        "sulfhydrogenase", "sulfur disproportionation", "sulfur dioxygenase",
        "sulfite dehydrogenase", "sulfide oxidation",
        "sulfite dehydrogenase (quinone)",
    ],
    "Hydrogen metabolism": [
        "NiFe hydrogenase", "membrane-bound hydrogenase",
        "ferredoxin hydrogenase", "hydrogen:quinone oxidoreductase",
        "NAD-reducing hydrogenase", "NADP-reducing hydrogenase",
        "NiFe hydrogenase Hyd-1",
    ],
    "Central carbon metabolism": [
        "glycolysis", "gluconeogenesis", "TCA Cycle",
        "Entner-Doudoroff Pathway",
    ],
    "Carbon fixation": [
        "CBB Cycle", "rTCA Cycle", "Wood-Ljungdahl",
        "3-Hydroxypropionate Bicycle", "4-Hydroxybutyrate/3-hydroxypropionate",
    ],
    "Photosynthesis": [
        "Photosystem II", "Photosystem I", "Cytochrome b6/f complex",
        "anoxygenic type-II reaction center", "anoxygenic type-I reaction center",
    ],
    "Oxidative phosphorylation": [
        "F-type ATPase", "V-type ATPase", "NADH-quinone oxidoreductase",
        "Cytochrome c oxidase, cbb3-type", "Cytochrome bd complex",
        "Cytochrome o ubiquinol oxidase", "Cytochrome c oxidase",
        "Ubiquinol-cytochrome c reductase",
    ],
    "Methanogenesis / methane ox": [
        "Methanogenesis via methanol", "Methanogenesis via dimethylamine",
        "Methanogenesis via methylamine", "Methanogenesis via trimethylamine",
        "Methanogenesis via acetate", "Methanogenesis via CO2",
        "Coenzyme M reduction to methane", "Soluble methane monooxygenase",
        "methanol dehydrogenase",
    ],
}


def load_table(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df


def filter_variance(df, min_unique=2):
    """Removes constant columns (all 0, all 1, or otherwise no variation)."""
    variable_cols = df.columns[df.nunique() >= min_unique]
    return df[variable_cols]


def filter_threshold(df, threshold=0.5):
    """Keeps only columns that exceed the threshold in at least one genome."""
    return df.loc[:, (df > threshold).any(axis=0)]


def filter_select(df, pathways=SELECTED_PATHWAYS):
    """Keeps only the pathways from the curated list (that exist in the df)."""
    present = [p for p in pathways if p in df.columns]
    missing = [p for p in pathways if p not in df.columns]
    if missing:
        print(f"[warning] {len(missing)} requested pathways not found in the TSV: {missing}")
    return df[present]


def filter_category_select(df, categories_to_keep, categories=CATEGORIES):
    """
    Keeps the INDIVIDUAL pathways (not aggregated) belonging to the
    categories chosen in `categories_to_keep`.
    """
    unknown = [c for c in categories_to_keep if c not in categories]
    if unknown:
        print(f"[warning] unrecognized categories (ignored): {unknown}")
        print(f"          available categories: {list(categories.keys())}")

    wanted_cols = []
    groups = []
    for cat in categories_to_keep:
        if cat not in categories:
            continue
        present_in_cat = []
        for col in categories[cat]:
            if col in df.columns and col not in wanted_cols:
                wanted_cols.append(col)
                present_in_cat.append(col)
            elif col not in df.columns:
                print(f"[warning] pathway '{col}' (category '{cat}') not found in the TSV")
        if present_in_cat:
            groups.append((cat, present_in_cat))

    if not wanted_cols:
        raise ValueError("No pathways found for the requested categories.")

    return df[wanted_cols], groups


def filter_variance_category(df, categories=CATEGORIES):
    """
    Filters pathways by variance first, then groups all remaining pathways
    into macro-categories. Pathways not mapped in CATEGORIES are grouped in 'Other'.
    """
    df_var = filter_variance(df)
    
    groups = []
    used_cols = set()

    for cat, cols in categories.items():
        present_in_cat = [c for c in cols if c in df_var.columns]
        if present_in_cat:
            groups.append((cat, present_in_cat))
            used_cols.update(present_in_cat)

    # Any pathway that passed variance filter but is not in CATEGORIES
    leftovers = [c for c in df_var.columns if c not in used_cols]
    if leftovers:
        groups.append(("Other", leftovers))

    wanted_cols = [c for _, cols in groups for c in cols]
    return df_var[wanted_cols], groups


def aggregate_categories(df, categories=CATEGORIES):
    """Aggregates pathways into macro-categories (mean per category)."""
    out = {}
    for cat, cols in categories.items():
        present = [c for c in cols if c in df.columns]
        if present:
            out[cat] = df[present].mean(axis=1)
    return pd.DataFrame(out)


def cluster_reorder(df):
    """Reorders rows and columns with Euclidean hierarchical clustering."""
    if len(df.index) < 2:
        return df
    row_linkage = hierarchy.linkage(distance.pdist(df.values), method="average")
    row_order = hierarchy.leaves_list(row_linkage)
    col_linkage = hierarchy.linkage(distance.pdist(df.T.values), method="average")
    col_order = hierarchy.leaves_list(col_linkage)
    return df.iloc[row_order, col_order]


def cluster_reorder_rows_only(df):
    """Like cluster_reorder but does NOT touch column order."""
    if len(df.index) < 2:
        return df
    row_linkage = hierarchy.linkage(distance.pdist(df.values), method="average")
    row_order = hierarchy.leaves_list(row_linkage)
    return df.iloc[row_order, :]


def plot_heatmap(df, outfile, figsize=None, dpi=300):
    sns.set_style("white")
    n_rows, n_cols = df.shape
    if figsize is None:
        figsize = (max(8, n_cols * 0.45), max(4, n_rows * 0.4))
    fig, ax = plt.subplots(figsize=figsize)
    ax = sns.heatmap(
        df, cmap=plt.cm.YlOrRd, linewidths=0.5, linecolor="grey",
        square=True, xticklabels=True, yticklabels=True,
        cbar_kws={"shrink": 0.6, "aspect": 15, "pad": 0.02},
        ax=ax,
    )
    ax.xaxis.tick_top()

    # Formattazione uniforme dei tick neri e ben definiti
    ax.tick_params(axis='x', which='both', top=True, bottom=False, length=4, width=1.0, color="black")
    ax.tick_params(axis='y', which='both', left=True, right=False, length=4, width=1.0, color="black")
    
    plt.setp(
        ax.get_xticklabels(),
        rotation=90,
        ha="center",
        va="bottom",
        fontsize=11
    )
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=10)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=11, length=4, pad=4)

    plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight", dpi=dpi)
    plt.close()


def plot_heatmap_grouped(df, groups, outfile, spacer_width=0.8, figsize=None, dpi=300):
    """
    Heatmap with pathways grouped by category, blocks separated by empty
    space, and a category label under each block.
    """
    import textwrap
    sns.set_style("white")

    col_inch = 0.45

    col_order = []
    col_positions = []   # x center of each real column
    group_bounds = []    # (x_start, x_end, category_name) for labels
    x = 0.0
    for gi, (cat, cols) in enumerate(groups):
        start_x = x
        for c in cols:
            col_order.append(c)
            col_positions.append(x + 0.5)
            x += 1.0
        group_bounds.append((start_x, x, cat))
        if gi < len(groups) - 1:
            x += spacer_width

    n_rows = df.shape[0]
    n_cols_total = x

    if figsize is None:
        figsize = (max(8, n_cols_total * col_inch), max(4, n_rows * 0.4) + 2.2)

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    data = df[col_order].values.astype(float)
    vmin, vmax = 0.0, 1.0
    cmap = plt.cm.YlOrRd

    for j, c in enumerate(col_order):
        x0 = col_positions[j] - 0.5
        col_vals = data[:, j]
        for i, val in enumerate(col_vals):
            color = cmap((val - vmin) / (vmax - vmin) if vmax > vmin else 0)
            rect = plt.Rectangle((x0, i), 1.0, 1.0, facecolor=color,
                                 edgecolor="grey", linewidth=0.5)
            ax.add_patch(rect)

    ax.set_xlim(0, n_cols_total)
    ax.set_ylim(n_rows, 0)
    ax.set_aspect("equal")

    # Rimuove la cornice esterna (spines)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Tick neri, spessi e ben visibili identici alla heatmap standard
    ax.tick_params(axis='x', which='both', top=True, bottom=False, length=4, width=1.0, color="black")
    ax.tick_params(axis='y', which='both', left=True, right=False, length=4, width=1.0, color="black")

    # Pathway ticks/labels in alto
    ax.set_xticks(col_positions)
    ax.set_xticklabels(
        col_order,
        rotation=90,
        ha="center",
        va="bottom",
        fontsize=11
    )
    ax.xaxis.tick_top()

    # Genome/sample ticks/labels a sinistra
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    ax.set_yticklabels(df.index, rotation=0, fontsize=10)

    # Category label under each block
    fontsize = 11
    char_width_in = fontsize * 0.62 / 72.0
    max_label_lines = 1
    for (x_start, x_end, cat) in group_bounds:
        mid = (x_start + x_end) / 2
        block_width_in = (x_end - x_start) * col_inch
        max_chars = max(6, int(block_width_in / char_width_in))

        if len(cat) > max_chars:
            label_text = "\n".join(textwrap.wrap(cat, width=max_chars))
            max_label_lines = max(max_label_lines, label_text.count("\n") + 1)
        else:
            label_text = cat

        ax.text(mid, n_rows + 0.5, label_text, ha="center", va="top",
                fontsize=fontsize, fontweight="bold", rotation=0, clip_on=False,
                linespacing=1.3)

    if max_label_lines > 1:
        fig.set_figheight(figsize[1] + 0.3 * (max_label_lines - 1))

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, aspect=15, pad=0.03)
    cbar.ax.tick_params(labelsize=11, length=4, pad=4)

    plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight", dpi=dpi)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Filter and visualize KEGG-decoder output for publication."
    )
    parser.add_argument("-i", "--input", required=True, help="Full TSV produced by KEGG-decoder")
    parser.add_argument("-o", "--output", required=True, help="Path for the filtered output TSV")
    parser.add_argument(
        "--mode", choices=["variance", "threshold", "select", "category", "category_select", "variance_category"],
        default="variance",
        help="Filtering strategy (default: variance)",
    )
    parser.add_argument("--threshold", type=float, default=0.5, help="Threshold for --mode threshold")
    parser.add_argument(
        "--categories", nargs="+",
        help="For --mode category_select: names of the categories to keep (see CATEGORIES in the script). "
             "Example: --categories \"Nitrogen cycle\" \"Sulfur cycle\"",
    )
    parser.add_argument("--no-cluster", action="store_true", help="Disable hierarchical clustering")
    parser.add_argument("--fig", help="Path for the heatmap image file (e.g. heatmap.png / heatmap.svg)")
    parser.add_argument(
        "--format", choices=["png", "svg"], default="png",
        help="Export format for the image (default: png)",
    )
    parser.add_argument(
        "--dpi", type=int, default=300,
        help="Image resolution in DPI for raster output (default: 300)",
    )
    args = parser.parse_args()

    df = load_table(args.input)
    print(f"Original table: {df.shape[0]} genomes x {df.shape[1]} pathways")

    groups = None  # used for grouped plotting modes

    if args.mode == "variance":
        df_filt = filter_variance(df)
    elif args.mode == "threshold":
        df_filt = filter_threshold(df, args.threshold)
    elif args.mode == "select":
        df_filt = filter_select(df)
    elif args.mode == "category":
        df_filt = aggregate_categories(df)
    elif args.mode == "category_select":
        if not args.categories:
            parser.error("--mode category_select requires --categories \"Name1\" \"Name2\" ...")
        df_filt, groups = filter_category_select(df, args.categories)
    elif args.mode == "variance_category":
        df_filt, groups = filter_variance_category(df)

    print(f"Filtered table ({args.mode}): {df_filt.shape[0]} genomes x {df_filt.shape[1]} pathways")

    if not args.no_cluster and df_filt.shape[0] >= 2:
        if args.mode in ["category_select", "variance_category"]:
            df_filt = cluster_reorder_rows_only(df_filt)
            if groups is not None:
                groups = [(cat, [c for c in cols if c in df_filt.columns]) for cat, cols in groups]
        else:
            df_filt = cluster_reorder(df_filt)

    df_filt.to_csv(args.output, sep="\t")
    print(f"Saved: {args.output}")

    if args.fig:
        base_fig_path, _ = os.path.splitext(args.fig)
        final_fig_path = f"{base_fig_path}.{args.format}"

        if args.mode in ["category_select", "variance_category"] and groups is not None:
            plot_heatmap_grouped(df_filt, groups, final_fig_path, dpi=args.dpi)
        else:
            plot_heatmap(df_filt, final_fig_path, dpi=args.dpi)
        print(f"Heatmap saved ({args.format.upper()}, {args.dpi} DPI): {final_fig_path}")


if __name__ == "__main__":
    main()