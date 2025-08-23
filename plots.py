# plot_composition_by_order.py
import os, re
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# === INPUTS ===
TAX_CSV = "Datasets/tax.csv"  # must contain columns: OTU, order (case-insensitive ok)
CSV_FILES = [
    'testing_outputs/01_Otu52_5_Sheep_Blood_Ery.csv',
    'testing_outputs/02_Otu13_MiPro_Kan.csv',
    'testing_outputs/03_Otu11_BHI_Cip.csv',
    'testing_outputs/04_Otu24_5_Sheep_Blood_Kan_Mucus_Amp.csv',
    'testing_outputs/05_Otu51_mGAM_Kan_BMM_raffinose_NA.csv',
    'testing_outputs/06_Otu67_Yeast_Casitone_NA_Postgate_Cip.csv',
    'testing_outputs/07_Otu23_mGAM_Cip_mGAM_Kan_MiPro_NA.csv',
    'testing_outputs/08_Otu39_BMM_Inulin_NA_5_Sheep_Blood_Kan_mGAM_NA.csv',
]

# === HELPERS ===
def parse_main_otu(path: str) -> str:
    m = re.search(r'(Otu\d+)', os.path.basename(path), flags=re.IGNORECASE)
    return m.group(1) if m else None

def growth_labels(cols):
    # If first col is 'input' (any case), label Input then Growth 1..N-1
    if len(cols) == 0:
        return []
    first = cols[0].strip().lower()
    if first == "input":
        return ["Input"] + [f"Growth {i}" for i in range(1, len(cols))]
    # Otherwise, just number all stages for clarity
    return [f"Stage {i+1}" for i in range(len(cols))]

def get_taxonomy_table(path: str) -> pd.DataFrame:
    tax = pd.read_csv(path)
    cols = {c.lower(): c for c in tax.columns}
    if 'otu' not in cols:
        raise ValueError("tax.csv must have an 'OTU' column.")
    if 'order' not in cols:
        raise ValueError("tax.csv must have an 'order' column.")
    tax = tax.rename(columns={cols['otu']: 'OTU', cols['order']: 'order'})
    tax['order'] = tax['order'].fillna('Unknown').replace('', 'Unknown')
    return tax[['OTU', 'order']]

def collect_all_orders(tax: pd.DataFrame, csv_paths: list[str]) -> list[str]:
    orders = set()
    otu_set = set()
    for p in csv_paths:
        df = pd.read_csv(p)
        otu_set |= set(df['OTU'].astype(str))
    orders |= set(tax.loc[tax['OTU'].astype(str).isin(otu_set), 'order'])
    # Always include Unknown if any OTUs not present in tax
    if len(otu_set - set(tax['OTU'].astype(str))) > 0:
        orders.add('Unknown')
    return sorted(orders)

def build_order_palette(order_list: list[str]):
    # Stable color mapping for orders using matplotlib tab20
    cmap = plt.get_cmap("tab20")
    colors = {}
    for i, ord_name in enumerate(order_list):
        colors[ord_name] = cmap(i % 20)
    return colors

def plot_one(csv_path: str, ax_top, ax_bottom, tax: pd.DataFrame, order_colors: dict):
    df = pd.read_csv(csv_path)
    if 'OTU' not in df.columns:
        raise ValueError(f"'OTU' column missing in {csv_path}")

    # numeric-ify conditions (leave OTU alone)
    conditions = [c for c in df.columns if c != 'OTU']
    df[conditions] = df[conditions].apply(pd.to_numeric, errors='coerce').fillna(0.0)

    # Focal OTU
    main_otu = parse_main_otu(csv_path)
    if main_otu not in set(df['OTU'].astype(str)):
        # fallback to most abundant by total across conditions
        row_idx = df[conditions].sum(axis=1).idxmax()
        main_otu = str(df.loc[row_idx, 'OTU'])

    # ----- TOP LINE: focal OTU over stages -----
    x = np.arange(len(conditions))
    y = df.loc[df['OTU'].astype(str) == main_otu, conditions].iloc[0].values.astype(float)
    ax_top.plot(x, y, marker='o', linewidth=1.6)
    ax_top.set_title(main_otu, fontsize=10, pad=1)
    ax_top.set_ylabel('Rel.\nabd', rotation=0, labelpad=18, va='center')
    ax_top.set_ylim(0, 1)           # y range from 0 to 1
    ax_top.set_yticks([0, 0.5, 1])  # ticks at 0, 0.5, and 1
    ax_top.set_xticks([])

    # ----- BOTTOM STACK: aggregated by ORDER -----
    # Merge order onto per-OTU rows
    merged = df.merge(tax, on='OTU', how='left')
    merged['order'] = merged['order'].fillna('Unknown').replace('', 'Unknown')

    # Sum by order for each condition
    by_order = (
        merged.groupby('order', as_index=True)[conditions]
        .sum()
        .reindex(order_colors.keys(), fill_value=0.0)  # ensure consistent order/color presence
    )

    # Determine focal order
    focal_order = merged.loc[merged['OTU'].astype(str) == main_otu, 'order']
    focal_order = focal_order.iloc[0] if len(focal_order) else 'Unknown'

    # Optional: order stacking by total at final stage with focal order last to get border fully visible
    final_col = conditions[-1]
    order_sequence = list(by_order.sort_values(final_col).index)
    if focal_order in order_sequence:
        order_sequence.remove(focal_order)
        order_sequence.append(focal_order)

    bottoms = np.zeros(len(conditions), dtype=float)
    for ord_name in order_sequence:
        vals = by_order.loc[ord_name, conditions].values.astype(float)
        if ord_name == focal_order:
            ax_bottom.bar(x, vals, bottom=bottoms, color=order_colors[ord_name],
                          edgecolor='black', linewidth=1.8)
        else:
            ax_bottom.bar(x, vals, bottom=bottoms, color=order_colors[ord_name])
        bottoms += vals

    ax_bottom.set_ylabel('Rel.\nabd', rotation=0, labelpad=18, va='center')
    ax_bottom.set_xticks(x)
    ax_bottom.set_xticklabels(growth_labels(conditions), rotation=45, ha='right')

def make_grid(csv_paths, tax_csv=TAX_CSV, out_path='compositional_by_order.pdf'):
    tax = get_taxonomy_table(tax_csv)
    all_orders = collect_all_orders(tax, csv_paths)
    order_colors = build_order_palette(all_orders)

    fig = plt.figure(figsize=(16, 12))
    outer = fig.add_gridspec(3, 3, wspace=0.35, hspace=0.5)

    # Positions for 8 plots (skip bottom-right for legend)
    positions = [(0,0), (0,1), (0,2),
                 (1,0), (1,1), (1,2),
                 (2,0), (2,1)]

    for i, (r, c) in enumerate(positions):
        csv_path = csv_paths[i]
        sub = outer[r, c].subgridspec(2, 1, height_ratios=[1, 2], hspace=0.08)
        ax_t = fig.add_subplot(sub[0, 0])
        ax_b = fig.add_subplot(sub[1, 0])
        plot_one(csv_path, ax_t, ax_b, tax, order_colors)

    # Bottom-right legend/info panel (the 9th slot)
    ax_leg = fig.add_subplot(outer[2, 2])
    ax_leg.axis('off')
    ax_leg.set_title('Order Legend', fontsize=11, pad=6)

    patches = [Patch(facecolor=order_colors[o], edgecolor='none', label=o) for o in all_orders]

    # Force multi-column legend inside this subplot
    ncol = 2   # or 3 depending on how wide you want
    ax_leg.legend(
        handles=patches,
        loc='center',
        frameon=False,
        ncol=ncol,
        fontsize=9,
        handlelength=1.2,
        columnspacing=0.8
    )

    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path

def plot_log10_histograms_subplots(save_prefix="fold_enrichment_log10"):
    # Filter only 1, 2, 3 culture enrichments & remove zeros
    filtered_data = MASTER_DATA[MASTER_DATA['num_growth_conditions'].isin([1, 2, 3])]
    filtered_data = filtered_data[filtered_data['fold_change'] > 0]

    colors = {1: "red", 2: "red", 3: "red"}

    # Global limits for consistency
    log_vals_all = np.log10(filtered_data['fold_change'])
    global_xmin, global_xmax = log_vals_all.min(), log_vals_all.max()

    # Bin width 0.15 for all
    bin_width = 0.15
    bins = np.arange(global_xmin, global_xmax + bin_width, bin_width)

    # Global y-limit for counts
    global_ymax = 0
    for condition in [1, 2, 3]:
        subset = filtered_data[filtered_data['num_growth_conditions'] == condition]
        log_vals = np.log10(subset['fold_change'])
        counts, _ = np.histogram(log_vals, bins=bins)  # counts, not density
        global_ymax = max(global_ymax, counts.max())

    # Create subplots — reduced width by 20%
    fig, axes = plt.subplots(1, 3, figsize=(5.5, 2.25), sharex=True, sharey=True)

    for idx, condition in enumerate([1, 2, 3]):
        ax = axes[idx]
        subset = filtered_data[filtered_data['num_growth_conditions'] == condition]
        log_vals = np.log10(subset['fold_change'])

        # Histogram with black borders (counts)
        ax.hist(log_vals, bins=bins, color=colors[condition], alpha=0.7,
                density=False, edgecolor="black", linewidth=0.5)

        # Labels & title
        ax.set_xlabel("log10(Fold Change)")
        ax.set_title(f"{condition} Culture Enrichment", fontsize=9)
        if idx == 0:
            ax.set_ylabel("Number of OTUs")

        # Horizontal grid only
        ax.grid(axis="y", linestyle="--", alpha=0.6)

        # Apply consistent limits
        ax.set_xlim(global_xmin, global_xmax)
        ax.set_ylim(0, global_ymax * 1.05)

    plt.tight_layout()

    # Save in both PDF (vector) and PNG (raster)
    plt.savefig(f"{save_prefix}_subplots.pdf", bbox_inches="tight")
    plt.savefig(f"{save_prefix}_subplots.png", dpi=300, bbox_inches="tight")
    plt.show()


def main():
    # plot_log10_histograms_subplots()
    make_grid(CSV_FILES, tax_csv=TAX_CSV, out_path='compositional_by_order.pdf')



if __name__ == "__main__":
    main()