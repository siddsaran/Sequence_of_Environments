import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import main

MASTER_DATA = pd.read_csv('allChanges.csv')

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
    plot_log10_histograms_subplots()


if __name__ == "__main__":
    main()