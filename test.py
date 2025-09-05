"""
Build per-scenario OTU dataframes from:
- microbiomeProportions_Human.csv
- cultureProportions_Human.csv
- testingset.csv

For each row in testingset.csv:
  - Start from "input" = mean relabd across all "Starting sample" rows per OTU
  - Apply each growth condition in order by multiplying by the matching column in cultureProportions
  - Renormalize after each step to keep relative abundances
  - Save a CSV with columns: OTU, input, A, A + B, A + B + C, ...
"""

import argparse
import re
from pathlib import Path
import pandas as pd


def clean_condition_label(label: str) -> str:
    """Strip quotes/whitespace from a condition label so it matches a CULTURES column name."""
    s = str(label).strip().strip('"').strip("'").strip()
    return s


def load_data(microbiome_csv: Path, cultures_csv: Path, testing_csv: Path):
    MICROBIOME = pd.read_csv(microbiome_csv)
    CULTURES = pd.read_csv(cultures_csv)
    TEST = pd.read_csv(testing_csv)
    return MICROBIOME, CULTURES, TEST


def compute_input_distribution(MICROBIOME: pd.DataFrame) -> pd.DataFrame:
    """
    Returns dataframe with columns ['OTU', 'input'] where 'input' is the mean relabd
    across rows whose 'Media' contains 'Starting sample' for that OTU.
    Missing OTUs get input=0.0 later when merged.
    """
    base = (
        MICROBIOME[MICROBIOME['Media'].str.contains('Starting sample', na=False)]
        .groupby('OTU', as_index=False)['relabd'].mean()
        .rename(columns={'relabd': 'input'})
    )
    return base


def cumulative_profiles(universe: pd.DataFrame, CULTURES: pd.DataFrame, conditions: list[str]) -> pd.DataFrame:
    """
    Given the 'universe' (['OTU','input']) and CULTURES table (with OTU + condition columns),
    build cumulative columns for the ordered list of conditions.
    """
    df = universe[['OTU', 'input']].copy()
    current = df[['OTU', 'input']].rename(columns={'input': 'abundance'}).copy()

    running_names = []
    for cond in conditions:
        cname = clean_condition_label(cond)
        running_names.append(cname)
        colname = ", ".join(running_names)
        if cname not in CULTURES.columns:
            raise KeyError(f"Condition '{cname}' not found in cultureProportions columns.")
        merged = pd.merge(current, CULTURES[['OTU', cname]], on='OTU', how='left')
        merged[cname] = merged[cname].fillna(0.0)
        merged['abundance'] = merged['abundance'] * merged[cname]
        total = merged['abundance'].sum()
        if total > 0:
            merged['abundance'] = merged['abundance'] / total
        step = merged[['OTU', 'abundance']].rename(columns={'abundance': colname})
        df = pd.merge(df, step, on='OTU', how='left')
        current = merged[['OTU', 'abundance']].copy()
    return df


def build_outputs(microbiome_csv: Path, cultures_csv: Path, testing_csv: Path, out_dir: Path) -> list[tuple[str, list[str], str]]:
    MICROBIOME, CULTURES, TEST = load_data(microbiome_csv, cultures_csv, testing_csv)
    base = compute_input_distribution(MICROBIOME)

    # Ensure we have one row per OTU found in CULTURES (the operating universe)
    universe = pd.merge(CULTURES[['OTU']], base, on='OTU', how='left')
    universe['input'] = universe['input'].fillna(0.0)

    out_dir.mkdir(parents=True, exist_ok=True)
    manifests: list[tuple[str, list[str], str]] = []

    for idx, row in TEST.iterrows():
        focal_otu = row['OTU']
        cond_field = str(row['conditions'])
        conditions = [c.strip() for c in cond_field.split(',') if str(c).strip() != ""]
        df = cumulative_profiles(universe, CULTURES, conditions)

        # Optional: sort OTUs alphabetically (or place focal OTU first — uncomment to place first)
        # df = df.set_index('OTU')
        # df = pd.concat([df.loc[[focal_otu]], df.drop(index=focal_otu)], axis=0).reset_index()

        # Safe filename
        slug = re.sub(r'[^A-Za-z0-9]+', '_', "__".join(conditions))[:120].strip("_")
        fname = f"{idx+1:02d}_{focal_otu}_{slug}.csv"
        fpath = out_dir / fname
        df.to_csv(fpath, index=False)
        manifests.append((focal_otu, conditions, str(fpath)))

    return manifests


def main():
    p = argparse.ArgumentParser(description="Build per-scenario OTU distributions from testing set.")
    p.add_argument("--microbiome", type=Path, default=Path("Datasets/microbiomeProportions_Human.csv"))
    p.add_argument("--cultures", type=Path, default=Path("Datasets/cultureProportions_Human (1).csv"))
    p.add_argument("--testing", type=Path, default=Path("testingset2.csv"))
    p.add_argument("--outdir", type=Path, default=Path("testing_outputs2"))
    args = p.parse_args()

    manifests = build_outputs(args.microbiome, args.cultures, args.testing, args.outdir)

    print(f"Wrote {len(manifests)} scenario files to '{args.outdir.resolve()}'.")
    for focal_otu, conds, path in manifests:
        print(f"- {focal_otu} | {', '.join(conds)} -> {path}")

if __name__ == "__main__":
    main()