import pandas as pd
from itertools import permutations

TAX = pd.read_csv("Datasets/tax.csv")
CULTURE = pd.read_csv("Datasets/cultureProportions_Human (1).csv")
MICROBIOME = pd.read_csv("Datasets/microbiomeProportions_Human.csv")
LIMIT = .01

def find_cols_to_keep(culture, bacteria):
    cols_to_keep = []
    bacteria_row = culture.loc[bacteria]
    for col in culture.columns:
        if bacteria_row[col] > LIMIT:
            cols_to_keep.append(col)
    return cols_to_keep

def find_n_sequences(bacteria, num_sequences):
    culture = CULTURE.set_index("OTU")
    cols_to_keep = find_cols_to_keep(culture, bacteria)
    culture = culture[cols_to_keep]
    # generate all possible combinations of n environments
    all_combos = list(permutations(culture.columns, num_sequences))
    return find_best_sequence(culture, bacteria, all_combos)

def find_best_sequence(culture, bacteria, combos):
    combo_dict = dict()
    bacteria_row = culture.loc[bacteria]
    max_value = -1
    key = 0
    max_key = 0
    # for each combo, calculate the relative abundance of the bacteria
    for combo in combos:
        combo = list(combo)
        filtered = culture[combo].copy()
        env_looked = []
        prev_label = None
        for env in combo:
            env_looked.append(env)
            label = ",".join(env_looked) + " normalized"
            if prev_label is None:
                filtered[label] = filtered[env] / filtered[env].sum()
            else:
                filtered[label] = filtered[prev_label] * filtered[env]
                filtered[label] /= filtered[label].sum()
            prev_label = label
        combo_dict[key] = filtered
        val = filtered.loc[bacteria, prev_label]
        if val > max_value:
            max_value = val
            max_key = key
        key += 1
    # return max abundance
    return combo_dict[max_key]


def main():
    culture = CULTURE.set_index("OTU")
    main_data = pd.DataFrame(columns=[
        'OTU',
        'num_growth_conditions',
        'original_abundance',
        'final_abundance',
        'fold_change',
        'conditions'
    ])

    # Progress var
    num_done = 1

    otus_to_skip = ["Otu3"]
    otus_done = set()
    for otu in MICROBIOME["OTU"]:
        bacteria = otu
        cols = find_cols_to_keep(culture, bacteria)
        if otu in otus_to_skip:
            continue
        if otu in otus_done:
            continue
        otus_done.add(otu)
        i = 4
        if len(cols) < 3:
            i = len(cols) + 1
        for num_sequence in range(1, i):
            num_sequences = num_sequence
            df = find_n_sequences(bacteria, num_sequences)
            merged = MICROBIOME.merge(df, on='OTU', how='left')
            merged = merged.fillna(0)
            normalized_cols = list(df.columns[num_sequences:])
            original = round(MICROBIOME.loc[MICROBIOME['OTU'] == bacteria, 'relabd'].values[0] * 100, 5)
            abd = round(merged.loc[merged['OTU'] == bacteria, normalized_cols[-1]].values[0] * 100,5)
            factor = abd / original
            new_row = {
                'OTU': bacteria,
                'num_growth_conditions': num_sequences,
                'original_abundance': str(original) + "%",
                'final_abundance': str(abd) + "%",
                'fold_change': factor,
                'conditions': normalized_cols[-1]
            }
            print(f"Adding {bacteria} with {num_sequence} growth conditions, {(num_done/137) * 100:.3f}% done!")
            main_data.loc[len(main_data) - 2] = new_row
        num_done += 1

    main_data.to_csv('allChanges_inputscaled.csv')
    # bacteria = input("Type a bacteria: ")
    # num_sequences = int(input("Type the number of environments you want to grow: "))
    # df = find_n_sequences(bacteria, num_sequences)
    # merged = MICROBIOME.merge(df, on='OTU', how='left')
    # merged = merged.fillna(0)
    # normalized_cols = list(df.columns[num_sequences:])
    # original = MICROBIOME.loc[MICROBIOME['OTU'] == bacteria, 'relabd'].values[0]
    #
    # for col in normalized_cols:
    #     print(df[col])
    # abd = merged.loc[merged['OTU'] == bacteria, normalized_cols[-1]].values[0]
    # print(f"Original Relative abundance of {bacteria} is {(original*100):.4g}%")
    # print(f"Final Relative abundance of {bacteria} is {(abd*100):.4g}%")
    # print(f"Enrichment Factor is  {float(abd / original):.3g}")




    

if __name__ == "__main__":
    main()