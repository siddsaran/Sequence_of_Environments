import pandas as pd
from itertools import permutations

TAX = pd.read_csv("Datasets/tax.csv")
CULTURE = pd.read_csv("Datasets/cultureProportions_Env.csv")
MICROBIOME = pd.read_csv("Datasets/microbiomeProportions_Env.csv")
LIMIT = 0.00001

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
    # for each combo, calculate the relative abundance of the bacteria
    for combo in combos:
        rel_abd = 0
        rolling_sum = 0
        for env in combo:
            rolling_sum += culture[env].sum()
            rel_abd = bacteria_row[env] / rolling_sum
        combo_dict[combo] = rel_abd
    # return max abundance
    max_key = max(combo_dict, key=combo_dict.get)
    max_value = combo_dict[max_key]
    return max_key, max_value


def main():
    bacteria = input("Type a bacteria: ")
    num_sequences = int(input("Type the number of environments you want to grow: "))
    print(find_n_sequences(bacteria, num_sequences))
    

if __name__ == "__main__":
    main()