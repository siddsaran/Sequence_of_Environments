from math import inf
import pandas as pd

TAX = pd.read_csv("Datasets/tax.csv")
CULTURE = pd.read_csv("Datasets/cultureProportions_Env.csv")
MICROBIOME = pd.read_csv("Datasets/microbiomeProportions_Env.csv")

def find_initial_environment(culture, bacteria):
    bacteria_row = culture.loc[bacteria]
    initial_max = bacteria_row.max()
    mask = bacteria_row == initial_max
    highest_initial = bacteria_row[mask]
    inital_culture = None
    if len(highest_initial) > 1:
        max_temp = float(inf)
        for name in highest_initial.index:
            sum = culture[name].sum() - initial_max
            if sum < max_temp:
                initial_culture = name
                max_temp = sum
    else:
        initial_culture = bacteria_row.idxmax()
    return initial_culture


def find_n_sequences(bacteria, num_sequences):
    environments = []
    culture = CULTURE.set_index("OTU")
    environments.append(find_initial_environment(culture, bacteria))
    return environments

def main():
    bacteria = input("Type a bacteria: ")
    num_sequences = input("Type the number of environments you want to grow ")

    print(find_n_sequences(bacteria, num_sequences))
    

if __name__ == "__main__":
    main()