import pandas as pd

TAX = pd.read_csv("Datasets/tax.csv")
CULTURE = pd.read_csv("Datasets/cultureProportions_Env.csv")
MICROBIOME = pd.read_csv("Datasets/microbiomeProportions_Env.csv")

def find_next_environment(culture, bacteria):
    bacteria_row = culture.loc[bacteria]
    return bacteria_row.idxmax()


def find_n_sequences(bacteria, num_sequences):
    environments = []
    culture = CULTURE.set_index("OTU")
    for _ in range(int(num_sequences)):
        curr_env = find_next_environment(culture, bacteria)
        environments.append(curr_env)
        culture = culture.drop(columns=curr_env)
    return environments

def main():
    bacteria = input("Type a bacteria: ")
    num_sequences = input("Type the number of environments you want to grow: ")

    print(find_n_sequences(bacteria, num_sequences))
    

if __name__ == "__main__":
    main()