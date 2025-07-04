import pandas as pd

TAX = pd.read_csv("Datasets/tax.csv")
CULTURE = pd.read_csv("Datasets/cultureProportions_Env.csv")
MICROBIOME = pd.read_csv("Datasets/microbiomeProportions_Env.csv")

def find_n_sequences(bacteria, num_sequences):
    culture = CULTURE.set_index("OTU")
    bacteria_row = culture.loc[bacteria]
    print(bacteria_row)
    print(bacteria_row.max())

    return []

def main():
    bacteria = input("Type a bacteria: ")
    num_sequences = input("Type a number of environments you want to grow it in: ")

    print(find_n_sequences(bacteria, num_sequences))
    

if __name__ == "__main__":
    main()