from topoly import alexander
from urllib import request
from tqdm import tqdm
import os
import pandas as pd
import wget


def download_structures(list_of_chains):
    if not os.path.exists("./data/PDB_structures"):
        os.mkdir("./data/PDB_structures")
    for id in tqdm(list_of_chains):
        url = f"https://files.rcsb.org/download/{id}.pdb.gz"
        request.urlretrieve(url, f"./data/PDB_structures/{id}.pdb.gz")


def simple_knot_calculation_pdb(family_id):
    file_h = open("./data/pdb_pfam_mapping.txt")
    file = file_h.readlines()
    file_h.close()
    file = [[i.split("\t")[i0] for i0 in [0, 1, 4]] for i in file[2:]]
    data = []

    list_of_pdbs = [(i[0],i[1]) for i in file if i[2] == family_id]
    #download_structures(list_of_pdbs)

    for id, chain in tqdm(list_of_pdbs):
        structure = f"./data/PDB_structures/{id}.pdb.gz"
        knot = alexander(structure)
        if "0_1" in knot:
            if knot["0_1"] > 0.5:
                print("unknot")
                data.append([id,chain, 0])
            else:
                tops = list(knot.keys())
                tops.sort(key=lambda i: knot[i], reverse=True)
                print("knot: " + tops[0])
                data.append([id,chain, 1])

        print(knot)
    pd.DataFrame(data).to_csv(f"{family_id}_knots.csv", header=["entry id", "chain", 'knotted'], index=False)


def simple_knot_calculation_AF(filename):
    df = pd.read_csv(f"{filename}_homologues.txt", sep=",", header=0)
    df = df.to_numpy()
    data = []

    for id, sequence in tqdm(df):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v3.cif"
        wget.download(url, out=f"./data/{id}.cif")

        knot = alexander(f"./data/{id}.cif")
        if "0_1" in knot:
            if knot["0_1"] > 0.5:
                print("unknot")
                data.append([id,sequence, 0])
            else:
                tops = list(knot.keys())
                tops.sort(key=lambda i: knot[i], reverse=True)
                print("knot: " + tops[0])
                data.append([id,sequence, 1])

        print(knot)
    pd.DataFrame(data).to_csv(f"{filename}_homologues_knots.csv", header=["entry id", "sequence", 'knotted'], index=False)

if __name__ == '__main__':
    #simple_knot_calculation("PF01699")
    simple_knot_calculation_AF("5hwx")
