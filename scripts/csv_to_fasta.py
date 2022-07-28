import pandas as pd
import numpy as np
from math import ceil
import re

def convert_to_fasta(file):
    data = pd.read_csv(f"{file}.csv", sep=",", header=0)
    data = data.to_numpy()
    output = ""

    for line in data:
        t = map(str, line[0:3])
        header = ">"+",".join(t)
        #header += str(line[4])
        output += (header+"\n")
        seq = line[3]
        size = len(seq)
        for i in 64*np.arange(ceil(size/64)):
            if i+64<size:
                output += (seq[i:i+64]+"\n")
            else:
                output += (seq[i:]+"\n")

    with open(f"{file}.fasta", "w") as output_file:
        output_file.write(output)

def convert_to_csv(file):
    header = ['id', 'sequence']
    data = []
    seq_temp = ""
    with open(f"{file}.fasta", "r") as input_file:
        for line in input_file.readlines():
            if line[0] == ">":
                if len(seq_temp) > 0:
                    data.append([entry, seq_temp])
                seq_temp = ""
                pattern = 'AFDB:AF-.*?-F1'
                match_results = re.search(pattern, line, re.IGNORECASE)
                entry = match_results.group()
                entry = re.sub("AFDB:AF-", "", entry)
                entry = re.sub("-F1", "", entry)
                #print(entry)

            else:
                seq_temp += line.strip()
        else:
            data.append([entry, seq_temp])
            print("Done")
    #print(data)
    pd.DataFrame(data).drop_duplicates(1).to_csv(f"{file}.txt", sep=",", header=header, index=False)

#convert_to_fasta("PF01699_data")
convert_to_csv("4kpp_homologues")
