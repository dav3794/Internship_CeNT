from urllib.request import urlopen
import pandas as pd
import numpy as np
import wget
import re
import os

def load_knotted_OTC():
    try:
        wget.download("https://knotprot.cent.uw.edu.pl/browse/?pfam=OTCace&set=True&bridgeType=&array=0&raw=1", out="raw_OTC1.csv")
        wget.download("https://knotprot.cent.uw.edu.pl/browse/?pfam=OTCace_N&set=True&bridgeType=&array=0&raw=1", out="raw_OTC2.csv")

        knotted_OTC = pd.concat(map(lambda x: pd.read_csv(x, sep=";"), ["raw_OTC1.csv","raw_OTC2.csv"]), ignore_index=True)
        os.remove("raw_OTC1.csv")
        os.remove("raw_OTC2.csv")
    except:
        print("Error")
        return 0

    knotted_OTC = knotted_OTC.drop_duplicates()
    knotted_OTC = knotted_OTC.to_numpy()
    size = len(knotted_OTC)
    print(knotted_OTC)
    return knotted_OTC

def load_ids(pfam):
    df_knots = pd.read_csv(f"{pfam}_knots.csv", sep=",", header=0)
    df_knots = df_knots.drop_duplicates()
    df_knots = df_knots.to_numpy()

    size = len(df_knots)
    df_knots = np.insert(df_knots, 3, np.full(size, ''), axis=1)
    df_knots = np.insert(df_knots, 4, np.full(size, ''), axis=1)
    print(df_knots)
    return df_knots

def get_seq(knotted_file, drop_dupl_seq=True):
    size = len(knotted_file)
    id_to_del = []

    for i, row in enumerate(knotted_file):
        print(f"{i+1} / {size}")
        id = row[0]
        chain = row[1]
        knot = row[2]

        try:
            url = f"https://genus.fuw.edu.pl/view/{id}/{chain}/"
            page = urlopen(url)
            html = page.read().decode("utf-8")

            seq_pattern = '<div class="sequence" id="sequence">.*?</div>'
            match_results = re.search(seq_pattern, html, re.IGNORECASE)
            seq = match_results.group()
            seq = re.sub("<.*?>", "", seq)
            knotted_file[i][3] = seq

            if knot:
                url = f"https://knotprot.cent.uw.edu.pl/view/{id}/{chain}/"
                page = urlopen(url)
                html = page.read().decode("utf-8")

                pattern_knotcore = '\/colour_sequence.*?;'
                subpattern = '"[0-9]+\-[0-9]+"'
                match_results = re.search(pattern_knotcore, html, re.IGNORECASE)
                loc_knot = match_results.group()
                match_results = re.search(subpattern, loc_knot, re.IGNORECASE)
                loc_knot = match_results.group()
                loc_knot = re.sub('"', "", loc_knot)
                knotted_file[i][4] = loc_knot

            #lower, upper = map(int, loc_knot.split('-'))
            #knot_seq = seq[lower-1 : upper]
            #knotted_file[i][4] = knot_seq

        except:
            print("Error for id:", id, chain)
            id_to_del.append(i)

    knotted_file = np.delete(knotted_file, id_to_del, axis=0)
    if drop_dupl_seq:
        knotted_file = pd.DataFrame(knotted_file).drop_duplicates([3])

    return knotted_file

def save_data(df, pfam):
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    header_csv = ['pdb_id', 'chain', "knotted", 'sequence', 'knot_loc']
    df.to_csv(f"{pfam}_data.csv", header=header_csv, index=False)
    print("Data saved.")

#data = load_knotted_OTC()
pfam = "PF01699"
data = load_ids(pfam)
data = get_seq(data)
save_data(data, pfam)
