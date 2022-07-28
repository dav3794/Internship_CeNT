import pandas as pd
import numpy as np
import glob
import os

def merge(family_id):
    df_0 = pd.read_csv(f"{family_id}_data.csv", sep=",", header=0)
    df_0 = df_0.to_numpy()

    files = os.path.join(".", "****_homologues_knots.csv")
    files = glob.glob(files)

    df_ext= pd.concat(map(lambda x: pd.read_csv(x, sep=","), files), ignore_index=True)
    df_ext = df_ext.drop_duplicates('sequence').to_numpy()

    for id, seq, knot in df_ext:
        df_0 = np.append(df_0, [[id, "", knot, seq, ""]], axis=0)
    print(df_0)

    pd.DataFrame(df_0).drop_duplicates(3).to_csv(f"{family_id}_extended.csv", header=["ID", 'chain', 'knotted', "sequence", 'loc_knot'], index=False)

merge("PF01699")
