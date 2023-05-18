import pandas as pd
import numpy as np
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt

csv_file = sys.argv[1]
segment= sys.argv[2]
output_file = sys.argv[3]


#PRIMER LENGHTS IN DICT - IMPORT INSTEAD? # MANGLER PROBE FOR TRIPLEX
primer_dict ={
    "InfB_MP_BMf_38_Ward" : 23,
    "InfB_MP_BMr_132_antisenseB_Ward" : 23,
    "InfB_MP_BM-probe_69_Ward" : 30,
    "InfA_MP_H1_H3_M52c" : 21,
    "InfA_MP_H1_H3_M149r": 23,
    "InfA_MP_H1_H3_probeM93c" : 23,
    "InfB_HA_BHA_188F" : 21,
    "InfB_HA_BHA270R" : 23,
    "InfB_HA_probe_VIC2" : 27,
    "InfB_HA_probe_VIC2" : 27,
    "InfB_HA_probe_VIC2" : 27,
    "InfB_HA_probe_YAM2" : 27,
    "InfB_HA_probe_YAM2" : 27,
    "InfB_HA_probe_YAM2" : 27,
    "InfB_HA_probe_YAM2" : 27,
    "InfA_MA_H1_H3_InfA_For1" : 25,
    "InfA_MA_H1_H3_InfA_For1" : 25,
    "InfA_MA_H1_H3_InfA_For2" : 25,
    "InfA_MA_H1_H3_InfA_For2" : 25,
    "InfA_MA_H1_H3_InfA_For2" : 25,
    "InfA_MA_H1_H3_InfA_For2" : 25,
    "InfA_MA_H1_H3_InfA_Rev1" : 23,
    "InfA_MA_H1_H3_InfA_Rev1" : 23,
    "InfA_MA_H1_H3_InfA_Rev1" : 23,
    "InfA_MA_H1_H3_InfA_Rev1" : 23,
    "InfA_MA_H1_H3_InfA_Rev1" : 23,
    "InfA_MA_H1_H3_InfA_Rev1" : 23,
    "InfA_MA_H1_H3_InfA_Rev2" : 23,
    "InfB_MA_InfA_For" : 22,
    "InfB_MA_InfA_For" : 22,
    "InfB_MA_InfA_Rev" : 22,
    "INFA_HA_H1_FluSwH1R236" : 25,
    "INFA_HA_H1_FluSwH1R318" : 25,
    "INFA_HA_H1_FluSwH1TM292" : 28
}


#IMPORT PRIMER ALIGNMENT
df1 = pd.read_csv(csv_file,  sep='\t', names=["Primer", "SequenceID", "Raw_Match(%)", "Alignment_Length", "Mismatches", "Gaps", "Primer_Start", "Primer_End", "Seq_Start", "Seq_End", "Expect", "Score"], index_col=False)


def check_subtype(Subtype):
    if "H3" in Subtype:
        return "H3"
    elif "H1" in Subtype:
        return "H1"
    elif "Bvic" in Subtype:
        return "Bvic"
    elif "Byam" in Subtype:
        return "Byam"
    else:
        return "Unknown subtype"

Discard = []
Discard.append(check_subtype(segment))


#DATA MINGELING NAD BASIC CALCULATIONS

df1 = df1.iloc[:,0:12]
df1['Primer_Length'] = df1['Primer'].map(primer_dict)
df1["Primer_Length/Alignment_length_dif"] = df1["Primer_Length"] - df1["Alignment_Length"]
df1["Primer_Length/Alignment_length_dif"] = df1["Primer_Length/Alignment_length_dif"].abs()
df1["ID"] = df1['Primer'].astype(str) +"_"+ df1["SequenceID"].astype(str)
df1["Length_Para"] = df1["Primer_Length"] - df1["Alignment_Length"]

df1 = df1.sort_values("Score", ascending=False).drop_duplicates('ID').sort_index()

df1["Match(%)"] = ((df1["Primer_Length"] - (df1["Mismatches"] + df1["Gaps"] + ((df1["Primer_Length"]) - df1["Alignment_Length"]))) / df1["Primer_Length"])
df1["Total_Mis"] = df1["Mismatches"] + df1["Gaps"] + df1["Primer_Length/Alignment_length_dif"]

mask = df1['Primer'].str.contains('|'.join(Discard))
df1 = df1.loc[mask]


#SAVE FINAL DATAFRAME TO CSV

df1.to_csv("{}_primer_alignment.csv".format(output_file), index=False)

# set SequenceID as the index
df1 = df1.set_index("SequenceID")

# pivot the dataframe
df_pivot = df1.pivot(columns="Primer", values="Total_Mis") 


# create the plot
plt.figure(figsize=(40, 20))
sns.heatmap(df_pivot.astype(float).fillna(np.nan), cmap="rocket_r", annot=True, fmt=".1f", linewidths=0.5, annot_kws={"size": 30}, square=True)
# set the axis labels
plt.xlabel("Primer", fontsize=20)
plt.ylabel("SequenceID", fontsize=20) 

# increase the fontsize of the xticks and yticks
plt.xticks(fontsize=20, rotation=90)
plt.yticks(fontsize=20, rotation=0)

# show the plot
plt.savefig("{}_primer_alignment.png".format(output_file), dpi=300)