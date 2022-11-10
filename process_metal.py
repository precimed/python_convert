import pandas as pd
import argparse
import os

# Input: metal output with UID in the MarkerName column.
# Add RSID from HRC reference panel and CHR and BP based on UID.

# Parse arguments --------------------------
parser = argparse.ArgumentParser(description="Process output of metal assuming MarkerName column contain UID. Add RSID from HRC reference panel and CHR and BP based on UID.")
parser.add_argument("--fname", required=True, help="Path to input metal sumstats file.")
parser.add_argument("--hrc", default="/cluster/projects/p33/users/alexeas/hrc/HRC.r1-1.GRCh37.wgs.mac5.sites.uid.tab.gz", help="Path to HRC reference file.")
parser.add_argument("--out", help="Output file name.")
args = parser.parse_args()
if not args.out:
    name, ext = os.path.splitext(args.fname)
    ext = ext if ext.endswith("gz") else f"{ext}.gz"
    args.out = f"{name}.processed{ext}"

fname = args.fname
df = pd.read_csv(fname, sep='\t')

hrc_fname = args.hrc
hrc = pd.read_csv(hrc_fname, sep='\t', usecols=["ID","UID"])
hrc.rename(columns={"ID":"RSID"}, inplace=True)
print(f"{hrc.shape[0]} variants in HRC reference.")

df = df.merge(hrc, left_on="MarkerName", right_on="UID", how='left')
df["CHR"] = [uid.split(":")[0] for uid in df.MarkerName]
df["BP"] = [uid.split(":")[1] for uid in df.MarkerName]
df.drop(columns=["UID"], inplace=True)
print(f"{df.shape[0]} variants in {fname}")
df.to_csv(args.out, sep='\t', na_rep="NA", index=False)

