# analysis/scripts/merge_logs.py
import pandas as pd
dfs = [pd.read_parquet(f) for f in snakemake.input]
combined = pd.concat(dfs, ignore_index=True)
combined.to_parquet(snakemake.output[0])