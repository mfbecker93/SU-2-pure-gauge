# analysis/scripts/parse_log.py
from analysis import load_single_log
df = load_single_log(snakemake.input.log)
df.to_parquet(snakemake.output.parquet)