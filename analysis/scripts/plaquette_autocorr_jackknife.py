from analysis import jackknife_autocorr_scan
import pandas as pd
df = pd.read_csv(snakemake.input.csv)
autocorr_df = jackknife_autocorr_scan(df, value='plaquette', groupby='beta', max_lag=snakemake.config["Wmax"], n_blocks=snakemake.config["n_blocks"])
autocorr_df.to_csv(snakemake.output.csv, index=False)