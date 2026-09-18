from analysis import integrated_autocorr_time_scan
import pandas as pd
df = pd.read_csv(snakemake.input.csv)
autocorr_df = integrated_autocorr_time_scan(df, value='plaquette', groupby='beta', Wmax=snakemake.config["Wmax"])
autocorr_df.to_csv(snakemake.output.csv, index=False)