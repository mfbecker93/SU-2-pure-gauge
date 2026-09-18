from analysis import thermalized
import pandas as pd
df = pd.read_csv(snakemake.input.csv)
df_therm = thermalized(df, trajectory_cut=snakemake.config["trajectory_cut"])
df_therm.to_csv(snakemake.output.csv, index=False)