from analysis import select_scan
import pandas as pd
df = pd.read_csv(snakemake.input.log)
scan_df = select_scan(
    df,
    beta=snakemake.config["betas"],
    mdsteps=snakemake.config["mdsteps"],
    trajL=snakemake.config["trajL"],
    lattice=snakemake.config["lattice"],
    trajectories=snakemake.config["trajectories"],
                      ) 
scan_df.to_csv(snakemake.output.csv, index=False)  