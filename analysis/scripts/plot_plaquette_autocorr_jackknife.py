from analysis import plot_autocorr_scan, savefig
import pandas as pd 
import matplotlib.pyplot as plt

df=pd.read_csv(snakemake.input.csv)

plot_autocorr_scan(
    df,
    betas=snakemake.config["betas"],
    max_lag=20, #TODO: add to config.yaml
    title="Autocorrelation function of plaquette for different β values" #TODO: add to config.yaml
)
savefig(plt.gcf(), snakemake.output.png)
