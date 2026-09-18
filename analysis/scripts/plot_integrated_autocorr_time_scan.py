from analysis import plot_int_autocorr_scan, savefig
import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv(snakemake.input.csv)
plot_int_autocorr_scan(df, betas=snakemake.config["betas"], Wmax=snakemake.config["Wmax"], ax=None)
savefig(plt.gcf(), snakemake.output.png)