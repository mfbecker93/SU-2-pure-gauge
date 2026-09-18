from analysis import plot_plaquette_history, savefig
import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv(snakemake.input.csv)
plot_plaquette_history(df, title="Plaquette History", ax=None)  
savefig(plt.gcf(), snakemake.output.png)

