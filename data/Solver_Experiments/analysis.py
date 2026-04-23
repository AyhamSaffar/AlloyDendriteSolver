# %%
import matplotlib.pyplot as plt
import pandas as pd
import pathlib as pl
import numpy as np

# %%
# must have more than 1 experiment
home_path = pl.Path(__file__).parent
experiments = {path.name: pd.read_csv(path/"data.csv") for path in home_path.glob(pattern='Fig11d*') if path.is_dir()}

# %%
fig, ax = plt.subplots(figsize=(8,6))
ax.set_yscale("log")
ax.set_ylim(1e-5, 1e1)
ax.set_xlim(0, 50)

for i, (name, data) in enumerate(experiments.items()):
	model = name.lstrip('Fig11d_')
	for C0 in data['C0'].unique():
		subset = data[data['C0']==C0]
		converged = subset[subset['converged']]
		ax.plot(converged['dT'], converged['V'], label=f'Sn-{C0:.1f}Ag {model}')
		ax.vlines(
			subset['V'][subset['diverged']], ymin=converged['V'].min(), ymax=converged['V'].max(),
			colors='red', alpha=0.1
		)
ax.set_xlabel('dT / K')
ax.set_ylabel('V / $ms^{-1}$')
ax.legend()

fig.savefig(home_path / 'plots.png')
