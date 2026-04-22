# %%
import matplotlib.pyplot as plt
import pandas as pd
import pathlib as pl
import numpy as np

# %%
# must have more than 1 experiment
home_path = pl.Path(__file__).parent
experiments = {path.name: pd.read_csv(path/"data.csv") for path in home_path.glob(pattern='AgSn*') if path.is_dir()}

# %%
fig, axes = plt.subplots(nrows=2, ncols=len(experiments), figsize=(16,10))
for i, (name, data) in enumerate(experiments.items()):
	axes[0, i].set_title(name)
	axes[1, i].set_title("Convergence")
	for C0 in data['C0'].unique():
		subset = data[data['C0']==C0]
		converged = subset[subset['converged']]
		# V plots
		axes[0, i].scatter(converged['V'], converged['dT'], label=f'{C0}wt% Ag')
		axes[0, i].set_xlim(0, 1.5)
		axes[0, i].vlines(subset['V'][subset['diverged']], ymin=0, ymax=60, colors='red', alpha=0.1)
		axes[0, i].legend()
		axes[0, i].set_xlabel('V / $ms^{-1}$')
		axes[0, i].set_ylabel('dT / K')

		# Convergence plots
		axes[1, i].scatter(subset['V'], subset['f1'].abs()+subset['f2'].abs(), label=f'{C0}wt% Ag f')
		axes[1, i].set_xlim(0, 1.5)
		all_fs = np.hstack([data['f1'].abs(), data['f2'].abs()])
		axes[1, i].vlines(subset['V'][subset['diverged']], ymin=all_fs.min(), ymax=all_fs.max(), colors='red', alpha=0.1)
		axes[1, i].set_yscale('log') # note ignores all f values that equal 0 (ie when solver diverges)
		axes[1, i].legend()
		axes[1, i].set_xlabel('V / $ms^{-1}$')
		axes[1, i].set_ylabel('abs(f)')

fig.savefig(home_path / 'plots.png')
