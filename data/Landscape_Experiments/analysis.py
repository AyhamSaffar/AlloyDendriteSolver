# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pathlib as pl

# %%
experiment_path = pl.Path(__file__).parent / input('Enter experiment name: ')
assert experiment_path.is_dir(), 'Given experiment name not found.'

# %%
raw_data = pd.read_csv(experiment_path / 'scan_data.csv')
Vs = np.sort(raw_data['V'].unique())
Rs = np.sort(raw_data['R'].unique())

data = {} # dT: {col: col_grid}
for dT in np.sort(raw_data['dT'].unique()):
	subset = raw_data[raw_data['dT']==dT].sort_values(by=['V', 'R'], ascending=[False, True]) # top row = highest V0
	grids = {col: subset[col].to_numpy().reshape([len(Vs), len(Rs)]) for col in ['f1', 'f2']} 
	data[dT] = grids

# %%
fig, axes = plt.subplots(ncols=len(data), figsize=(7*len(data), 5))

for i, (dT, grids) in enumerate(data.items()):

	axes[i].set_title(f'Undercooling = {dT:.0f}K')
	axes[i].set_xlabel(r'$R_0$ / $log_{10}(m)$')
	axes[i].set_ylabel(r'$V_0$ / $log_{10}(m/s)$')

	log_f1, log_f2 = np.log10(np.abs(grids['f1'])), np.log10(np.abs(grids['f2']))
	F_mean_im = axes[i].imshow(
		(log_f1 + log_f2) / 2,
		extent=(np.log10(Rs.min()), np.log10(Rs.max()), np.log10(Vs.min()), np.log10(Vs.max())),
		interpolation=None,
		aspect='auto',
		cmap='coolwarm',
	)

	fig.colorbar(F_mean_im, label=r'Mean|F| / $log_{10}$')

fig.savefig(experiment_path / "plots.png")

