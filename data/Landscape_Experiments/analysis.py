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
fig, axes = plt.subplots(ncols=len(data), nrows=3, figsize=(5*len(data), 12))

for col, (dT, grids) in enumerate(data.items()):
	log_f1, log_f2 = np.log10(np.abs(grids['f1'])), np.log10(np.abs(grids['f2']))
	axes[0, col].set_title(f'Undercooling = {dT:.0f}K', fontsize='x-large')

	for row in range(3):
		axes[row, col].set_xlabel(r'$R_0$ / $log_{10}(m)$')
		axes[row, col].set_ylabel(r'$V_0$ / $log_{10}(m/s)$')

		if col==0:
			label = f'$log_{{10}}$ |{'f1+f2' if row==0 else 'f1' if row==1 else 'f2'}|'
			axes[row, col].annotate(
				label, xy=(-0.3, 0.5), xycoords='axes fraction', rotation='vertical',
				fontsize='x-large', verticalalignment='center'
			)

		im = axes[row, col].imshow(
			log_f1+log_f2 if row==0 else log_f1 if row==1 else log_f2,
			extent=(np.log10(Rs.min()), np.log10(Rs.max()), np.log10(Vs.min()), np.log10(Vs.max())),
			interpolation=None,
			aspect='auto',
			cmap='coolwarm',
		)
		fig.colorbar(im)

fig.tight_layout(rect=[0.01, 0, 1, 1]) # ensures row labels are not cut off on left edge
fig.savefig(experiment_path / "plots.png")
