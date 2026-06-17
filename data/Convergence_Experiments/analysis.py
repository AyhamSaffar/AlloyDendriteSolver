# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pathlib as pl

# %%
home_path = pl.Path(__file__).parent
experiment_path = home_path / 'NiB_Global_Newton'
assert experiment_path.is_dir()

# %%
raw_data = pd.read_csv(experiment_path / 'solver_data.csv')
V0s = np.sort(raw_data['V0'].unique())
R0s = np.sort(raw_data['R0'].unique())

data = {} # dT: {col: col_grid}
for dT in np.sort(raw_data['dT'].unique()):
	subset = raw_data[raw_data['dT']==dT].sort_values(by=['V0', 'R0'], ascending=[False, True]) # top row = highest V0
	grids = {col: subset[col].to_numpy().reshape([len(V0s), len(R0s)]) for col in raw_data.columns} 
	grids['valid'] = (grids['converged']) & (grids['R']>0) & (grids['V']>0)
	grids['f'] = np.abs(grids['f1']) + np.abs(grids['f2'])
	data[dT] = grids

for dT, subset in data.items():
	steps = np.mean(subset['steps'][subset['converged']])
	print(f'dT={dT}, average steps={steps:.1f}')

approx_data = pd.read_csv(experiment_path / 'approx_data.csv')


# %%
fig, axes = plt.subplots(nrows=2, ncols=len(data), figsize=(7*len(data),10))

for i, (dT, grids) in enumerate(data.items()):
	approx_row = approx_data[approx_data['dT']==dT]

	axes[0, i].set_title(f'Undercooling = {dT:.0f}K')
	for row in [0, 1]:
		axes[row, i].set_xlabel(r'$R_0$ / $log_{10}(m)$')
		axes[row, i].set_ylabel(r'$V_0$ / $log_{10}(m/s)$')
		axes[row, i].scatter(np.log10(approx_row['R']), np.log10(approx_row['V']), color='black', marker='X')

	R_im = axes[0, i].imshow(
		np.log10(grids['R'],
		out=np.full(shape=grids['R'].shape, fill_value=np.nan), where=grids['valid']),
		extent=(np.log10(R0s.min()), np.log10(R0s.max()), np.log10(V0s.min()), np.log10(V0s.max())),
		aspect='auto',
	)
	fig.colorbar(R_im, label=r'$R$ / $log_{10}(m)$')

	neg_R_grid = np.full(shape=grids['R'].shape, fill_value=np.nan)
	neg_R_grid[(grids['converged']) & (grids['R']<0)] = 1.0
	neg_R_grid[0, 0] = 0.0

	neg_R_im = axes[0, i].imshow(
		neg_R_grid,
		extent=(np.log10(R0s.min()), np.log10(R0s.max()), np.log10(V0s.min()), np.log10(V0s.max())),
		cmap='hsv',
		aspect='auto',
	)

	V_im = axes[1, i].imshow(
		np.log10(grids['V'],
		out=np.full(shape=grids['V'].shape, fill_value=np.nan), where=grids['valid']),
		extent=(np.log10(R0s.min()), np.log10(R0s.max()), np.log10(V0s.min()), np.log10(V0s.max())),
		aspect='auto',
	)
	fig.colorbar(V_im, label=r'$V$ / $log_{10}(m/s)$')

	neg_V_grid = np.full(shape=grids['V'].shape, fill_value=np.nan)
	neg_V_grid[(grids['converged']) & (grids['V']<0)] = 1.0
	neg_V_grid[0, 0] = 0.0

	neg_V_im = axes[0, i].imshow(
		neg_V_grid,
		extent=(np.log10(R0s.min()), np.log10(R0s.max()), np.log10(V0s.min()), np.log10(V0s.max())),
		cmap='hsv',
		aspect='auto',
	)

	# steps_im = axes[1, i].imshow(
	# 	np.log10(grids['steps']+1,
	# 	out=np.full(shape=grids['R'].shape, fill_value=np.nan), where=(grids['diverged'] | ~grids['converged'])),
	# 	extent=(np.log10(R0s.min()), np.log10(R0s.max()), np.log10(V0s.min()), np.log10(V0s.max())),
	# 	cmap='Spectral',
	# 	aspect='auto',
	# )
	# fig.colorbar(steps_im, label=r'$log_{10}$ steps')


	fig.savefig(experiment_path / "plots.png")

