# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pathlib as pl

# %%
home_path = pl.Path(__file__).parent
experiment_paths = [path for path in home_path.glob(pattern='*') if path.is_dir()]

# %%
def read_data(path: pl.Path):
	data = pd.read_csv(path / 'solver_data.csv')
	V0s = np.sort(data['V0'].unique())
	dTs = np.sort(data['dT'].unique())

	grids = {col: data[col].to_numpy().reshape([len(dTs), len(V0s)]).T[::-1] for col in data.columns} # top row should be highest V0 value
	grids['f'] = np.abs(grids['f1']) + np.abs(grids['f2'])
	grids['f_log'] = np.full(shape=[len(V0s), len(dTs)], fill_value=np.nan)
	np.log10(grids['f']+1e-20, out=grids['f_log'], where=~grids['diverged'])
	
	approx_data = pd.read_csv(path / 'approx_data.csv')
	approx_v = approx_data['V'].to_numpy()[np.newaxis, :]
	percent_error = 100 * np.abs(grids['V'] -  approx_v) / approx_v
	percent_error[grids['diverged']] = np.nan

	return approx_data, V0s, dTs, grids, percent_error

# %%
fig, axes = plt.subplots(nrows=3, ncols=len(experiment_paths), figsize=(28,10))
fig.suptitle('Sn-3.5Ag Non Diverged Data Points')

for i, path in enumerate(experiment_paths):
	axes[0, i].set_title(path.stem)
	approx_data, V0s, dTs, grids, percent_error = read_data(path)

	for row in (0, 1, 2):
		axes[row, i].plot(approx_data['dT'], np.log10(approx_data['V']), color='blue')
		axes[row, i].set_xlabel(r'$\Delta T$ / $K$')
		axes[row, i].set_ylabel(r'$log_{10}V0$ / $ms^{-1}$')

	f_im = axes[0, i].imshow(
		grids['f_log'], vmax=0, vmin=-12, extent=(dTs.min(), dTs.max(), np.log10(V0s.min()), np.log10(V0s.max())),
		aspect='auto', alpha=0.7, cmap='RdYlGn_r'
	)
	fig.colorbar(f_im, label=r'$log_{10}$Convergence')

	error_im = axes[1, i].imshow(
		np.log10(percent_error),
		extent=(dTs.min(), dTs.max(), np.log10(V0s.min()), np.log10(V0s.max())),
		vmin=0, aspect='auto', alpha=0.7, cmap='RdYlGn_r',
		# vmax=3
	)
	fig.colorbar(error_im, label=r'$log_{10}$ V % Error Compared to V approx')

	R_im = axes[2, i].imshow(
		np.log10(grids['R'], out=np.full(shape=grids['R'].shape, fill_value=np.nan), where=(grids['converged'])&(grids['R']>0)),
		extent=(dTs.min(), dTs.max(), np.log10(V0s.min()), np.log10(V0s.max())),
		aspect='auto',
	)
	fig.colorbar(R_im, label=r'$log_{10}$ R / m')

	neg_R_grid = np.full(shape=grids['R'].shape, fill_value=np.nan)
	neg_R_grid[(grids['converged']) & (grids['R']<0)] = 1.0
	neg_R_grid[0, 0] = 0.0

	neg_R_im = axes[2, i].imshow(
		neg_R_grid, extent=(dTs.min(), dTs.max(), np.log10(V0s.min()), np.log10(V0s.max())), aspect='auto', cmap='hsv',
	)

	axes[0, 0].annotate('non diverged', xy=(-0.4, 0.5), xycoords='axes fraction')
	axes[1, 0].annotate('non diverged', xy=(-0.4, 0.5), xycoords='axes fraction')
	axes[2, 0].annotate('converged', xy=(-0.4, 0.5), xycoords='axes fraction')
	axes[2, 0].annotate('red = R < 0', xy=(-0.4, 0.3), xycoords='axes fraction')

	fig.savefig(home_path / "plots.png")

