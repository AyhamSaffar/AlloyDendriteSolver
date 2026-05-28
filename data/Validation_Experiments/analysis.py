# %%
import matplotlib.pyplot as plt
import pandas as pd
import pathlib as pl
import numpy as np

# %%
home_path = pl.Path(__file__).parent 
experiment_path = home_path / '27_May_test'
experiments = {path.stem: pd.read_csv(path) for path in experiment_path.glob(pattern='*.csv')}

# %%
data = experiments['SucAce_LGK']
fig, axes = plt.subplots(nrows=2, figsize=(5,7))
suc_Mr, ace_Mr = 80.09, 58.08

axes[0].set_ylabel(r'$ V (10^{-5} ms^{-1}) $')
axes[0].set_ylim(0, 22)
axes[0].set_yticks(range(0, 21, 4))
axes[0].set_xlabel(r'$ C_{Ace} (mol.\%) $')
axes[0].set_xlim(0, 1)

axes[1].set_ylabel('$ R(10^{-6} m) $')
axes[1].set_ylim(0, 30)
axes[1].set_yticks(range(0, 29, 4))
axes[1].set_xlabel(r'$ C_{Ace} (mol.\%) $')
axes[1].set_xlim(0, 1)

for dT in data['dT'].unique():
    subset = data[data['dT']==dT]
    C0_mol_percent = 100 * (subset['C0']/ace_Mr) / (subset['C0']/ace_Mr + (100-subset['C0'])/suc_Mr)
    axes[0].plot(C0_mol_percent, subset['V']*1e5, label=f'$ ΔT = {dT} K $')
    axes[1].plot(C0_mol_percent, subset['R']*1e6, label=f'$ ΔT = {dT} K $')

axes[0].legend()
axes[1].legend()

fig.savefig(home_path / 'SucAce_LGK.png')

# %%
data = experiments['AlFe_LGK']
fig, ax = plt.subplots(figsize=(6, 6))

ax.set_xlabel('Undercooling (K)')
ax.set_xscale('log')
ax.set_xlim(1, 10**2.5)
ax.set_xticks([1, 10, 100])
ax.set_ylabel('Growth Velocity (cm/s)')
ax.set_yscale('log')
ax.set_ylim(1e-4, 1e3)

for C0 in data['C0'].unique():
    subset = data[data['C0']==C0]
    ax.plot(subset['dT'], subset['V']*100, label=f'Al-{C0:.1f}wt% Fe')

ax.legend()
fig.savefig(home_path / 'AlFe_LGK.png')

# %%
data = experiments['NiSn_LGK']
fig, axes = plt.subplots(ncols=2, figsize=(12,8))

axes[0].set_xlabel('BULK UNDERCOOLING K')
axes[0].set_xscale('log')
axes[0].set_xlim(1, 10**3.5)
axes[0].set_ylabel('DENDRITE TIP VELOCITY, m/s')
axes[0].set_yscale('log')
axes[0].set_ylim(1e-4, 1e3)

axes[1].set_xlabel('BULK UNDERCOOLING, K')
axes[1].set_xlim(0, 350)
axes[1].set_xticks(range(0, 301, 100))
axes[1].set_ylabel('DENDRITE TIP RADIUS, m')
axes[1].set_yscale('log')
axes[1].set_ylim(1e-10, 1e-2)

axes[0].plot(data['dT'], data['V'])
axes[1].plot(data['dT'], data['R'])

fig.savefig(home_path / 'NiSn_LGK.png')
