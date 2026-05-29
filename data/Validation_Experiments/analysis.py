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

# %%
data_gamma = experiments['FeCoGamma_LKT_BCT']
data_delta = experiments['FeCoDelta_LKT_BCT']
fig, axes = plt.subplots(ncols=3, figsize=(16, 4))

for ax in axes:
    ax.set_xlabel('ΔT/K')
    ax.set_xlim(0, 350)
    ax.set_xticks(range(0, 301, 100))
    ax.set_ylabel('$ V/m-s^{-1} $')
    ax.set_ylim(0, 30)
    ax.set_yticks(range(0, 31, 5))

for i, C0 in enumerate(data_gamma['C0'].unique()):
    axes[i].set_title(f'Fe-{C0:.0f} at.% Co')
    subset_gamma = data_gamma[data_gamma['C0']==C0]
    axes[i].plot(subset_gamma['dT'], subset_gamma['V'], color='black', linestyle='--')
    subset_delta = data_delta[data_delta['C0']==C0]
    T_offset = 10 if C0==30 else 24 if C0==40 else 40 # delta phase nucleates after gamma as it has a lower solidus T
    axes[i].vlines(T_offset, ymin=0, ymax=30)
    axes[i].plot(subset_delta['dT']+T_offset, subset_delta['V'], color='black', linestyle='-')


fig.savefig(home_path / 'FeCo_LKT_BCT.png')

# %%
data_LGK, data_LKT_BCT = experiments['SnAg_LGK'], experiments['SnAg_LKT_BCT']
fig, ax = plt.subplots(figsize=(8,6))

ax.set_ylabel('Growth velocity (m/s)')
ax.set_yscale("log")
ax.set_ylim(1e-5, 1e1)
ax.set_yticks([10**power for power in range(-5, 2)], [f'1E{power:+.0f}' for power in range(-5, 2)])
ax.set_xlabel('Growth undercooling (K)')
ax.set_xlim(0, 50)
ax.set_xticks(range(0, 51, 10))

for model, data in [('LGK', data_LGK), ('LKT_BCT', data_LKT_BCT)]:
    linestyle = '-' if model=='LGK' else '--'
    for color, C0 in [('gray', 3.5), ('red', 5.0)]:
        subset = data[data['C0']==C0]
        ax.plot(subset['dT'], subset['V'], label=f'Sn-{C0:.1f}Ag {model}', color=color, linestyle=linestyle)

ax.legend()

fig.savefig(home_path / 'SnAg_LKT_BCT.png')
