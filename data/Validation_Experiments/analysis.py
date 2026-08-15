# %%
import matplotlib.pyplot as plt
import pandas as pd
import pathlib as pl
import numpy as np

# %%
experiment_path = pl.Path(__file__).parent / input('Enter experiment name: ')
assert experiment_path.is_dir(), 'Given experiment name not found.'
experiments = {path.stem: pd.read_csv(path) for path in experiment_path.glob(pattern='*.csv')}

# %%
data = experiments['SucAce_LGK']
fig, axes = plt.subplots(nrows=2, figsize=(5,7))

bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
for ax in axes.flatten():
    ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

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

suc_Mr, ace_Mr = 80.09, 58.08
for dT in data['dT'].unique():
    subset = data[data['dT']==dT]
    C0_mol_percent = 100 * (subset['C0']/ace_Mr) / (subset['C0']/ace_Mr + (100-subset['C0'])/suc_Mr)
    axes[0].plot(C0_mol_percent, subset['V']*1e5, label=f'$ ΔT = {dT} K $')
    axes[1].plot(C0_mol_percent, subset['R']*1e6, label=f'$ ΔT = {dT} K $')

axes[0].legend()
axes[1].legend()

fig.tight_layout()
fig.savefig(experiment_path / 'SucAce_LGK.png')

# %%
data = experiments['AlFe_LGK']
fig, ax = plt.subplots(figsize=(6, 6))

bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

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
fig.tight_layout()
fig.savefig(experiment_path / 'AlFe_LGK.png')

# %%
data = experiments['NiSn_LGK']
fig, axes = plt.subplots(ncols=2, figsize=(12,8))

bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
for ax in axes.flatten():
    ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

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

fig.tight_layout()
fig.savefig(experiment_path / 'NiSn_LGK.png')

# %%
data = experiments['AgCu_LKT_BCT']
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
for ax in axes.flatten():
    ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

for ax in axes.flatten():
    ax.set_xlabel('Undercooling K')
    ax.set_xscale('log')
    ax.set_xlim(1, 600)

axes[0,0].set_ylabel('V     cm/s')
axes[0,0].set_yscale('log')
axes[0,0].set_ylim(1e-5, 1e6)
axes[0,0].scatter(data['dT'], data['V']*100, color='black')

axes[0,1].set_ylabel('Distribution Coefficient')
axes[0,1].set_ylim(0.3, 1.01)
axes[0,1].scatter(data['dT'], data['k'], color='black')

axes[1,0].set_ylabel('Solid Concentration   wt% Cu')
axes[1,0].set_ylim(4, 18)
axes[1,0].scatter(data['dT'], data['Cs'], color='black')

axes[1,1].set_ylabel('Radius    cm')
axes[1,1].set_yscale('log')
axes[1,1].set_ylim(1e-6, 10**(-2.5))
axes[1,1].scatter(data['dT'], data['R']*100, color='black')

fig.tight_layout()
fig.savefig(experiment_path / 'AgCu_LKT_BCT.png')

# %%
data_gamma = experiments['FeCoGamma_LKT_BCT']
data_delta = experiments['FeCoDelta_LKT_BCT']
fig, axes = plt.subplots(ncols=3, figsize=(16, 4))

data = pd.concat([data_gamma, data_delta])
bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
for ax in axes.flatten():
    ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

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

fig.tight_layout()
fig.savefig(experiment_path / 'FeCo_LKT_BCT.png')

# %%
data_LGK, data_LKT_BCT = experiments['SnAg_LGK'], experiments['SnAg_LKT_BCT']
fig, ax = plt.subplots(figsize=(8,6))

data = pd.concat([data_LGK, data_LKT_BCT])
bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

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

fig.tight_layout()
fig.savefig(experiment_path / 'SnAg_LKT_BCT.png')

# %%
data = experiments['NiB_LKT_BCT']
fig, axes = plt.subplots(nrows=3, figsize=(7,13))

bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
for ax in axes.flatten():
    ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

axes[0].set_xlabel('Undercooling    ΔT (K)')
axes[0].set_xlim(0, 320)
axes[0].set_xticks(range(0, 301, 100))
axes[0].set_ylabel('Dendrite growth velocity V (m/s)')
axes[0].set_ylim(0, 30)
axes[0].set_yticks(range(0, 31, 10))

axes[0].axvline(214, color='orange', linestyle=':')
axes[0].axvline(267, color='green', linestyle=':')

for ax in [axes[1], axes[2]]:
    ax.set_xlabel('Undercooling ΔT (K)')
    ax.set_xlim(0, 400)
    ax.set_xticks(range(0, 401, 100))
    ax.axvline(x=267, linestyle=':')

axes[1].set_ylabel('Concentration (at% B)')
axes[1].set_ylim(0, 7)
axes[1].set_yticks(range(8))
axes[1].axhline(y=1, label='C0', linestyle=':')

axes[2].set_ylabel('Dendrite tip radius R (m)')
axes[2].set_yscale('log')
axes[2].set_ylim(1e-8, 1e-5)
axes[2].set_yticks([1e-5, 1e-6, 1e-7, 1e-8])

for C0 in [0, 0.7, 1]:
    subset = data[data['C0']==C0]
    axes[0].plot(subset['dT'], subset['V'], label=f'x = {C0:.1f}')
    if C0==0:
        continue
    axes[2].plot(subset['dT'], subset['R'], label=f'$ Ni_{{{100-C0}}}B_{{{C0}}} $')
    if C0==0.7:
        continue
    axes[1].plot(subset['dT'], subset['Cs'], label='$ C_s^* $')
    axes[1].plot(subset['dT'], subset['Cl'], label='$ C_l^* $')
    m = -14.3 # NiB liquidus slope
    axes[1].plot(subset['dT'], subset['C0']-subset['dT']/m, linestyle=':', label='$ C_l^{eq.} $')

axes[0].legend(title='$Ni_{100-x}B_x$', loc='upper left')
axes[1].legend(title='$Ni_{99}B_1$', loc='upper left')
axes[2].legend()

fig.tight_layout()
fig.savefig(experiment_path / 'NiB_LKT_BCT.png')

# %%
data = experiments['CoCu_CLW']
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(12,8))

bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
for ax in axes.flatten():
    ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

# velocity plots
for col in [0, 1]:
    ax = axes[0, col]
    ax.set_xlabel('ΔT (K)')
    ax.set_ylabel('V (m/s)')
    subset = data[data['C0'] == (20 if col==0 else 60)]
    ax.plot(subset['dT'], subset['V'], color='black')
    

axes[0, 0].set_xlim(0, 300)
axes[0, 0].set_xticks(range(0, 301, 100))
axes[0, 0].set_ylim(0, 40)
axes[0, 0].set_yticks(range(0, 41, 10))

axes[0, 1].set_xlim(50, 120)
axes[0, 1].set_xticks(range(60, 121, 20))
axes[0, 1].set_ylim(-0.025, 0.2)
axes[0, 1].set_yticks(np.arange(0, 0.21, 0.05))


# undercooling plots
for col in [0, 1]:
    ax = axes[1, col]
    ax.set_xlabel('Bulk undercooling (K)')
    ax.set_ylabel('Partial undercoolings (K)')
    ax.set_yscale('log')
    subset = data[data['C0'] == (20 if col==0 else 60)]
    ax.plot(subset['dT'], subset['dTt'], label='$ ΔT_t $')
    ax.plot(subset['dT'], subset['dTc'], label='$ ΔT_c $')
    ax.plot(subset['dT'], subset['dTr'], label='$ ΔT_r $')
    ax.plot(subset['dT'], subset['dTk'], label='$ ΔT_k $')
    ax.legend()


axes[1, 0].set_xlim(0, 350)
axes[1, 0].set_xticks(range(0, 301, 100))
axes[1, 0].set_ylim(0.1, 250)
axes[1, 0].set_yticks([1, 10, 100])

axes[1, 1].set_xlim(0, 120)
axes[1, 1].set_xticks(range(0, 121, 20))
axes[1, 1].set_ylim(0.01, 120)
axes[1, 1].set_yticks([0.1, 1, 10, 100])


fig.tight_layout()
fig.savefig(experiment_path / 'CoCu_CLW.png')

# %%
data = experiments['FeSb_CLW']
fig, axes = plt.subplots(nrows=3, sharex=True, figsize=(6, 7))

bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
for ax in axes.flatten():
    ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

axes[2].set_xlabel("ΔT (K)")
axes[2].set_xlim(0, 450)
axes[2].set_xticks(range(0, 401, 100))

axes[1].set_ylabel("V (m/s)")
axes[1].set_ylim(0, 2)
axes[1].set_yticks(range(3))

axes[0].set_ylabel("V (m/s)")
axes[0].set_yscale("log")
axes[0].set_ylim(2, 1e3)
axes[0].set_yticks([2, 10, 100, 1000])

for i in range(2):
    axes[i].plot(data['dT'], data['V'], color='black', linestyle='-.')

axes[2].tick_params(labelleft=False)
ax = axes[2].twinx()
ax.set_ylabel("k")
ax.set_ylim(0, 1)
ax.set_yticks(np.arange(0, 1.01, 0.25))

ax.plot(data['dT'], data['k0'], color='black', linestyle='-.')
ax.plot(data['dT'], data['kv'], color='black', linestyle='-')

fig.tight_layout()
fig.savefig(experiment_path / 'FeSb_CLW.png')

# %%
data = experiments['NiB_WLCYZ']
fig, axes = plt.subplots(ncols=3, figsize=(12, 3.5))

bad_rows = (~data['converged']) | (data['V']<0) | (data['R']<0)
for ax in axes.flatten():
    ax.vlines(data.loc[bad_rows, 'dT'], ymin=0, ymax=1e10, color='red', alpha=0.2)

axes[0].set_xlabel('bath undercooling ΔT (K)')
axes[0].set_xlim(0, 350)
axes[0].set_xticks(range(0, 350, 100))
axes[0].set_ylabel('The dendrite tip velocity V (m/s)')
axes[0].set_ylim(-1, 40)
axes[0].set_yticks(range(0, 41, 10))

axes[0].plot(data['dT'], data['V'], color='black')

axes[1].set_xlabel('bath undercooling ΔT (K)')
axes[1].set_xlim(0, 400)
axes[1].set_xticks([0, 100, 300])
axes[1].set_ylabel('The dendrite tip radius R (m)')
axes[1].set_yscale('log')
axes[1].set_ylim(1e-9, 1e-6)
axes[1].set_yticks([10**i for i in range(-9, -5)])

axes[1].plot(data['dT'], data['R'], color='black')

axes[2].set_xlabel('bath undercooling ΔT (K)')
axes[2].set_xlim(0, 400)
axes[2].set_xticks([0, 100, 300])
axes[2].set_ylabel('undercooling contribution (K)')
axes[2].set_ylim(-5, 225)
axes[2].set_yticks(range(0, 210, 50))

axes[2].plot(data['dT'], data['dTr'], color='black', linestyle='-', label='dTr')
axes[2].plot(data['dT'], data['dTc'], color='black', linestyle='--', label='dTc')
axes[2].plot(data['dT'], data['dTk'], color='black', linestyle=':', label='dTk')
axes[2].plot(data['dT'], data['dTt'], color='black', linestyle='-.', label='dTt')

axes[2].legend()

fig.tight_layout()
fig.savefig(experiment_path / 'NiB_WLCYZ.png')
