# %%
import pandas as pd
import numpy as np
import pathlib as pl
import matplotlib.pyplot as plt

# %%
home_path = pl.Path(__file__).parent / input("Enter phase diagram folder name: ")
assert home_path.is_dir(), f"Could not find folder at {home_path}"
data_path = home_path / "data.csv"
assert data_path.is_file(), f"Could not find data.csv file in {home_path}"

# %%
data = pd.read_csv(data_path)
for col in ["T / K", "Cl", "Cs"]:
    assert col in data.columns, f"Could not find {col} Column in {data.columns.to_list()} columns"

# %%
orders = range(2, 8)
fig, axes = plt.subplots(ncols=len(orders), nrows=2, figsize=(6*len(orders),10))
data['k0'] = data['Cs'] / data['Cl']
Tl_fits, all_Tl_fit_errors, k0_fits, all_k0_fit_errors = [], [], [], []

for i, order in enumerate(orders):
    axes[0, i].set_title(f'{order} order Fits')

    # plotting Cl to Tl fit
    axes[0, i].set_xlabel("$ C_L $")
    axes[0, i].set_ylabel("T / K")
    axes[0, i].scatter(data['Cl'], data['T / K'], color='black', marker='x')
    Tl_fit = np.polynomial.Polynomial.fit(data['Cl'], data['T / K'], deg=order)
    axes[0, i].plot(np.arange(1, 101), Tl_fit(np.arange(1, 101)), color='red', linestyle='--')
    twin_ax = axes[0, i].twinx()
    twin_ax.set_ylabel('error')
    T_errors = np.abs(data['T / K']-Tl_fit(data['Cl']))
    twin_ax.bar(data['Cl'], T_errors, color='grey', alpha=0.2)
    Tl_fits.append(Tl_fit)
    all_Tl_fit_errors.append(T_errors)

    # plotting T to k0 fit
    axes[1, i].set_xlabel("T / K")
    axes[1, i].set_ylabel("k0")
    axes[1, i].set_ylim(0, 1)
    axes[1, i].scatter(data['T / K'], data['k0'], color='black', marker='x')
    k0_fit = np.polynomial.Polynomial.fit(data['T / K'], data['k0'], deg=order)
    fit_T, fit_k0 = k0_fit.linspace()
    axes[1, i].plot(fit_T, fit_k0,  color='red', linestyle='--')
    twin_ax = axes[1, i].twinx()
    twin_ax.set_ylabel('error')
    k0_errors = np.abs(data['k0']-k0_fit(data['T / K']))
    twin_ax.bar(data['T / K'], k0_errors, color='grey', alpha=0.5)
    k0_fits.append(k0_fit)
    all_k0_fit_errors.append(k0_errors)

fig.tight_layout()
fig.savefig(home_path / "fits.png")

# %%
all_fit_stats = []
for fits, error_arrays in ((Tl_fits, all_Tl_fit_errors), (k0_fits, all_k0_fit_errors)):
    for fit, error_array in zip(fits, error_arrays):
        fit_stats = {"kind": "Cl to Tl" if fits==Tl_fits else "T to k0"}
        coefs = fit.convert().coef
        fit_stats['order'] = len(coefs) - 1 # -1 as 0th order coefficient exists 
        for order in range(max(orders)+1):
            fit_stats[f'{order} order coefficient'] = coefs[order] if order<len(coefs) else 0
        fit_stats['mean error'] = np.mean(error_array)
        fit_stats['max error'] = np.max(error_array)
        all_fit_stats.append(fit_stats)

pd.DataFrame(all_fit_stats).to_csv(home_path / "fits.csv", index=False)

