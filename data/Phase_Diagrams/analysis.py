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
for col in ["T", "Cl", "Cs"]:
    assert col in data.columns, f"Could not find {col} Column in {data.columns.to_list()} columns"

# %%
orders = range(3, 8)
fits = [('Cl', 'T'), ('T', 'Cl'), ('T', 'Cs')] # parameter pairs to fit (x axis, y axis)
fig, axes = plt.subplots(ncols=len(orders), nrows=len(fits), figsize=(6*len(orders), 5*len(fits)))
stat_data = []

for i, order in enumerate(orders):
    axes[0, i].set_title(f'{order} order Fits')

    for j, (x, y) in enumerate(fits):
        ax = axes[j, i]

        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.scatter(data[x], data[y], color='black', marker='x')

        fit = np.polynomial.Polynomial.fit(data[x], data[y], deg=order)
        y_fit = fit(data[x])
        ax.plot(data[x], y_fit, color='red', linestyle='--')
        ax = ax.twinx()
        ax.set_ylabel(f'Error')
        y_errors = np.abs(data[y]-y_fit)
        ax.bar(data[x], y_errors, color='grey', alpha=0.5)

        stat_row = dict()
        stat_row['kind'] = f'{x} to {y}'
        stat_row['order'] = order
        coefs = fit.convert().coef
        for coef in range(max(orders)+1):
            stat_row[f'{coef} order coefficient'] = coefs[coef] if coef<len(coefs) else 0
        stat_row['mean error'] = np.mean(y_errors)
        stat_row['max error'] = np.max(y_errors)
        stat_data.append(stat_row)


fig.tight_layout()
fig.savefig(home_path / "fits.png")
df = pd.DataFrame(stat_data).sort_values(by=['kind', 'order'])
df.to_csv(home_path / "fits.csv", index=False)

