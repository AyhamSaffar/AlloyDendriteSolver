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
for col in ["T/K", "Cl/%", "Cs/%"]:
    assert col in data.columns, f"Could not find {col} Column in {data.columns.to_list()} columns"

Cl_splits = input(r'Split up liquidus fits? Enter temperatures / K seperated by a space or press enter: ')
Cl_splits = [data['T/K'].min()] + sorted([float(i) for i in Cl_splits.split(' ') if i!='']) + [data['T/K'].max()]
Cl_ranges = [[Cl_splits[i], Cl_splits[i+1]] for i in range(len(Cl_splits)-1)]

Cs_splits = input(r'Split up solidus fits? Enter temperatures / K seperated by a space or press enter: ')
Cs_splits = [data['T/K'].min()] + sorted([float(i) for i in Cs_splits.split(' ') if i!='']) + [data['T/K'].max()]
Cs_ranges = [[Cs_splits[i], Cs_splits[i+1]] for i in range(len(Cs_splits)-1)]

# %%
orders = range(1, 5)
fits = [('Cl/%', 'T/K'), ('T/K', 'Cl/%'), ('T/K', 'Cs/%')] # parameter pairs to fit (x axis, y axis)
fig, axes = plt.subplots(ncols=len(orders), nrows=len(fits), figsize=(6*len(orders), 5*len(fits)))
stat_data = []

for i, order in enumerate(orders):
    axes[0, i].set_title(f'{order} Order Fits')

    for j, (x, y) in enumerate(fits):
        left_ax = axes[j, i]
        right_ax = left_ax.twinx()
        left_ax.set_xlabel(x)
        left_ax.set_ylabel(y)
        left_ax.scatter(data[x], data[y], color='black', marker='x')
        
        ranges = Cl_ranges if 'Cl/%' in [x, y] else Cs_ranges
        for [T_min, T_max] in ranges:
            T_lowest = data.loc[np.argmin(np.abs(data['T/K']-T_min)), 'T/K']
            T_highest = data.loc[np.argmin(np.abs(data['T/K']-T_max)), 'T/K']
            subset = data[(data['T/K']>=T_lowest) & (data['T/K']<=T_highest)]
            left_ax.axvline(subset[x].min(), color='black', linestyle='--')
            left_ax.axvline(subset[x].max(), color='black', linestyle='--')
            
            fit = np.polynomial.Polynomial.fit(subset[x], subset[y], deg=order)
            y_fit = fit(subset[x])
            left_ax.plot(subset[x], y_fit, color='red', linestyle='--')
            right_ax.set_ylabel(f'Error')
            y_errors = np.abs(subset[y]-y_fit)
            average_x_size = (subset[x].max()-subset[x].min())/len(subset[x])
            right_ax.bar(subset[x], y_errors, color='grey', alpha=0.5, width=average_x_size)

            stat_row = dict()
            stat_row['kind'] = f'{x} to {y}'
            stat_row['order'] = order
            stat_row['x lowest'] = subset[x].min()
            stat_row['x highest'] = subset[x].max()
            coefs = fit.convert().coef
            for coef in range(max(orders)+1):
                stat_row[f'{coef} order coefficient'] = f'{coefs[coef]:.5E}' if coef<len(coefs) else '0'
            stat_row['mean error'] = np.mean(y_errors)
            stat_row['max error'] = np.max(y_errors)
            stat_data.append(stat_row)


fig.tight_layout()
fig.savefig(home_path / "fits.png")
df = pd.DataFrame(stat_data).sort_values(by=['kind', 'x lowest', 'order'])
df.to_csv(home_path / "fits.csv", index=False)

