import numpy as np
import pandas as pd
from matplotlib import cm


def get_entropy_and_free_energy(df: pd.DataFrame) -> pd.DataFrame:
    cs = df["c"].unique()

    over_all_entropy = []

    for c in cs:
        dfi = df[df["c"] == c]
        temps = dfi["temp"].values
        cs = dfi["c"].values
        # max_entropy = -(c * np.log(c) + (1 - c) * np.log(1 - c))
        intergral = 0
        entropy = [intergral]
        for i in range(0, len(cs) - 1):
            if temps[i] == 0:
                print("skipped one step in intgration for entropy")
                entropy.append(intergral)
                continue
            intergral += (
                (temps[i + 1] - temps[i])
                / 2
                * (cs[i + 1] / temps[i + 1] + cs[i] / temps[i])
            )
            entropy.append(intergral)
        over_all_entropy += entropy

    df["entropy"] = over_all_entropy

    df["free energy"] = df["energy"] - df["temp"] * df["entropy"]
    return df


def get_chemical_potential(df: pd.DataFrame) -> pd.DataFrame:
    temps = df["temp"].unique()

    over_all_mu = []
    for t in temps:
        dfi = df[df["temp"] == t]
        energies = dfi["energy"].values
        cs = dfi["c"].values
        over_all_mu.append(0)
        for i in range(0, len(cs) - 1):
            energies[i] 
        over_all_entropy += entropy

    df["mu"] = over_all_mu
    return df


def plot_grid_3d(xs, ys, zs, ax_3d, cmap=cm.coolwarm):
    n_xs = len(xs.unique())
    ax_3d.plot_surface(
        np.array(xs).reshape((n_xs, -1)),
        np.array(ys).reshape((n_xs, -1)),
        np.array(zs).reshape((n_xs, -1)),
        cmap=cmap,
        linewidth=0,
        antialiased=False,
    )
