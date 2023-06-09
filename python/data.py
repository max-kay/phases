import numpy as np
from matplotlib import cm
import pandas as pd
import numpy as np
from scipy.integrate import cumulative_trapezoid


def prepare_data(path: str, **kwarg) -> pd.DataFrame:
    df = pd.read_csv(path, header=1, **kwarg)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.sort_values(["c", "temp"], ascending=[True, True], inplace=True)
    add_entropy_and_free_energy(df)
    add_chemical_potential(df)
    return df


def add_entropy_and_free_energy(df: pd.DataFrame):
    df["entropy"] = np.nan
    df["free energy"] = np.nan

    def fn(group):
        entropy = cumulative_trapezoid(
            np.where(
                np.bitwise_or(group["temp"].values == 0, group["heat capacity"].isna()),
                0,
                group["heat capacity"].values / group["temp"].values,
            ),
            x=group["temp"].values,
            initial=0,
        )
        concentration = group.name
        entropy_ideal = -(
            concentration * np.log(concentration)
            + (1 - concentration) * np.log(1 - concentration)
        )
        constant = entropy_ideal - entropy[-1]
        df.loc[group.index, "entropy"] = entropy# + constant
        df.loc[group.index, "free energy"] = group["energy"] - group["temp"] * entropy

    df.groupby("c").apply(fn)


def add_chemical_potential(df: pd.DataFrame):
    df["mu a"] = np.nan

    def fn(group):
        mu_a = np.gradient(group["free energy"], group["c"])
        df.loc[group.index, "mu a"] = mu_a

    df.groupby("temp").apply(fn)


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


if __name__ == "__main__":
    df = prepare_data("out/logs_bin/data_2023-04-16_00-05.csv")
    print(df)
