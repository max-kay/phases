#!/opt/homebrew/Caskroom/miniconda/base/envs/datasc
import pretty_errors
import numpy as np
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt

import data

if len(argv) == 2:
    df = pd.read_csv(argv[1])
else:
    df = pd.read_csv("logs_bin/data_2023-04-16_00-05.csv")

df.sort_values(["c", "temp"], ascending=[True, True], inplace=True)
df = data.get_entropy_and_free_energy(df)


fig = plt.figure()

ax = fig.add_subplot(231, projection="3d")
data.plot_grid_3d(df["c"], df["temp"], df["energy"], ax)

ax.set_title("Internal Energy")
ax.set_xlabel("$X_a$")
ax.set_ylabel("T")
ax.set_zlabel("U")

ax = fig.add_subplot(232, projection="3d")
data.plot_grid_3d(df["c"], df["entropy"], df["energy"], ax)

ax.set_title("Internal Energy")
ax.set_xlabel("$X_a$")
ax.set_ylabel("S")
ax.set_zlabel("U")


ax = fig.add_subplot(233, projection="3d")
data.plot_grid_3d(df["c"], df["temp"], df["free energy"], ax)

ax.set_xlabel("$X_a$")
ax.set_ylabel("T")
ax.set_zlabel("F")
ax.set_title("Free Energy")


ax = fig.add_subplot(234, projection="3d")
data.plot_grid_3d(df["c"], df["temp"], np.log(df["heat capacity"]), ax)

ax.set_xlabel("$X_a$")
ax.set_ylabel("T")
ax.set_zlabel("ln(C)")
ax.set_title("Heat Capacity")


ax = fig.add_subplot(235, projection="3d")
data.plot_grid_3d(df["c"], df["temp"], df["entropy"], ax)

ax.set_xlabel("$X_a$")
ax.set_ylabel("T")
ax.set_zlabel("S")
ax.set_title("Entropy")


plt.show()
