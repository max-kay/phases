#!/opt/homebrew/Caskroom/miniconda/base/envs/datasc
import numpy as np
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt

import data

if len(argv) == 2:
    name = argv[1]
else:
    name = "b_2_f_2023-04-17_14-21"

df = data.prepare_data(f"out/logs/{name}.csv")

df = df[df["temp"] < 100]

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

ax.set_title("Free Energy")
ax.set_xlabel("$X_a$")
ax.set_ylabel("T")
ax.set_zlabel("F")


ax = fig.add_subplot(234, projection="3d")
data.plot_grid_3d(df["c"], df["temp"], np.log(df["heat capacity"]), ax)

ax.set_title("Heat Capacity")
ax.set_xlabel("$X_a$")
ax.set_ylabel("T")
ax.set_zlabel("ln(C)")


ax = fig.add_subplot(235, projection="3d")
data.plot_grid_3d(df["c"], df["temp"], df["entropy"], ax)

ax.set_title("Entropy")
ax.set_xlabel("$X_a$")
ax.set_ylabel("T")
ax.set_zlabel("S")

ax = fig.add_subplot(236, projection="3d")
data.plot_grid_3d(df["c"], df["temp"], df["mu a"], ax)

ax.set_title("Chemical Potential")
ax.set_xlabel("$X_a$")
ax.set_ylabel("T")
ax.set_zlabel("$\mu _a$")

plt.get_current_fig_manager().full_screen_toggle()

plt.show()
