#!/opt/homebrew/Caskroom/miniconda/base/envs/datasc python
import pretty_errors
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from sys import argv
from math import prod


if len(argv) >= 2:
    name = argv[1]
else:
    name = "b_2_t_2023-04-17_10-03"

with open(f"out/systems/{name}.txt") as file:
    lines = file.read_lines()
    energies = eval(lines[1])
    sites = prod([int(x) for x in lines[3].split(",")])

df = pd.read_csv(f"out/logs/{name}.csv", dtype=float, header=1)

idxs = np.linspace(0, len(df) / sites, len(df))
df["energy"] = df["energy"] / sites

fig = plt.figure()

ax = fig.add_subplot(221)
ax.plot(idxs, df["energy"])
ax.set_xlabel("Step")
ax.set_ylabel("E")

ax = fig.add_subplot(222)
ax.plot(df["temp"], df["energy"])
ax.set_xlabel("T")
ax.set_ylabel("E")

ax = fig.add_subplot(223)
ax.plot(idxs, df["temp"])
ax.set_xlabel("Step")
ax.set_ylabel("T")

plt.get_current_fig_manager().full_screen_toggle()

plt.show()
