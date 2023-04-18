#!/opt/homebrew/Caskroom/miniconda/base/envs/datasc python
import pretty_errors
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from sys import argv
from math import prod


def get_title(name: str) -> str:
    if "b_2_" in name:
        return "Cooling Curves for a 2D System"
    if "b_3_" in name:
        return "Cooling Curves for a 3D System"


if len(argv) >= 2:
    name = argv[1]
else:
    name = "b_2_t_2023-04-18_15-22"

with open(f"out/systems/{name}.txt") as file:
    lines = file.read()


df = pd.read_csv(f"out/logs/{name}.csv", dtype=float, header=1)

df["energy"] = df["energy"]

fig = plt.figure(figsize=(7 * 1.5, 5 * 1.5))

fig.suptitle(get_title(name))

ax = fig.add_subplot(221)
ax.plot(df["step"], df["energy"])
ax.set_xlabel("Step")
ax.set_ylabel("E")

ax = fig.add_subplot(222)
ax.plot(df["temp"], df["energy"])
ax.set_xlabel("T")
ax.set_ylabel("E")

ax = fig.add_subplot(223)
ax.plot(df["step"], df["temp"])
ax.set_xlabel("Step")
ax.set_ylabel("T")

fig.text(0.55, 0.4, f"{lines}", ha='left', va='top', fontsize=9)


plt.savefig(f"out/figs/{name}.png", dpi=300)
