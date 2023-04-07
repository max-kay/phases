import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

# df = pd.read_csv("logs/data_2023-04-06_15-52.csv")
# df = pd.read_csv("logs/data_2023-04-06_20-18.csv")
# df = pd.read_csv("logs/data_2023-04-06_22-57.csv") ##
# df = pd.read_csv("logs/data_2023-04-06_22-59.csv")
# df = pd.read_csv("logs/data_2023-04-06_23-01.csv")
# df = pd.read_csv("logs/data_2023-04-07_00-13.csv")
# df = pd.read_csv("logs/data_2023-04-07_01-36.csv")
# df = pd.read_csv("logs/data_2023-04-07_09-08.csv")
# temperature was chosen in the wrong direction until here
# df = pd.read_csv("logs/data_2023-04-07_09-53.csv")
# df = pd.read_csv("logs/data_2023-04-07_10-04.csv")
df = pd.read_csv("logs/data_2023-04-07_12-36.csv")


n_concentrations = len(df["concentration a"].unique())

fig = plt.figure()
ax1 = fig.add_subplot(121, projection="3d")
ax1.plot_surface(
    df["concentration a"].values.reshape((n_concentrations, -1)),
    df["temperature"].values.reshape((n_concentrations, -1)),
    df["internal energy U"].values.reshape((n_concentrations, -1)),
    cmap=cm.coolwarm,
    linewidth=0,
    antialiased=False,
)


ax1.set_xlabel("$X_a$")
ax1.set_ylabel("T")
ax1.set_zlabel("U")
ax1.set_title("Internal energy")

ax2 = fig.add_subplot(122, projection="3d")

ax2.plot_surface(
    df["concentration a"].values.reshape((n_concentrations, -1)),
    df["temperature"].values.reshape((n_concentrations, -1)),
    df["heat capacity"].values.reshape((n_concentrations, -1)),
    cmap=cm.coolwarm,
    linewidth=0,
    antialiased=False,
)

ax2.set_xlabel("$X_a$")
ax2.set_ylabel("T")
ax2.set_zlabel("C")
ax2.set_title("heat capacity")

plt.show()
