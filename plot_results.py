import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


df = pd.read_csv (r'results.csv')
states = df.to_numpy()

time = states[:,0]
phi = states[:,1]
theta = states[:,2]
psi = states[:,3]

plt.figure("Attitude Estimation (using EKF)", figsize=(16,9))
plt.subplot(3,1,1)
plt.plot(time, phi)
plt.xlabel("Time (sec)")
plt.ylabel("phi (deg)")
plt.grid("true")
plt.subplot(3,1,2)
plt.plot(time, theta)
plt.xlabel("Time (sec)")
plt.ylabel("theta (deg)")
plt.grid("true")
plt.subplot(3,1,3)
plt.plot(time, psi)
plt.xlabel("Time (sec)")
plt.ylabel("psi (deg)")
plt.grid("true")
plt.show()
