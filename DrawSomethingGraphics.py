import numpy as np
import matplotlib.pyplot as plt

f = open("Res.txt", "r")

t = np.zeros(1)
y = np.zeros(1)
for l in f:
    q = l.split()
    t = np.append(t, float(q[0]))
    y = np.append(y, float(q[1]))
t = np.delete(t, 0)
y = np.delete(y, 0)


plt.plot(t, y, 'g')
plt.savefig("RungeKuttaForOrbit.png", dpi = 200)
plt.show()

