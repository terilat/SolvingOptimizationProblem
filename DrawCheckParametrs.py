import numpy as np
import matplotlib.pyplot as plt

f = open("CheckParametrs.txt", "r")

t = np.zeros(1)
x0 = np.zeros(1)
x1 = np.zeros(1)
x2 = np.zeros(1)
x3 = np.zeros(1)

for l in f:
    q = l.split()
    t = np.append(t, float(q[0]))
    x0 = np.append(x0, float(q[1]))
    x1 = np.append(x1, float(q[2]))
    x2 = np.append(x2, float(q[3]))
    x3 = np.append(x3, float(q[4]))
t = np.delete(t, 0)
x0 = np.delete(x0, 0)
x1 = np.delete(x1, 0)
x2 = np.delete(x2, 0)
x3 = np.delete(x3, 0)

plt.plot(t, x0, 'g')
plt.title("x(t)")
plt.xlabel("t")
plt.ylabel("x")
plt.savefig("x0.png", dpi = 200)
plt.show()

plt.plot(t, x1, 'g')
plt.title("y(t)")
plt.xlabel("t")
plt.ylabel("y")
plt.savefig("x1.png", dpi = 200)
plt.show()

plt.plot(t, x2, 'g')
plt.title("$p_x(t)$")
plt.xlabel("t")
plt.ylabel("$p_x$")
plt.savefig("x2.png", dpi = 200)
plt.show()

plt.plot(t, x3, 'g')
plt.title("$p_y(t)$")
plt.xlabel("t")
plt.ylabel("$p_y$")
plt.savefig("x3.png", dpi = 200)
plt.show()
