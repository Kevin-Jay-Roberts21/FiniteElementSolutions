import math

import matplotlib.pyplot as plt
import numpy as np
import matplotlib

# defining the global vars
L = 5 # length in cm
dx = 0.2 # space step
M = int(L/dx) # total number of space steps
Y_0 = 0
Y_M = 2

# Exact Analytical Solution
###########################
x_exact = np.linspace(0, L, 10000)  # 10,000 points → visually continuous
exact_Y = 2 * np.sin(2 * x_exact) / np.sin(2 * L)


# Finite Difference Approximation
#################################

# initializing the solution vector of length M
fd_Y = np.zeros(M+1)
fd_Y[0] = Y_0
fd_Y[M] = Y_M

a = 4*dx*dx - 2

A = np.zeros((M-1, M-1))
b = np.zeros(M-1)
b[0] = -Y_0
b[-1] = -Y_M
A[0][0] = a
A[0][1] = 1
A[M-2][M-3] = 1
A[M-2][M-2] = a

for i in range(1, M-2):
    A[i][i] = a
    A[i][i-1] = 1
    A[i][i+1] = 1

fd_Y_inner = np.linalg.solve(A, b)
fd_Y[1:-1] = fd_Y_inner


# Finite Element Approximation
##############################








# Plotting the solution and approximations
x = np.linspace(0, L, M+1)
plt.plot(x_exact, exact_Y, label="Exact Solution")
plt.plot(x, fd_Y, label="Finite Difference Solution")



plt.xlabel('x in cm')
plt.ylabel('y(x) in cm')
plt.title('Comparison of Solutions')
plt.legend()
plt.grid(True)
plt.show()