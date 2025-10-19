import math

import matplotlib.pyplot as plt
import numpy as np

# defining the global vars
L = 10 # length in cm
dx = 0.05 # space step
M = int(round(L/dx)) # total number of space steps
dx = L/M # ensuring spacing matches the coordinates
Y_0 = 0
Y_M = 2
FE_x = np.arange(0, L + 1e-12, dx)
FD_x = FE_x

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

# defining the M, K matrices and F vector
FE_M = np.zeros((M-1, M-1))
FE_K = np.zeros((M-1, M-1))
FE_F = (8*dx/L)*FE_x[1:-1]

FE_K[0][0] = 2
FE_K[0][1] = -1
FE_K[M-2][M-3] = -1
FE_K[M-2][M-2] = 2
FE_M[0][0] = 4
FE_M[0][1] = 1
FE_M[M-2][M-3] = 1
FE_M[M-2][M-2] = 4

for i in range(1, M-2):
    FE_K[i][i-1] = -1
    FE_K[i][i] = 2
    FE_K[i][i+1] = -1
    FE_M[i][i-1] = 1
    FE_M[i][i] = 4
    FE_M[i][i+1] = 1

FE_K = (1/dx)*FE_K
FE_M = (dx/6)*FE_M

# solving for Y interior
FE_Y = np.linalg.solve((FE_K - 4*FE_M), FE_F)

# defining the boundary function
def y_b(x):
    return 2*x/L

def phi_i(x, i):
    if i < 0 or i > M:
        return 0
    x_i = FE_x[i]

    if i > 0:
        x_i_minus_1 = FE_x[i-1]
        if x_i_minus_1 <= x < x_i:
            return (x - x_i_minus_1)/dx
    if i < M:
        x_i_plus_1 = FE_x[i+1]
        if x_i <= x <= x_i_plus_1:
            return (x_i_plus_1 - x)/dx
    return 0



# defining the finite element summation
def y_dx(x):
    sum = 0
    for i in range(1, M):
        sum += FE_Y[i-1] * phi_i(x, i)
    return sum

# defining the final FE approximation
fe_Y = np.zeros(M+1)
fe_Y[0] = Y_0
fe_Y[M] = Y_M
for i in range(1, M):
    x_i = i*dx
    fe_Y[i] = y_dx(x_i) + y_b(x_i)


# Plotting the solution and approximations
plt.plot(x_exact, exact_Y, label="Exact Solution")
plt.plot(FD_x, fd_Y, label="Finite Difference Solution")
plt.plot(FE_x, fe_Y, label="Finite Element Solution")


plt.xlabel('x in cm')
plt.ylabel('y(x) in cm')
plt.title('Comparison of Solutions')
plt.legend()
plt.grid(True)
plt.show()