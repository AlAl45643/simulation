import numpy as np
import matplotlib.pyplot as mp

np.random.seed(150)


N_S0 = 10000
s = 0.7

N_P0 = 1
p = 0.7
p_kappa = 0.2
p_theta = 0.7

N_A0 = 1
a = 0.7
a_kappa = 0.1
a_theta = 0.3

N_Y0 = 1
y = 0.9
y_kappa = 0.3


beta = 0.3
gamma = 0.3
phi = 0.3

total = N_S0 + N_P0 + N_A0 + N_Y0

steps = 25
cycles = 100

Time = np.zeros((cycles, steps+1))
N_S = np.zeros((cycles, steps+1))
N_P = np.zeros((cycles, steps+1))
N_A = np.zeros((cycles, steps+1))
N_Y = np.zeros((cycles, steps+1))
N_R = np.zeros((cycles, steps+1))

N_S[:, 0] = N_S0
N_P[:, 0] = N_P0
N_A[:, 0] = N_A0
N_Y[:, 0] = N_Y0

for i in range(cycles):
    for j in range(steps):
        S_rate = N_S[i, j] / total
        P_rate = N_P[i, j] / total
        A_rate = N_A[i, j] / total
        Y_rate = N_Y[i, j] / total
        R_rate = N_R[i, j] / total
        T = (s * S_rate) + (p * P_rate) + \
            (a * A_rate) + (y * Y_rate) + (R_rate)

        # S = S - 1
        # P = P + 1
        P = p_theta * ((s * (p + a + y) / T)) * beta * (S_rate) * \
            ((p_kappa * P_rate) + (a_kappa * A_rate) + (y_kappa * Y_rate))
        # S = S - 1
        # A = A + 1
        A = a_theta * ((s * (p + a + y) / T)) * beta * (S_rate) \
            * ((p_kappa * P_rate) + (a_kappa * A_rate)
               + (y_kappa * Y_rate))

        # P = P - 1
        # Y = Y + 1
        Y = phi * P_rate

        # A = A - 1
        # R = R + 1
        R_1 = gamma * A_rate

        # Y = Y - 1
        # R = R + 1
        R_2 = gamma * Y_rate

        Event_rate = P + A + Y + R_1 + R_2
        # Virus died
        if Event_rate == 0:
            break

        u1 = np.random.random()
        tau = 1/Event_rate * np.log(1/u1)
        Time[i, j+1] = Time[i, j] + tau

        event = np.random.choice(["P", "A", "Y", "R_1", "R_2"], p=[
                                 P/Event_rate, A/Event_rate, Y/Event_rate, R_1/Event_rate, R_2/Event_rate])
        if event == "P":
            N_S[i, j+1] = N_S[i, j] - 1
            N_P[i, j+1] = N_P[i, j] + 1
            N_A[i, j+1] = N_A[i, j]
            N_Y[i, j+1] = N_Y[i, j]
            N_R[i, j+1] = N_R[i, j]
        elif event == "A":
            N_S[i, j+1] = N_S[i, j] - 1
            N_P[i, j+1] = N_P[i, j]
            N_A[i, j+1] = N_A[i, j] + 1
            N_Y[i, j+1] = N_Y[i, j]
            N_R[i, j+1] = N_R[i, j]
        elif event == "Y":
            N_S[i, j+1] = N_S[i, j]
            N_P[i, j+1] = N_P[i, j] - 1
            N_A[i, j+1] = N_A[i, j]
            N_Y[i, j+1] = N_Y[i, j] + 1
            N_R[i, j+1] = N_R[i, j]
        elif event == "R_1":
            N_S[i, j+1] = N_S[i, j]
            N_P[i, j+1] = N_P[i, j]
            N_A[i, j+1] = N_A[i, j] - 1
            N_Y[i, j+1] = N_Y[i, j]
            N_R[i, j+1] = N_R[i, j] + 1
        else:
            N_S[i, j+1] = N_S[i, j]
            N_P[i, j+1] = N_P[i, j]
            N_A[i, j+1] = N_A[i, j]
            N_Y[i, j+1] = N_Y[i, j] - 1
            N_R[i, j+1] = N_R[i, j] + 1

