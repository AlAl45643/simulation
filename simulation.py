import numpy as np
import matplotlib.pyplot as plt
np.random.seed(150)

N_S0 = 500
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

N_R0 = 0

beta = 1.5
gamma = 0.2
phi = 0.3

total = N_S0 + N_P0 + N_A0 + N_Y0

steps = 1500
cycles = 50

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
N_R[:, 0] = N_R0

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
            N_S[i, j+1:] = N_S[i, j]
            N_P[i, j+1:] = N_P[i, j]
            N_A[i, j+1:] = N_A[i, j]
            N_Y[i, j+1:] = N_Y[i, j]
            N_R[i, j+1:] = N_R[i, j]
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

avg_steps = 100
Time_max = Time.max()

Time_avg = np.linspace(0, Time_max, avg_steps)
N_S_avg = np.zeros(avg_steps)
N_P_avg = np.zeros(avg_steps)
N_A_avg = np.zeros(avg_steps)
N_Y_avg = np.zeros(avg_steps)
N_R_avg = np.zeros(avg_steps)

N_S_avg[0] = N_S0
N_P_avg[0] = N_P0
N_A_avg[0] = N_A0
N_Y_avg[0] = N_Y0
N_R_avg[0] = N_R0

for i in range(0, avg_steps):
    time_max = Time_avg[i]
    S_sum = 0
    P_sum = 0
    A_sum = 0
    Y_sum = 0
    R_sum = 0
    t_count = 0

    for j in range(cycles):
        for k in range(steps):
            if Time[j, k] <= time_max and Time[j, k + 1] > time_max:
                t_count += 1
                S_sum += N_S[j, k]
                P_sum += N_P[j, k]
                A_sum += N_A[j, k]
                Y_sum += N_Y[j, k]
                R_sum += N_R[j, k]

        if t_count == 0:
            N_S_avg[i] = N_S_avg[i-1]
            N_P_avg[i] = N_P_avg[i-1]
            N_A_avg[i] = N_A_avg[i-1]
            N_Y_avg[i] = N_Y_avg[i-1]
            N_R_avg[i] = N_R_avg[i-1]
        else:
            N_S_avg[i] = S_sum / t_count
            N_P_avg[i] = P_sum / t_count
            N_A_avg[i] = A_sum / t_count
            N_Y_avg[i] = Y_sum / t_count
            N_R_avg[i] = R_sum / t_count


plt.ylabel("Population")
plt.xlabel("Time")

plt.plot(Time_avg, N_S_avg, marker="",
         color="red", linewidth=1.9, alpha=0.9)

plt.plot(Time_avg, N_P_avg, marker="",
         color="blue", linewidth=1.9, alpha=0.9)

plt.plot(Time_avg, N_A_avg, marker="",
         color="orange", linewidth=1.9, alpha=0.9)

plt.plot(Time_avg, N_Y_avg, marker="",
         color="yellow", linewidth=1.9, alpha=0.9)

plt.plot(Time_avg, N_R_avg, marker="",
         color="pink", linewidth=1.9, alpha=0.9)

plt.legend(["Susceptible",
            "Pre-symptomatic", "Asymptomatic", "Symptomatic", "Recovered"])

plt.savefig("plot.png")

# Peak number of infections at a time
peak_infections = 0
for i in range(0, avg_steps):
    num_infect = N_P_avg[i] + N_A_avg[i] + N_Y_avg[i]
    if num_infect > peak_infections:
        peak_infections = num_infect


print(f"Peak infections at a point in time: {peak_infections}")


