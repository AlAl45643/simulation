import numpy as np
import matplotlib.pyplot as plt

def simulation(time: list, N_S: list, N_P: list, N_A: list, N_Y: list, N_R: list, s: float, p: float, a: float, y: float, p_theta: float, a_theta: float, p_kappa: float, a_kappa: float, y_kappa: float, total: float, gamma: float, beta: float, phi: float, steps: int, cycles: int):
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

            event_rate = P + A + Y + R_1 + R_2
            # Virus died
            if event_rate == 0:
                N_S[i, j+1:] = N_S[i, j]
                N_P[i, j+1:] = N_P[i, j]
                N_A[i, j+1:] = N_A[i, j]
                N_Y[i, j+1:] = N_Y[i, j]
                N_R[i, j+1:] = N_R[i, j]
                break

            u1 = np.random.random()
            tau = 1/event_rate * np.log(1/u1)
            time[i, j+1] = time[i, j] + tau

            event = np.random.choice(["P", "A", "Y", "R_1", "R_2"], p=[
                P/event_rate, A/event_rate, Y/event_rate, R_1/event_rate, R_2/event_rate])
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


def get_average(time: list, N_S0: int, N_P0: int, N_A0: int, N_Y0: int, N_R0: int, N_S: list, N_P: list, N_A: list, N_Y: list, N_R: list, avg_steps: int, cycles: int, steps: int) -> (list, list, list, list, list, list):
    time_max = time.max()
    Time_avg = np.linspace(0, time_max, avg_steps)
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
                if time[j, k] <= time_max and time[j, k + 1] > time_max:
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

    return (Time_avg, N_S_avg, N_P_avg, N_A_avg, N_Y_avg, N_R_avg)


def visualizer(Time_avg: list, N_S_avg: list, N_P_avg: list, N_A_avg: list, N_Y_avg: list, N_R_avg: list, plot_name):
    plt.ylabel("Population")
    plt.xlabel("time")

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

    plt.savefig(plot_name)


def data_collector(N_P_avg: list, N_A_avg: list, N_Y_avg: list, avg_steps):
    # Peak number of infections at a time
    peak_infections = 0
    for i in range(0, avg_steps):
        num_infect = N_P_avg[i] + N_A_avg[i] + N_Y_avg[i]
        if num_infect > peak_infections:
            peak_infections = num_infect
    print(peak_infections)


def main():
    np.random.seed = 312342
    N_S0 = 500
    N_P0 = 1
    N_A0 = 1
    N_Y0 = 1
    N_R0 = 0
    total = N_S0 + N_P0 + N_A0 + N_Y0 + N_R0
    s = 0.7
    p = 0.7
    a = 0.7
    y = 0.9
    p_kappa = 0.2
    a_kappa = 0.7
    y_kappa = 0.9
    p_theta = 0.7
    a_theta = 0.3
    beta = 1.5
    gamma = 0.2
    phi = 0.3
    steps = 1500
    cycles = 50

    time = np.zeros((cycles, steps+1))
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

    simulation(time, N_S, N_P, N_A, N_Y, N_R, s, p, a, y, p_theta, a_theta,
               p_kappa, a_kappa, y_kappa, total, gamma, beta, phi, steps, cycles)

    avg_steps = 50
    Time_avg, N_S_avg, N_P_avg, N_A_avg, N_Y_avg, N_R_avg = get_average(time, N_S0, N_P0, N_A0, N_Y0, N_R0, N_S,
                                                                        N_P, N_A, N_Y, N_R, avg_steps, cycles, steps)

    visualizer(Time_avg, N_S_avg, N_P_avg, N_A_avg,
               N_Y_avg, N_R_avg, "plot.png")
    data_collector(N_P_avg, N_A_avg, N_Y_avg, avg_steps)


if __name__ == '__main__':
    main()
