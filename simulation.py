import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys


def simulation(time, N_S, N_P, N_A, N_Y, N_R, s, p, a, y, p_theta, a_theta, p_kappa, a_kappa, y_kappa, total, gamma, beta, phi, steps, cycles):
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


def get_average(time, N_S0, N_P0, N_A0, N_Y0, N_R0, N_S, N_P, N_A, N_Y, N_R, avg_steps, cycles, steps, conf_level):
    time_max = time.max()
    Time_avg = np.linspace(0, time_max, avg_steps)
    N_S_avg = np.zeros(avg_steps)
    N_P_avg = np.zeros(avg_steps)
    N_A_avg = np.zeros(avg_steps)
    N_Y_avg = np.zeros(avg_steps)
    N_R_avg = np.zeros(avg_steps)

    N_S_conf = np.zeros(avg_steps, dtype=object)
    N_P_conf = np.zeros(avg_steps, dtype=object)
    N_A_conf = np.zeros(avg_steps, dtype=object)
    N_Y_conf = np.zeros(avg_steps, dtype=object)
    N_R_conf = np.zeros(avg_steps, dtype=object)

    N_S_avg[0] = N_S0
    N_P_avg[0] = N_P0
    N_A_avg[0] = N_A0
    N_Y_avg[0] = N_Y0
    N_R_avg[0] = N_R0

    for i in range(0, avg_steps):
        time_max = Time_avg[i]
        S_slice = []
        P_slice = []
        A_slice = []
        Y_slice = []
        R_slice = []
        total_count = 0

        for j in range(cycles):
            S_sum = 0
            P_sum = 0
            A_sum = 0
            Y_sum = 0
            R_sum = 0
            step_count = 0

            for k in range(steps):
                if time[j, k] <= time_max and time[j, k + 1] > time_max:
                    step_count += 1
                    total_count += 1
                    S_sum += N_S[j, k]
                    P_sum += N_P[j, k]
                    A_sum += N_A[j, k]
                    Y_sum += N_Y[j, k]
                    R_sum += N_R[j, k]

            
            if step_count != 0:
                S_slice.append(S_sum)
                P_slice.append(P_sum)
                A_slice.append(A_sum)
                Y_slice.append(Y_sum)
                R_slice.append(R_sum)

        if total_count == 0:
            N_S_avg[i] = N_S_avg[i-1]
            N_P_avg[i] = N_P_avg[i-1]
            N_A_avg[i] = N_A_avg[i-1]
            N_Y_avg[i] = N_Y_avg[i-1]
            N_R_avg[i] = N_R_avg[i-1]

            N_S_conf[i] = N_S_conf[i-1]
            N_P_conf[i] = N_P_conf[i-1]
            N_A_conf[i] = N_A_conf[i-1]
            N_Y_conf[i] = N_Y_conf[i-1]
            N_R_conf[i] = N_R_conf[i-1]
        else:
            N_S_avg[i] = np.mean(S_slice)
            N_P_avg[i] = np.mean(P_slice)
            N_A_avg[i] = np.mean(A_slice)
            N_Y_avg[i] = np.mean(Y_slice)
            N_R_avg[i] = np.mean(R_slice)

            N_S_conf[i] = stats.t.interval(conf_level, df=len(
                S_slice)-1, loc=np.mean(S_slice), scale=np.std(S_slice, ddof=1) / np.sqrt(len(S_slice)))
            N_P_conf[i] = stats.t.interval(conf_level, df=len(
                P_slice)-1, loc=np.mean(P_slice), scale=np.std(P_slice, ddof=1) / np.sqrt(len(P_slice)))
            N_A_conf[i] = stats.t.interval(conf_level, df=len(
                A_slice)-1, loc=np.mean(A_slice), scale=np.std(A_slice, ddof=1) / np.sqrt(len(A_slice)))
            N_Y_conf[i] = stats.t.interval(conf_level, df=len(
                Y_slice)-1, loc=np.mean(Y_slice), scale=np.std(Y_slice, ddof=1) / np.sqrt(len(Y_slice)))
            N_R_conf[i] = stats.t.interval(conf_level, df=len(
                R_slice)-1, loc=np.mean(R_slice), scale=np.std(R_slice, ddof=1) / np.sqrt(len(R_slice)))


    return (Time_avg, (N_S_avg, N_S_conf), (N_P_avg, N_P_conf), (N_A_avg, N_A_conf), (N_Y_avg, N_Y_conf), (N_R_avg, N_R_conf))


def visualizer(Time_avg, N_S_res, N_P_res, N_A_res, N_Y_res, N_R_res, plot_name):
    plt.ylabel("Population")
    plt.xlabel("time")

    plt.plot(Time_avg, N_S_res[0], marker="",
             color="red", linewidth=1.9, alpha=0.9)

    plt.plot(Time_avg, N_P_res[0], marker="",
             color="blue", linewidth=1.9, alpha=0.9)

    plt.plot(Time_avg, N_A_res[0], marker="",
             color="orange", linewidth=1.9, alpha=0.9)

    plt.plot(Time_avg, N_Y_res[0], marker="",
             color="yellow", linewidth=1.9, alpha=0.9)

    plt.plot(Time_avg, N_R_res[0], marker="",
             color="pink", linewidth=1.9, alpha=0.9)

    plt.legend(["Susceptible",
                "Pre-symptomatic", "Asymptomatic", "Symptomatic", "Recovered"])

    plt.savefig(plot_name)


def data_collector(N_P_avg, N_A_avg, N_Y_avg, avg_steps):
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
    conf_level = 0.95

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
    Time_avg, N_S_res, N_P_res, N_A_res, N_Y_res, N_R_res = get_average(time, N_S0, N_P0, N_A0, N_Y0, N_R0, N_S,
                                                                        N_P, N_A, N_Y, N_R, avg_steps, cycles, steps, conf_level)

    # visualizer(Time_avg, N_S_res, N_P_res, N_A_res,
               # N_Y_res, N_R_res, "plot.png")
    # data_collector(N_P_avg, N_A_avg, N_Y_avg, avg_steps)


if __name__ == '__main__':
    main()
