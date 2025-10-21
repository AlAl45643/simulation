import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
import sys
np.set_printoptions(threshold=sys.maxsize)


def resize_np_array(array):
    pad_size = (int((array.shape[1] - 1) * 2))
    return np.pad(array, ((0, 0), (0, (pad_size))), 'constant', constant_values=0)


def simulation(time, N_S, N_P, N_A, N_Y, N_R, s, p, a, y, p_theta, a_theta, p_kappa, a_kappa, y_kappa, total, gamma, beta, phi, steps, cycles):
    for i in range(cycles):
        j = 0
        while N_P[i, j] + N_Y[i, j] + N_A[i, j] != 0:
            if j > N_S.shape[1] - 2:
                N_S = resize_np_array(N_S)
                N_P = resize_np_array(N_P)
                N_A = resize_np_array(N_A)
                N_Y = resize_np_array(N_Y)
                N_R = resize_np_array(N_R)
                time = resize_np_array(time)
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
            j += 1
    return (N_S, N_P, N_A, N_Y, N_R, time)


def clean_up_conf(confs, avgs):
    lower_conf = []
    higher_conf = []
    for idx, conf in enumerate(confs):
        if np.isnan(conf[0]):
            lower_conf.append(avgs[idx])
        else:
            lower_conf.append(conf[0])

    for idx, conf in enumerate(confs):
        if np.isnan(conf[1]):
            higher_conf.append(avgs[idx])
        else:
            higher_conf.append(conf[1])

    return (lower_conf, higher_conf)


def get_avg_and_conf(time, N_S0, N_P0, N_A0, N_Y0, N_R0, N_S, N_P, N_A, N_Y, N_R, avg_steps, cycles, steps, conf_level):
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

            for k in range(N_S.shape[1] - 1):
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

    N_S_lower_conf, N_S_higher_conf = clean_up_conf(N_S_conf, N_S_avg)
    N_P_lower_conf, N_P_higher_conf = clean_up_conf(N_P_conf, N_P_avg)
    N_A_lower_conf, N_A_higher_conf = clean_up_conf(N_A_conf, N_A_avg)
    N_Y_lower_conf, N_Y_higher_conf = clean_up_conf(N_Y_conf, N_Y_avg)
    N_R_lower_conf, N_R_higher_conf = clean_up_conf(N_R_conf, N_R_avg)

    return (Time_avg, (N_S_avg, N_S_lower_conf, N_S_higher_conf), (N_P_avg, N_P_lower_conf, N_P_higher_conf), (N_A_avg, N_A_lower_conf, N_A_higher_conf), (N_Y_avg, N_Y_lower_conf, N_Y_higher_conf), (N_R_avg, N_R_lower_conf, N_R_higher_conf))


def plot_avg_confs(axs, title, x, y, low_conf_y, high_conf_y, color):
    axs.set_title(title)
    axs.plot(x, y, marker="", color=color, linewidth=0.5, alpha=0.9)
    axs.plot(x, low_conf_y, marker=".",
             color="black", linewidth=0.5, alpha=0.5)
    axs.plot(x, high_conf_y, marker=".",
             color="black", linewidth=0.5, alpha=0.5)
    axs.set_xlim(0, max(x))


def visualizer(Time_avg, N_S_res, N_P_res, N_A_res, N_Y_res, N_R_res, plot_name):
    plt.ylabel("Population")
    plt.xlabel("time")

    fig, axs = plt.subplots(5, 1, figsize=(6, 20))
    plot_avg_confs(axs[0], "Susceptible", Time_avg,
                   N_S_res[0], N_S_res[1], N_S_res[2], "red")
    plot_avg_confs(axs[1], "Pre-symptomatic", Time_avg,
                   N_P_res[0], N_P_res[1], N_P_res[2], "blue")
    plot_avg_confs(axs[2], "Asymptomatic", Time_avg,
                   N_A_res[0], N_A_res[1], N_A_res[2], "orange")
    plot_avg_confs(axs[3], "Symptomatic", Time_avg,
                   N_Y_res[0], N_Y_res[1], N_Y_res[2], "yellow")
    plot_avg_confs(axs[4], "Recovered", Time_avg,
                   N_R_res[0], N_R_res[1], N_R_res[2], "pink")

    plt.savefig(plot_name)


def data_collector(Time_avg, N_S_avg, N_P_avg, p_kappa, p_theta, N_A_avg, a_kappa, a_theta, N_Y_avg, N_R_avg, y_kappa, avg_steps, s, p, a, y, gamma, beta, phi, total):
    peak_infections = 0
    peak_time = 0
    for i in range(0, avg_steps):
        num_infect = N_P_avg[i] + N_A_avg[i] + N_Y_avg[i]
        if num_infect > peak_infections:
            peak_infections = num_infect
            peak_time = Time_avg[i]

    R0 = ((p_kappa * p + a_kappa * a + y_kappa * y) * beta) / \
        ((p_theta * (gamma + phi)) + (a_theta * (gamma)))

    attack_rate = N_R_avg[len(N_R_avg) - 1] / N_S_avg[0]
    return (("peak_infections", peak_infections), ("peak_time", peak_time), ("R0", R0), ("attack_rate", attack_rate))


def export_csv(metrics_name, metrics, *arrays):
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(metrics_name)
    for array in arrays:
        array_df = pd.DataFrame(array[0])
        array_df.to_csv(array[1])


def main():
    np.random.seed = 312342
    N_S0 = 5000
    N_P0 = 1
    N_A0 = 0
    N_Y0 = 0
    N_R0 = 0
    total = N_S0 + N_P0 + N_A0 + N_Y0 + N_R0
    s = 0.5
    p = 0.7
    a = 0.7
    y = 0.7
    p_kappa = 0.7
    a_kappa = 0.7
    y_kappa = 0.7
    p_theta = 0.8
    a_theta = 0.2
    beta = 0.6
    gamma = 0.2
    phi = 0.3
    steps = 10
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

    N_S, N_P, N_A, N_Y, N_R, time = simulation(time, N_S, N_P, N_A, N_Y, N_R, s, p, a, y, p_theta, a_theta,
                                               p_kappa, a_kappa, y_kappa, total, gamma, beta, phi, steps, cycles)

    avg_steps = 50
    Time_avg, N_S_res, N_P_res, N_A_res, N_Y_res, N_R_res = get_avg_and_conf(time, N_S0, N_P0, N_A0, N_Y0, N_R0, N_S,
                                                                             N_P, N_A, N_Y, N_R, avg_steps, cycles, 200, conf_level)

    visualizer(Time_avg, N_S_res, N_P_res, N_A_res,
               N_Y_res, N_R_res, "plot.png")
    metrics = data_collector(
        Time_avg, N_S_res[0], N_P_res[0], p_kappa, p_theta, N_A_res[0], a_kappa, a_theta, N_Y_res[0], N_R_res[0], y_kappa, avg_steps, s, p, a, y, gamma, beta, phi, total)
    export_csv("metrics.csv", metrics, (N_S, "N_S.csv"), (N_P, "N_P.csv"),
               (N_A, "N_A.csv"), (N_Y, "N_Y.csv"), (N_R, "N_R.csv"))


if __name__ == '__main__':
    main()
