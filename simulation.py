import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
import argparse


class Simulation:
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC7825872/
    p_kappa = 0.07
    a_kappa = 0.01
    y_kappa = 0.06
    # https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1003987&utm_source=nl_landingpage&utm_medium=email&utm_campaign=timestop10_daily_newsletter
    p_theta = 0.81
    a_theta = 0.19
    # https://www.nature.com/articles/s41467-021-20990-2/figures/1
    beta = 3.333333
    # https://www.sciencedirect.com/science/article/pii/S2468042721000038
    gamma = 0.2
    # https://academic.oup.com/cid/article/74/9/1678/6359063
    phi = 0.7
    # https://www.sciencedirect.com/science/article/pii/S2468042721000038 (again)
    zeta = 0.28
    start_steps = 10
    conf_level = 0.95

    def __init__(self, seed, N_S0, N_P0, N_A0, N_Y0, N_R0, s, p, a, y, cycles, avg_steps):
        np.random.seed = seed
        self.N_S0 = N_S0
        self.N_P0 = N_P0
        self.N_A0 = N_A0
        self.N_Y0 = N_Y0
        self.N_R0 = N_R0
        self.total = N_S0 + N_P0 + N_A0 + N_Y0 + N_R0
        self.s = s
        self.p = p
        self.a = a
        self.y = y
        self.time = np.zeros((cycles, self.start_steps+1))
        self.N_S = np.zeros((cycles, self.start_steps+1))
        self.N_P = np.zeros((cycles, self.start_steps+1))
        self.N_A = np.zeros((cycles, self.start_steps+1))
        self.N_Y = np.zeros((cycles, self.start_steps+1))
        self.N_R = np.zeros((cycles, self.start_steps+1))
        self.N_S[:, 0] = N_S0
        self.N_P[:, 0] = N_P0
        self.N_A[:, 0] = N_A0
        self.N_Y[:, 0] = N_Y0
        self.N_R[:, 0] = N_R0
        self.cycles = cycles
        self.avg_steps = avg_steps

    def __resize_np_array(self, array):
        pad_size = (int((array.shape[1] - 1) * 2))
        return np.pad(array, ((0, 0), (0, (pad_size))), 'constant', constant_values=0)

    def simulation(self):
        for i in range(self.cycles):
            j = 0
            while self.N_P[i, j] + self.N_Y[i, j] + self.N_A[i, j] != 0:
                if j > self.N_S.shape[1] - 2:
                    self.N_S = self.__resize_np_array(self.N_S)
                    self.N_P = self.__resize_np_array(self.N_P)
                    self.N_A = self.__resize_np_array(self.N_A)
                    self.N_Y = self.__resize_np_array(self.N_Y)
                    self.N_R = self.__resize_np_array(self.N_R)
                    self.time = self.__resize_np_array(self.time)
                S_rate = self.N_S[i, j] / self.total
                P_rate = self.N_P[i, j] / self.total
                A_rate = self.N_A[i, j] / self.total
                Y_rate = self.N_Y[i, j] / self.total
                R_rate = self.N_R[i, j] / self.total
                T = (self.s * S_rate) + (self.p * P_rate) + \
                    (self.a * A_rate) + (self.y * Y_rate) + (R_rate)

                # S = S - 1
                # P = P + 1
                P = self.p_theta * ((self.s * (self.p + self.a + self.y) / T)) * self.beta * (S_rate) * \
                    ((self.p_kappa * P_rate) +
                     (self.a_kappa * A_rate) + (self.y_kappa * Y_rate))
                # S = S - 1
                # A = A + 1
                A = self.a_theta * ((self.s * (self.p + self.a + self.y) / T)) * self.beta * (S_rate) \
                    * ((self.p_kappa * P_rate) + (self.a_kappa * A_rate)
                       + (self.y_kappa * Y_rate))

                # P = P - 1
                # Y = Y + 1
                Y = self.phi * P_rate

                # A = A - 1
                # R = R + 1
                R_1 = self.gamma * A_rate

                # Y = Y - 1
                # R = R + 1
                R_2 = self.zeta * Y_rate

                event_rate = P + A + Y + R_1 + R_2
                # Virus died
                if event_rate == 0:
                    self.N_S[i, j+1:] = self.N_S[i, j]
                    self.N_P[i, j+1:] = self.N_P[i, j]
                    self.N_A[i, j+1:] = self.N_A[i, j]
                    self.N_Y[i, j+1:] = self.N_Y[i, j]
                    self.N_R[i, j+1:] = self.N_R[i, j]
                    break

                u1 = np.random.random()
                tau = 1/event_rate * np.log(1/u1)
                self.time[i, j+1] = self.time[i, j] + tau

                event = np.random.choice(["P", "A", "Y", "R_1", "R_2"], p=[
                    P/event_rate, A/event_rate, Y/event_rate, R_1/event_rate, R_2/event_rate])
                if event == "P":
                    self.N_S[i, j+1] = self.N_S[i, j] - 1
                    self.N_P[i, j+1] = self.N_P[i, j] + 1
                    self.N_A[i, j+1] = self.N_A[i, j]
                    self.N_Y[i, j+1] = self.N_Y[i, j]
                    self.N_R[i, j+1] = self.N_R[i, j]
                elif event == "A":
                    self.N_S[i, j+1] = self.N_S[i, j] - 1
                    self.N_P[i, j+1] = self.N_P[i, j]
                    self.N_A[i, j+1] = self.N_A[i, j] + 1
                    self.N_Y[i, j+1] = self.N_Y[i, j]
                    self.N_R[i, j+1] = self.N_R[i, j]
                elif event == "Y":
                    self.N_S[i, j+1] = self.N_S[i, j]
                    self.N_P[i, j+1] = self.N_P[i, j] - 1
                    self.N_A[i, j+1] = self.N_A[i, j]
                    self.N_Y[i, j+1] = self.N_Y[i, j] + 1
                    self.N_R[i, j+1] = self.N_R[i, j]
                elif event == "R_1":
                    self.N_S[i, j+1] = self.N_S[i, j]
                    self.N_P[i, j+1] = self.N_P[i, j]
                    self.N_A[i, j+1] = self.N_A[i, j] - 1
                    self.N_Y[i, j+1] = self.N_Y[i, j]
                    self.N_R[i, j+1] = self.N_R[i, j] + 1
                else:
                    self.N_S[i, j+1] = self.N_S[i, j]
                    self.N_P[i, j+1] = self.N_P[i, j]
                    self.N_A[i, j+1] = self.N_A[i, j]
                    self.N_Y[i, j+1] = self.N_Y[i, j] - 1
                    self.N_R[i, j+1] = self.N_R[i, j] + 1
                j += 1

    def __clean_up_conf(self, confs, avgs):
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

    def __get_avg_and_conf(self):
        time_max = self.time.max()
        Time_avg = np.linspace(0, time_max, self.avg_steps)
        N_S_avg = np.zeros(self.avg_steps)
        N_P_avg = np.zeros(self.avg_steps)
        N_A_avg = np.zeros(self.avg_steps)
        N_Y_avg = np.zeros(self.avg_steps)
        N_R_avg = np.zeros(self.avg_steps)

        N_S_conf = np.zeros(self.avg_steps, dtype=object)
        N_P_conf = np.zeros(self.avg_steps, dtype=object)
        N_A_conf = np.zeros(self.avg_steps, dtype=object)
        N_Y_conf = np.zeros(self.avg_steps, dtype=object)
        N_R_conf = np.zeros(self.avg_steps, dtype=object)

        N_S_avg[0] = self.N_S0
        N_P_avg[0] = self.N_P0
        N_A_avg[0] = self.N_A0
        N_Y_avg[0] = self.N_Y0
        N_R_avg[0] = self.N_R0

        for i in range(0, self.avg_steps):
            time_max = Time_avg[i]
            S_slice = []
            P_slice = []
            A_slice = []
            Y_slice = []
            R_slice = []
            total_count = 0

            for j in range(self.cycles):
                S_sum = 0
                P_sum = 0
                A_sum = 0
                Y_sum = 0
                R_sum = 0
                step_count = 0

                for k in range(self.N_S.shape[1] - 1):
                    if self.time[j, k] <= time_max and self.time[j, k + 1] > time_max:
                        step_count += 1
                        total_count += 1
                        S_sum += self.N_S[j, k]
                        P_sum += self.N_P[j, k]
                        A_sum += self.N_A[j, k]
                        Y_sum += self.N_Y[j, k]
                        R_sum += self.N_R[j, k]

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

                N_S_conf[i] = stats.t.interval(self.conf_level, df=len(
                    S_slice)-1, loc=np.mean(S_slice), scale=np.std(S_slice, ddof=1) / np.sqrt(len(S_slice)))
                N_P_conf[i] = stats.t.interval(self.conf_level, df=len(
                    P_slice)-1, loc=np.mean(P_slice), scale=np.std(P_slice, ddof=1) / np.sqrt(len(P_slice)))
                N_A_conf[i] = stats.t.interval(self.conf_level, df=len(
                    A_slice)-1, loc=np.mean(A_slice), scale=np.std(A_slice, ddof=1) / np.sqrt(len(A_slice)))
                N_Y_conf[i] = stats.t.interval(self.conf_level, df=len(
                    Y_slice)-1, loc=np.mean(Y_slice), scale=np.std(Y_slice, ddof=1) / np.sqrt(len(Y_slice)))
                N_R_conf[i] = stats.t.interval(self.conf_level, df=len(
                    R_slice)-1, loc=np.mean(R_slice), scale=np.std(R_slice, ddof=1) / np.sqrt(len(R_slice)))

        N_S_lower_conf, N_S_higher_conf = self.__clean_up_conf(
            N_S_conf, N_S_avg)
        N_P_lower_conf, N_P_higher_conf = self.__clean_up_conf(
            N_P_conf, N_P_avg)
        N_A_lower_conf, N_A_higher_conf = self.__clean_up_conf(
            N_A_conf, N_A_avg)
        N_Y_lower_conf, N_Y_higher_conf = self.__clean_up_conf(
            N_Y_conf, N_Y_avg)
        N_R_lower_conf, N_R_higher_conf = self.__clean_up_conf(
            N_R_conf, N_R_avg)

        return (Time_avg, (N_S_avg, N_S_lower_conf, N_S_higher_conf), (N_P_avg, N_P_lower_conf, N_P_higher_conf), (N_A_avg, N_A_lower_conf, N_A_higher_conf), (N_Y_avg, N_Y_lower_conf, N_Y_higher_conf), (N_R_avg, N_R_lower_conf, N_R_higher_conf))

    def __plot_avg_confs(self, axs, title, x, y, low_conf_y, high_conf_y, color):
        axs.set_title(title)
        axs.plot(x, y, marker="", color=color, linewidth=0.5, alpha=0.9)
        axs.plot(x, low_conf_y, marker=".",
                 color="black", linewidth=0.5, alpha=0.5)
        axs.plot(x, high_conf_y, marker=".",
                 color="black", linewidth=0.5, alpha=0.5)
        axs.set_xlim(0, max(x))

    def __visualizer(self, Time_avg, N_S_res, N_P_res, N_A_res, N_Y_res, N_R_res, plot_name):
        plt.ylabel("Population")
        plt.xlabel("time")
        fig, axs = plt.subplots(5, 1, figsize=(6, 20))
        self.__plot_avg_confs(axs[0], "Susceptible", Time_avg,
                              N_S_res[0], N_S_res[1], N_S_res[2], "red")
        self.__plot_avg_confs(axs[1], "Pre-symptomatic", Time_avg,
                              N_P_res[0], N_P_res[1], N_P_res[2], "blue")
        self.__plot_avg_confs(axs[2], "Asymptomatic", Time_avg,
                              N_A_res[0], N_A_res[1], N_A_res[2], "orange")
        self.__plot_avg_confs(axs[3], "Symptomatic", Time_avg,
                              N_Y_res[0], N_Y_res[1], N_Y_res[2], "yellow")
        self.__plot_avg_confs(axs[4], "Recovered", Time_avg,
                              N_R_res[0], N_R_res[1], N_R_res[2], "pink")

        plt.savefig(plot_name)

    def __data_collector(self, Time_avg, N_S_avg, N_P_avg, N_A_avg, N_Y_avg, N_R_avg):
        peak_infections = 0
        peak_time = 0
        for i in range(0, self.avg_steps):
            num_infect = N_P_avg[i] + N_A_avg[i] + N_Y_avg[i]
            if num_infect > peak_infections:
                peak_infections = num_infect
                peak_time = Time_avg[i]

        R0 = ((self.p * self.a * self.y) * self.beta) / \
            ((self.p_theta * (self.zeta + self.phi)) + (self.a_theta * (self.gamma)))

        attack_rate = N_R_avg[len(N_R_avg) - 1] / N_S_avg[0]
        return (("peak_infections", peak_infections), ("peak_time", peak_time), ("R0", R0), ("attack_rate", attack_rate))

    def export_data(self, plot_file_name, metrics_file_name, N_S_name, N_P_name, N_A_name, N_Y_name, N_R_name):
        Time_avg, N_S_res, N_P_res, N_A_res, N_Y_res, N_R_res = self.__get_avg_and_conf()
        self.__visualizer(Time_avg, N_S_res, N_P_res, N_A_res,
                          N_Y_res, N_R_res, plot_file_name)
        metrics = self.__data_collector(
            Time_avg, N_S_res[0], N_P_res[0], N_A_res[0], N_Y_res[0], N_R_res[0])
        metrics_df = pd.DataFrame(metrics)
        metrics_df.to_csv(metrics_file_name)
        N_S_df = pd.DataFrame(self.N_S)
        N_S_df.to_csv(N_S_name)
        N_P_df = pd.DataFrame(self.N_P)
        N_P_df.to_csv(N_P_name)
        N_A_df = pd.DataFrame(self.N_A)
        N_A_df.to_csv(N_A_name)
        N_Y_df = pd.DataFrame(self.N_Y)
        N_Y_df.to_csv(N_Y_name)
        N_R_df = pd.DataFrame(self.N_R)
        N_R_df.to_csv(N_R_name)


def main(seed, N_S0, N_P0, N_A0, N_Y0, N_R0, s, p, a, y, cycles, avg_steps, plot_name, metrics_name, N_S_name, N_P_name, N_A_name, N_Y_name, N_R_name):
    simulation = Simulation(seed, N_S0, N_P0, N_A0, N_Y0, N_R0, s, p, a, y, cycles, avg_steps)
    simulation.simulation()
    simulation.export_data("plot1.png", "metrics.csv",
                           "N_S.csv", "N_P.csv", "N_A.csv", "N_Y.csv", "N_R.csv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="Behavioral Effect on Covid Simulation",
        description="Simulate the effect of behavioral changes on the COVID-19 epidemic.")

    parser.add_argument(
        "seed", type=int, help="value used to initialize random number generator")
    parser.add_argument(
        "N_S0", type=int, help="initial amount of susceptible population")
    parser.add_argument(
        "N_P0", type=int, help="initial amount of pre-symptomatic population")
    parser.add_argument(
        "N_A0", type=int, help="initial amount of asymptomatic population")
    parser.add_argument(
        "N_Y0", type=int, help="initial amount of symptomatic population")
    parser.add_argument(
        "N_R0", type=int, help="initial amount of recoverd population")

    parser.add_argument(
        "s", type=float, help="the fraction of contact that susceptible members will reduce")
    parser.add_argument(
        "p", type=float, help="the fraction of contact that pre-symptomatic members will reduce")
    parser.add_argument(
        "a", type=float, help="the fraction of contact that asymptomatic members will reduce")
    parser.add_argument(
        "y", type=float, help="the fraction of contact that symptomatic members will reduce")

    parser.add_argument("cycles", type=int,
                        help="the number of cycles the simulation will run")
    parser.add_argument("avg_steps", type=int,
                        help="the number of equally distant in time averages we will be computing over the simulation time")

    parser.add_argument("plot_name", type=str, help="plot file name to output")
    parser.add_argument("metrics_name", type=str,
                        help="metrics csv file name to output")
    parser.add_argument("N_S_name", type=str,
                        help="Susceptible array csv file name to output")
    parser.add_argument("N_P_name", type=str,
                        help="Pre-symptomatic array csv file name to output")
    parser.add_argument("N_A_name", type=str,
                        help="Asymptomatic array csv file name to output")
    parser.add_argument("N_Y_name", type=str,
                        help="Symptomatic array csv file name to output")
    parser.add_argument("N_R_name", type=str,
                        help="Recovered array csv file name to output")

    args = parser.parse_args()
    main(args.seed, args.N_S0, args.N_P0, args.N_A0, args.N_Y0, args.N_R0, args.s, args.p, args.a, args.y, args.cycles, args.avg_steps,
         args.plot_name, args.metrics_name, args.N_S_name, args.N_P_name, args.N_A_name, args.N_Y_name, args.N_R_name)
