import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
import argparse


class Simulation:
    # Constant variables
    # Real world values for the following variables are taken from their above citation:

    # Xueting Qui. “The role of asymptomatic and pre-symptomatic infection in
    # SARS-CoV-2 transmission—a living systematic review”. In: Clin Microbial
    # Infect (Jan. 2021). doi: 10.1016/j.cmi.2021.01.011.

    p_kappa = 0.07  # pre-symptomatic rate of infection on contact
    a_kappa = 0.01  # asymptomatic rate of infection on contact
    y_kappa = 0.06  # symptomatic rate of infection on contact

    # ======================
    # Diana Buitrago-Garcia et al. “Occurrence and transmission potential of
    # asymptomatic and presymptomatic SARS-CoV-2 infections: Update of a
    # living systematic review and meta-analysis”. In: PLOS Medicine 19 (Apr.
    # 2022). doi: 10.1371/journal.pmed.1003987.

    p_theta = 0.81  # the percent of infections that become pre-symptomatic
    a_theta = 0.19  # the percent of infections that become asymptomatic

    # ======================
    # Ayesha S. Mahmud Dennis M. Feehan. “Quantifying population contact
    # patterns in the United States during the COVID-19 pandemic”. In: Nat
    # Commun 12 (Sept. 2021). doi: 10.1038/s41467-021-20990-2.

    beta = 3.333333  # population general rate of contact

    # ======================
    # Yue Xiang et al. “COVID-19 epidemic prediction and the impact of pub-
    # lic health interventions: A review of COVID-19 epidemic models”. In: In-
    # fectious Disease Modelling 6 (2021), pp. 324–342. issn: 2468-0427. doi:
    # https://doi.org/10.1016/j.idm.2021.01.001. url: https://www.
    # sciencedirect.com/science/article/pii/S2468042721000038.

    gamma = 0.2  # rate of asymptomatic to removed

    # ======================
    # Hualei Xin et al. “Estimating the Latent Period of Coronavirus Disease
    # 2019 (COVID-19)”. In: Clinical Infectious Diseases 74.9 (Sept. 2021), pp. 1678–
    # 1681. issn: 1058-4838. doi: 10 . 1093 / cid / ciab746. eprint: https : / /
    # academic.oup.com/cid/article-pdf/74/9/1678/43525252/ciab746.
    # pdf. url: https://doi.org/10.1093/cid/ciab746.

    phi = 0.7  # rate of pre-symptomatic to symptomatic

    # ======================
    # Yue Xiang et al. “COVID-19 epidemic prediction and the impact of pub-
    # lic health interventions: A review of COVID-19 epidemic models”. In: In-
    # fectious Disease Modelling 6 (2021), pp. 324–342. issn: 2468-0427. doi:
    # https://doi.org/10.1016/j.idm.2021.01.001. url: https://www.
    # sciencedirect.com/science/article/pii/S2468042721000038.

    zeta = 0.28  # rate of symptomatic to removed

    # ======================

    start_steps = 10  # the initial size of the array
    conf_level = 0.95  # confidence level of confidence interval

    def __init__(self, seed, N_S0, N_P0, N_A0, N_Y0, N_R0, s, p, a, y, cycles, avg_steps):
        """Initialize variables for simulation and calculating metrics."""
        np.random.seed = seed
        # initial values for each compartment
        self.N_S0 = N_S0
        self.N_P0 = N_P0
        self.N_A0 = N_A0
        self.N_Y0 = N_Y0
        self.N_R0 = N_R0

        # total population
        self.total = N_S0 + N_P0 + N_A0 + N_Y0 + N_R0

        # the fraction of contact that each compartment will reduce during a pandemic
        self.s = s
        self.p = p
        self.a = a
        self.y = y

        # create two dimensional array for the number of cycles and length start_steps for each compartment
        self.time = np.zeros((cycles, self.start_steps+1))
        self.N_S = np.zeros((cycles, self.start_steps+1))
        self.N_P = np.zeros((cycles, self.start_steps+1))
        self.N_A = np.zeros((cycles, self.start_steps+1))
        self.N_Y = np.zeros((cycles, self.start_steps+1))
        self.N_R = np.zeros((cycles, self.start_steps+1))

        # initalize each cycle to the initial value of each compartment
        self.N_S[:, 0] = N_S0
        self.N_P[:, 0] = N_P0
        self.N_A[:, 0] = N_A0
        self.N_Y[:, 0] = N_Y0
        self.N_R[:, 0] = N_R0

        self.cycles = cycles  # the number of cycles to to repeat the simulation
        # the number of equally distant steps in time to average each compartment
        self.avg_steps = avg_steps

    def __resize_np_array(self, array):
        """Resize array to double its length inheriting its current values and initializing indexes to 0."""
        pad_size = (int((array.shape[1] - 1) * 2))
        return np.pad(array, ((0, 0), (0, (pad_size))), 'constant', constant_values=0)

    def simulation(self):
        """Simulate the epidemic."""
        # repeat the simulation for each cycle
        for i in range(self.cycles):
            # store each arrays current index
            j = 0

            # run simulation until infectious compartments have reached 0
            while self.N_P[i, j] + self.N_Y[i, j] + self.N_A[i, j] != 0:
                # check if arrays have run out of space
                if j > self.N_S.shape[1] - 2:
                    # resize arrays
                    self.N_S = self.__resize_np_array(self.N_S)
                    self.N_P = self.__resize_np_array(self.N_P)
                    self.N_A = self.__resize_np_array(self.N_A)
                    self.N_Y = self.__resize_np_array(self.N_Y)
                    self.N_R = self.__resize_np_array(self.N_R)
                    self.time = self.__resize_np_array(self.time)

                # Compute the current rate of each compartment
                S_rate = self.N_S[i, j] / self.total
                P_rate = self.N_P[i, j] / self.total
                A_rate = self.N_A[i, j] / self.total
                Y_rate = self.N_Y[i, j] / self.total
                R_rate = self.N_R[i, j] / self.total

                # compute the total rate of contact of components
                T = (self.s * S_rate) + (self.p * P_rate) + \
                    (self.a * A_rate) + (self.y * Y_rate) + (R_rate)

                # calculate event rate:
                # S = S - 1
                # P = P + 1
                # self.p_theta
                # -   rate of infections that are pre-symptomatic
                # ((self.s * (self.p + self.a + self.y) / T))
                # -   fraction of contacts made by susceptible members with infectious members
                # self.beta
                # -   general rate of contact
                # (S_rate)
                # -   current rate of susceptible population
                # ((self.p_kappa * P_rate) + (self.a_kappa * A_rate)
                #       + (self.y_kappa * Y_rate))
                # -   current rate of infectious population
                P = self.p_theta * (self.s * (self.p + self.a + self.y) / T) * self.beta * (S_rate) * \
                    ((self.p_kappa * P_rate) +
                     (self.a_kappa * A_rate) + (self.y_kappa * Y_rate))

                # calculate event rate:
                # S = S - 1
                # A = A + 1
                # self.a_theta
                # -    rate of infections that are asymptomatic
                # ((self.s * (self.p + self.a + self.y) / T))
                # -    fraction of contacts made by susceptible members with infectious members
                # self.beta
                # -    general rate of contact
                # (S_rate)
                # -    current rate of susceptible population
                # ((self.p_kappa * P_rate) + (self.a_kappa * A_rate)
                #       + (self.y_kappa * Y_rate))
                # -    current rate of infectious population
                A = self.a_theta * ((self.s * (self.p + self.a + self.y) / T)) * self.beta * (S_rate) \
                    * ((self.p_kappa * P_rate) + (self.a_kappa * A_rate)
                       + (self.y_kappa * Y_rate))

                # calculate event rate:
                # P = P - 1
                # Y = Y + 1
                # self.phi
                # -    rate of pre-symptomatic to symptomatic
                # P_rate
                # -    current rate of pre-symptomatic population
                Y = self.phi * P_rate

                # calculate event rate:
                # A = A - 1
                # R = R + 1
                # self.gamma
                # -    rate of asymptomatic to removed
                # A_rate
                # -   current rate of asymptomatic population
                R_1 = self.gamma * A_rate

                # calculate event probability:
                # Y = Y - 1
                # R = R + 1
                # self.zeta
                # -    rate of symptomatic to removed
                # Y_rate
                # -    current rate of symptomatic
                R_2 = self.zeta * Y_rate

                # total event rate
                event_rate = P + A + Y + R_1 + R_2

                # if virus is dead replace rest of the array with last event value
                # infectious arrays do not need to be modified as no events means
                # infectious population is already zero
                if event_rate == 0:
                    self.N_S[i, j+1:] = self.N_S[i, j]
                    self.N_R[i, j+1:] = self.N_R[i, j]
                    self.time[i, j+1:] = self.time[i, j]
                    break

                # compute time to next event
                u1 = np.random.random()
                tau = 1/event_rate * np.log(1/u1)
                self.time[i, j+1] = self.time[i, j] + tau

                # randomly choose the event based on the proportion of each event to the total of all events
                event = np.random.choice(["P", "A", "Y", "R_1", "R_2"], p=[
                    P/event_rate, A/event_rate, Y/event_rate, R_1/event_rate, R_2/event_rate])
                # update arrays for each event
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

                # update array index
                j += 1

    def __clean_up_conf(self, confs, avgs):
        """Clean up NaNs in confidence intervals by replacing each NaN with its avg."""
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
        """Calculate averages and confidences of results."""
        # Create avg arrays from 0 to max_time with length avg_steps
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

        for i in range(0, self.avg_steps):
            time_max = Time_avg[i]
            # slices to store each cycles value for populations for current avg_step
            S_slice = []
            P_slice = []
            A_slice = []
            Y_slice = []
            R_slice = []
            # count number of values added for each current avg_step
            total_count = 0

            for j in range(self.cycles):
                length = self.N_S.shape[1]
                for k in range(length - 1):
                    # store all populations at time_max
                    if self.time[j, k] <= time_max and self.time[j, k + 1] > time_max:
                        total_count += 1
                        S_slice.append(self.N_S[j, k])
                        P_slice.append(self.N_P[j, k])
                        A_slice.append(self.N_A[j, k])
                        Y_slice.append(self.N_Y[j, k])
                        R_slice.append(self.N_R[j, k])

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
        """Plot the average and confidence interval of component."""
        axs.set_title(title)
        axs.plot(x, y, marker="", color=color, linewidth=0.5, alpha=0.9)
        axs.plot(x, low_conf_y, marker=".",
                 color="black", linewidth=0.5, alpha=0.5)
        axs.plot(x, high_conf_y, marker=".",
                 color="black", linewidth=0.5, alpha=0.5)
        axs.set_xlim(0, max(x))

    def __visualizer(self, Time_avg, N_S_res, N_P_res, N_A_res, N_Y_res, N_R_res, plot_name):
        """Create plot."""
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
        """Calculate metrics."""
        # Calculate peak infections and peak time
        peak_infections = 0
        peak_time = 0
        for i in range(0, self.avg_steps):
            num_infect = N_P_avg[i] + N_A_avg[i] + N_Y_avg[i]
            if num_infect > peak_infections:
                peak_infections = num_infect
                peak_time = Time_avg[i]

        # R0 calculation is based on:
        # Fred Brauer. “A simple model for behaviour change in epidemics”. In: BMC
        # Public Health 11.S1 (Feb. 2011). doi: 10.1186/1471-2458-11-s1-s3.
        R0 = ((self.p * self.a * self.y) * self.beta) / \
            ((self.p_theta * (self.zeta + self.phi)) + (self.a_theta * (self.gamma)))

        # Calculate attack rate
        attack_rate = N_R_avg[len(N_R_avg) - 1] / N_S_avg[0]

        return (("peak_infections", peak_infections), ("peak_time", peak_time), ("R0", R0), ("attack_rate", attack_rate))

    def __export_to_csv(self, array, name):
        """Export csv of array as name."""
        array_df = pd.DataFrame(array)
        array_df.to_csv(name)

    def export_data(self, plot_file_name, metrics_file_name, N_S_name, N_P_name, N_A_name, N_Y_name, N_R_name, time_name):
        """Export plot and metrics."""

        # Create plot and compute metrics
        Time_avg, N_S_res, N_P_res, N_A_res, N_Y_res, N_R_res = self.__get_avg_and_conf()
        self.__visualizer(Time_avg, N_S_res, N_P_res, N_A_res,
                          N_Y_res, N_R_res, plot_file_name)
        metrics = self.__data_collector(
            Time_avg, N_S_res[0], N_P_res[0], N_A_res[0], N_Y_res[0], N_R_res[0])

        self.__export_to_csv(metrics, metrics_file_name)
        self.__export_to_csv(self.N_S, N_S_name)
        self.__export_to_csv(self.N_P, N_P_name)
        self.__export_to_csv(self.N_A, N_A_name)
        self.__export_to_csv(self.N_Y, N_Y_name)
        self.__export_to_csv(self.N_R, N_R_name)
        self.__export_to_csv(self.time, time_name)


def main(seed, N_S0, N_P0, N_A0, N_Y0, N_R0, s, p, a, y, cycles, avg_steps, plot_name, metrics_name, N_S_name, N_P_name, N_A_name, N_Y_name, N_R_name, time_name):
    simulation = Simulation(seed, N_S0, N_P0, N_A0, N_Y0,
                            N_R0, s, p, a, y, cycles, avg_steps)
    simulation.simulation()
    simulation.export_data(plot_name, metrics_name,
                           N_S_name, N_P_name, N_A_name, N_Y_name, N_R_name, time_name)


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
    parser.add_argument("time_name", type=str,
                        help="Time array csv file name to output")

    args = parser.parse_args()
    main(args.seed, args.N_S0, args.N_P0, args.N_A0, args.N_Y0, args.N_R0, args.s, args.p, args.a, args.y, args.cycles, args.avg_steps,
         args.plot_name, args.metrics_name, args.N_S_name, args.N_P_name, args.N_A_name, args.N_Y_name, args.N_R_name, args.time_name)
