from Distiller import Distiller
import scipy.optimize
from dotenv import load_dotenv
import os
import csv
import numpy as np

# Input Distiller Variables

# q''sun in W/m2
qsun = 1000

# Tinf in K
tinf = 298.15

# Insulation thickness in m
t = 0.01

# Stage length in m
a = 0.148

# Stage width in m
c = 0.115

# Stage thickness in m
b = 0.00215

# Insulation conductivity in W/mk
k = 0.25

# Solar absorber emissivity
epsilon = 0.03

# Diffusivity of water in air at 298.15 K, m2/s
dwa298 = 0.000026

# Type of heat sink. 1 for no heat sink, 2 for rectangular, 3 for pin
htsnk = 3

# Stages to Simulate
N_start = 1
N_end = 10

# Guess for front wall temperature in Kelvin
T_guess = 340

# Parameters for "rectangular fin", [tf -> fin thickness(m), L -> fin length (m), n -> number of fins, ks -> conductivity of material (W/mK)]


# Parameters for "pin fin", [d -> fin diameter (m), L -> fin length (m), n -> number of fins on side a, m -> number of fins on side c, z -> fin spacing (m), ks -> conductivity of material]


def simulate(param, N):

    # Build the desalinator Setup
    if htsnk == 1:
        param1 = param
        desalinator_setup = Distiller(qsun, tinf, t, a, b, k, epsilon, dwa298, N, htsnk, param1, c)
    elif htsnk == 2:
        param2 = param
        desalinator_setup = Distiller(qsun, tinf, t, a, b, k, epsilon, dwa298, N, htsnk, param2, c)
    elif htsnk == 3:
        param3 = param
        desalinator_setup = Distiller(qsun, tinf, t, a, b, k, epsilon, dwa298, N, htsnk, param3, c)
    else:
        print("Invalid heat sink")
        exit()


    # Import delta
    load_dotenv()
    delta = float(os.environ.get("delta"))

    # Solve the setup
    scipy.optimize.newton(desalinator_setup.solve, T_guess, rtol=delta)

    # Loop through results, and sum up qevap from each stage.
    # The sum of the qevap are used to compute for the total efficiency
    qevap_tot = 0
    for stage_num, stage in enumerate(desalinator_setup.stages):
        qevap_tot += stage.qevap
    return [desalinator_setup.stages[-1].qout, qevap_tot / qsun]

def main():
    for N in range(N_start, N_end + 1):
        print(f"Simulating N = {N} stages")

        results = []

        # Open parameters file
        with open("parameters.csv") as csv_file:
            csv_reader = csv.reader(csv_file)

            # Loop Through parameters
            current_n = 1
            for row in csv_reader:

                # Convert data to float
                param = list(map(float, row))

                if htsnk == 2 and param[2] != current_n:
                    results.append(["",""])
                    current_n = param[2]
                if htsnk == 3 and (param[2] ** 2 + param[3] != current_n):
                    results.append(["",""])
                    current_n = param[2] ** 2 + param[3]

                # Simulate
                try:
                    results.append(simulate(param, N))
                except:
                    print(f"error")
                    results.append(["err", "err"])
                    continue

        np.savetxt(f"results(N={N}).csv", results, delimiter=",", fmt='%s')


main()
