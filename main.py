from Distiller import Distiller
import scipy.optimize
from dotenv import load_dotenv
import os

# Input Distiller Variables

# q''sun in W/m2
qsun = 1000

# Tinf in K
tinf = 298.15

# Insulation thickness in m
t = 0.011

# Stage length in m
a = 0.12

# Stage thickness in m
b = 0.00215

# Insulation conductivity in W/mk
k = 0.25

# Solar absorber emissivity
epsilon = 0.03

# Diffusivity of water in air at 298.15 K, m2/s
dwa298 = 0.000026

# Number of stages
N = 1

# Type of heat sink. 1 for no heat sink, 2 for rectangular, 3 for pin
htsnk = 2

# Parameters for "no heat sink". [--blank--]
param1 = []

# Parameters for "rectangular fin", [tf -> fin thickness(m), L -> fin length (m), n -> number of fins, ks -> conductivity of material (W/mK)]
param2 = [0.01, 0.01, 3, 0.05]

# Parameters for "pin fin", [d -> fin diameter (m), L -> fin length (m), n -> number of fins, , ks -> conductivity of material]
param3 = [0.01, 0.1, 4, 0.05]

if __name__ == "__main__":

    if htsnk == 1:
        desalinator_setup = Distiller(qsun, tinf, t, a, b, k, epsilon, dwa298, N, htsnk, param1)
    elif htsnk == 2:
        desalinator_setup = Distiller(qsun, tinf, t, a, b, k, epsilon, dwa298, N, htsnk, param2)
    else:
        desalinator_setup = Distiller(qsun, tinf, t, a, b, k, epsilon, dwa298, N, htsnk, param2)



    # Import delta
    load_dotenv()
    delta = float(os.environ.get("delta"))

    # Solve the setup
    scipy.optimize.newton(desalinator_setup.solve, 350, rtol=delta)

    # Print results
    print("_________________________________________")


    # Loop through results, and sum up qevap from each stage.
    # The sum of the qevap are used to compute for the total efficiency
    qevap_tot = 0
    for stage_num, stage in enumerate(desalinator_setup.stages):

        qevap_tot += stage.qevap

        print(f"Stage {stage_num + 1}\n")

        print(f"Tf = {stage.tf}")
        print(f"Tb = {stage.tb}\n")

        print(f"q''in = {stage.qin}")
        print(f"q''out = {stage.qout}")
        print(f"q''cond = {stage.qcond}")
        print(f"q''side = {stage.qside}")
        print(f"q''evap = {stage.qevap}")
        if stage_num == 0:
            print(f"q''rad = {stage.qrad}")

        print(f"\nJ''evap = {stage.jevap}\n")

        print("_________________________________________")

    print(f"ntot = {qevap_tot / qsun}")