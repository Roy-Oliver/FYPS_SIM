from CoolProp.CoolProp import PropsSI, HAPropsSI
from dotenv import load_dotenv
import os
import math


class stage1:
    '''A class that models the first stage of the desalinator'''
    def __init__(self, qsun, tinf, tf, t, tm, a, b, k, km, epsilon, epsilonm, dwa298, r):
        # Initializes the stage object

        self.qsun = qsun
        self.tinf = tinf
        self.tf = tf
        self.t = t
        self.tm = tm
        self.a = a
        self.b = b
        self.k = k
        self.km = km
        self.epsilon = epsilon
        self.epsilonm = epsilonm
        self.dwa298 = dwa298
        self.r = r

    def __psi(self):
        # Solves for q''side

        #Solves for tbar
        self.tbar = (self.tf + self.tb) / 2

        # Guess tw as tbar
        self.tw = self.tbar

        # Get P, g, delta
        load_dotenv()
        p = int(os.environ.get("P"))
        g = float(os.environ.get("g"))
        delta = float(os.environ.get("delta"))

        while True:

            # Compute for tfilmside
            tfilmside = (self.tw + self.tinf) / 2

            # Compute for air properties
            myu = PropsSI("V", "T", tfilmside, "P", p, "Air")
            rho = PropsSI("DMASS", "T", tfilmside, "P", p, "Air")
            cp = PropsSI("CPMASS", "T", tfilmside, "P", p, "Air")
            kaside = PropsSI("CONDUCTIVITY", "T", tfilmside, "P", p, "Air")
            beta = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", tfilmside, "P", p, "Air")

            # Update Pr
            pr = cp * myu / kaside

            # Update Gr
            gr = (self.b ** 3) * (rho ** 2) * g * beta * (self.tw - self.tinf) / (myu ** 2)

            # Update Numturb
            numturb = (0.13 * (pr ** 0.22) / ((1 + 0.61 * (pr ** 0.22)) ** 0.42)) * (((gr * pr) ** (1/3))/(1 + 1400000000 / gr))

            # Update cl
            cl = 0.671 / ((1 + (0.492 / pr) ** (9/16)) ** (4/9))

            # Update Numthin
            numthin = cl * pr * ((gr * pr) ** (1/4))

            # Update Numlam
            numlam = 2.0 / math.log(1 + 2 / numthin)

            # Update haside
            haside = (kaside / self.b) * (((numlam ** 6) + (numturb ** 6)) ** (1/6))

            # Update q ''side
            qside = (self.tbar - self.tinf) / (1 / haside + self.t / self.k)

            # Update tw
            tw_2 = -(qside * (self.t / self.k) - self.tbar)

            # Decision
            if abs(tw_2 - self.tw) > delta:
                self.tw = tw_2
                continue
            else:
                break

        # output qside
        self.qside = qside

    def __qrad(self):
        # Computes for qrad

        # Get sigma
        load_dotenv()
        sigma = float(os.environ.get("sigma"))

        # Output qrad
        self.qrad = self.epsilon * sigma * ((self.tf ** 4) - (self.tinf ** 4))

    def __Jevap(self):

        # Solve for pbar

        pf = PropsSI("P", "T", self.tf, "Q", 0, "Water")
        pb = PropsSI("P", "T", self.tb, "Q", 0, "Water")
        self.pbar = (pf + pb) / 2

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        # Get P, MH2O and R
        load_dotenv()
        mh2o = float(os.environ.get("MH2O"))
        R = float(os.environ.get("R"))
        p = float(os.environ.get("P"))
        pi = float(os.environ.get("pi"))

        # Solve for dwa
        dwa = self.dwa298 * ((self.tbar / 298.15) ** 1.75)

        # Compute for first addend in 1/K
        first_addend = (p - self.pbar) * ((2 - self.epsilonm) ** 2) * R * self.tbar * self.tm / ((self.epsilonm ** 2) * p * dwa * mh2o)

        # Compute for second addend in 1/K
        second_addend = 3 * R * self.tbar * self.tm * ((2 - self.epsilonm) ** 2) / (2 * mh2o * (self.epsilonm ** 2) * self.r * ((8 * R * self.tbar / (pi * mh2o)) ** (1/2)))

        # Compute for third addent in 1 / K
        third_addend = (p - self.pbar) * R * self.tbar * (self.b - self.tm) / (p * dwa * mh2o)

        # Compute for K
        K = (first_addend + second_addend + third_addend) ** (-1)

        # Return J''evap
        self.jevap = K * (pf - pb) # In kg/sm2

    def __qcond(self):
        # Get P
        load_dotenv()
        p = float(os.environ.get("P"))

        # get hlf and hvf
        hlf = PropsSI("HMASS", "T", self.tinf, "P", p, "Water")
        self.hvf = PropsSI("HMASS", "T", self.tf, "Q", 1, "Water")

        #Return qcond
        self.qcond = self.qsun + hlf * self.jevap - self.qrad - self.hvf * self.jevap

    def __tb(self):

        # Solve for pbar

        pf = PropsSI("P", "T", self.tf, "Q", 0, "Water")
        pb = PropsSI("P", "T", self.tb, "Q", 0, "Water")
        self.pbar = (pf + pb) / 2

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        # Get P
        load_dotenv()
        p = float(os.environ.get("P"))

        # Solve for ka
        ka = HAPropsSI("k", "T", self.tbar, "P_w", self.pbar, "P", p)

        # Solve for tb
        self.tb = -(self.qcond * ((self.b - self.tm) / ka + self.tm/self.km) - self.tf)

    def __jside(self):
        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        # Solve for hvside and hlb
        hvside = PropsSI("HMASS", "T", self.tbar, "Q", 1, "Water")
        self.hlb = PropsSI("HMASS", "T", self.tb, "Q", 0, "Water")

        self.jside = self.qside / (hvside - self.hlb)

    def __qout(self):

        self.qout = ((self.a ** 2) * (self.qcond + self.hvf * self.jevap) - 4 * self.a * self.b * self.qside - self.hlb * ((self.a ** 2) * self.jevap - 4 * self.a * self.b * self.jside)) / (self.a ** 2)

    def __ntot(self):
        self.ntot = ((self.a ** 2) * (self.qout - self.qcond) + 4 * self.a * self.b * self.qside) / ((self.a ** 2) * self.qsun)

    def phi(self):

        # Guess Tb
        self.tb = (self.tinf + self.tf) / 2

        # Get delta
        load_dotenv()
        delta = float(os.environ.get("delta"))

        while True:

            # Initial value of tb
            tb_1 = self.tb

            # Update values
            self.__qrad()
            self.__Jevap()
            self.__qcond()
            self.__tb()

            if abs(tb_1 - self.tb) > delta:
                continue
            else:
                break

        self.__psi()
        self.__jside()
        self.__qout()
        self.__ntot()




