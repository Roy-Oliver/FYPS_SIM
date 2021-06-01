from CoolProp.CoolProp import PropsSI, HAPropsSI
from dotenv import load_dotenv
import os
import math


class Stage1:
    """A class that models the first stage of the desalinator"""

    def __init__(self, qsun, tinf, tb, t, a, b, k, epsilon, dwa298):
        # Initializes the stage object

        self.qsun = qsun
        self.tinf = tinf
        self.tb = tb
        self.t = t
        self.a = a
        self.b = b
        self.k = k
        self.epsilon = epsilon
        self.dwa298 = dwa298

    def __psi(self):
        # Solves for q''side

        # Solves for tbar
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
            numturb = (0.13 * (pr ** 0.22) / ((1 + 0.61 * (pr ** 0.22)) ** 0.42)) * (
                        ((gr * pr) ** (1 / 3)) / (1 + 1400000000 / gr))

            # Update cl
            cl = 0.671 / ((1 + (0.492 / pr) ** (9 / 16)) ** (4 / 9))

            # Update Numthin
            numthin = cl * ((gr * pr) ** (1 / 4))

            # Update Numlam
            numlam = 2.0 / math.log(1 + 2 / numthin)

            # Update haside
            haside = (kaside / self.b) * (((numlam ** 6) + (numturb ** 6)) ** (1 / 6))

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
        # Computes for q''rad

        # Get sigma
        load_dotenv()
        sigma = float(os.environ.get("sigma"))

        # Output qrad
        self.qrad = self.epsilon * sigma * ((self.tf ** 4) - (self.tinf ** 4))

    def __Jevap(self):

        # Get P
        load_dotenv()
        p = float(os.environ.get("P"))

        # get hlf and hvf
        hlf = PropsSI("HMOLAR", "T", self.tinf, "P", p, "Water")
        self.hvf = PropsSI("HMOLAR", "T", self.tf, "Q", 1, "Water")

        # Return jevap
        self.jevap = (self.qsun - self.qrad - self.qcond) / (self.hvf - hlf)  # mol/sm2

    def __qcond(self):

        # Solve for pbar
        pf = PropsSI("P", "T", self.tf, "Q", 1, "Water")
        pb = PropsSI("P", "T", self.tb, "Q", 1, "Water")
        self.pbar = (pf + pb) / 2

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        # Get P
        load_dotenv()
        p = float(os.environ.get("P"))

        # Solve for ka
        ka = HAPropsSI("K", "T", self.tbar, "P", p, "P_w", self.pbar)

        # Solve for qcond
        self.qcond = ka * (self.tf - self.tb) / self.b

    def __tf(self):

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        # Solve for cb
        cb = PropsSI("DMOLAR", "T", self.tb, "Q", 1, "Water")  # mol/m3

        # Solve for dwa
        dwa = self.dwa298 * ((self.tbar / 298.15) ** 1.75)

        # Return Tf
        cf = self.jevap * self.b / dwa + cb  # mol/m3
        self.tf = PropsSI("T", "DMOLAR", cf, "Q", 1, "Water")

    def __jside(self):

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        # Solve for hvside and hlb
        hvside = PropsSI("HMOLAR", "T", self.tbar, "Q", 1, "Water")  # J/mol
        self.hlb = PropsSI("HMOLAR", "T", self.tb, "Q", 0, "Water")  # J/mol

        self.jside = self.qside / (hvside - self.hlb)  # mol/sm2

    def __qout(self):

        self.qout = ((self.a ** 2) * (
                    self.qcond + self.hvf * self.jevap) - 4 * self.a * self.b * self.qside - self.hlb * (
                                 (self.a ** 2) * self.jevap - 4 * self.a * self.b * self.jside)) / (self.a ** 2)

    def __ntot(self):
        self.ntot = ((self.a ** 2) * (self.qout - self.qcond) + 4 * self.a * self.b * self.qside) / (
                    (self.a ** 2) * self.qsun)

    def phi(self):

        # Guess Tf
        self.tf = (self.tb - self.tinf) + self.tb

        # Get delta
        load_dotenv()
        delta = float(os.environ.get("delta"))

        while True:
            # Initial value of tb
            tf_1 = self.tf

            # Update values
            self.__qrad()
            self.__qcond()
            self.__Jevap()
            self.__tf()

            if abs(tf_1 - self.tf) > delta:
                continue
            else:
                break

        # Update Values
        self.__psi()
        self.__jside()
        self.__qout()
        self.__ntot()
