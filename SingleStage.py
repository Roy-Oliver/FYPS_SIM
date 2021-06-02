from CoolProp.CoolProp import PropsSI, HAPropsSI
from dotenv import load_dotenv
import os
import math


class Stagei:
    """A class that models the ith stage of the desalinator"""

    def __init__(self, qin, tinf, tb, t, a, b, k, epsilon, dwa298):
        # Initializes the stage object

        self.qin = qin
        self.tinf = tinf
        self.tb = tb
        self.t = t
        self.a = a
        self.b = b
        self.k = k
        self.epsilon = epsilon
        self.dwa298 = dwa298

    def _psi(self):
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

        try:
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
        except:
            print("Warning: Tbar < Tinf")
            self.tw = 0
            self.qside = (self.tbar - self.tinf) / (1 / 2.487 + self.t / self.k) # haside = 2.487


    def _Jevap(self):

        # Get P
        load_dotenv()
        p = float(os.environ.get("P"))

        try:
            # get hlf and hvf
            hlf = PropsSI("HMOLAR", "T", self.tinf, "P", p, "Water")
            self.hvf = PropsSI("HMOLAR", "T", self.tf, "Q", 1, "Water")

            # Return jevap
            self.jevap = (self.qin - self.qcond) / (self.hvf - hlf)  # mol/sm2
        except ValueError:
            print("Warning: Value Error for hlf and/or hvf. Switching to default heat of vaporization")
            self.jevap = (self.qin - self.qcond) / 43988

    def _qcond(self):

        # Solve for pf
        try:
            pf = PropsSI("P", "T", self.tf, "Q", 1, "Water")
        except ValueError:
            print("Warning: Value Error for pf. Switching to Antoine Equation")
            pf = (10 ** (4.6543 - 1435.264 / (self.tf - 64.848))) * 100000 # Vapor pressure in Pascal

        # Solve for pb
        try:
            pb = PropsSI("P", "T", self.tb, "Q", 1, "Water")
        except ValueError:
            print("Warning: Value Error for pb. Switching to Antoine Equation")
            pb = (10 ** (4.6543 - 1435.264 / (self.tb - 64.848))) * 100000  # Vapor pressure in Pascal

        # Solve for pbar
        self.pbar = (pf + pb) / 2

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        # Get P
        load_dotenv()
        p = float(os.environ.get("P"))

        # Solve for ka
        try:
            ka = HAPropsSI("K", "T", self.tbar, "P", p, "P_w", self.pbar)
        except ValueError:
            print("Warning: Value Error for ka. Switching to ka = 0.026")
            ka = 0.026

        # Solve for qcond
        self.qcond = ka * (self.tf - self.tb) / self.b

    def _tf(self):

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        try:
            # Solve for cb
            cb = PropsSI("DMOLAR", "T", self.tb, "Q", 1, "Water")  # mol/m3
        except ValueError:
            print("Warning: Value Error for cb. Switching to Antoine Equation and Ideal Gas Equation")
            # Solve for vapor pressure
            pb = (10 ** (4.6543 - 1435.264 / (self.tb - 64.848))) * 100000 # pascals
            cb = pb / (8.314 * self.tb) # ideal gas law

        # Solve for dwa
        dwa = self.dwa298 * ((self.tbar / 298.15) ** 1.75)

        # solve for cf:
        cf = self.jevap * self.b / dwa + cb  # mol/m3

        # Solve for Tf
        try:
            self.tf = PropsSI("T", "DMOLAR", cf, "Q", 1, "Water")
        except ValueError:
            print("Warning: Value Error for tf. Switching to Antoine Equation and Ideal Gas Law")
            # Get tolerance
            load_dotenv()
            delta = float(os.environ.get("delta"))
            R = float(os.environ.get("R"))

            # Iteratively solve for Tf based on Ideal gas law and Antoine Equation
            self.tf = self.tb # Guess an initial value of tf
            while True:
                tf_1 = self.tf # store an initial value of tf

                # COmpute for vapor pressure via antoine equation
                pf = (10 ** (4.6543 - 1435.264 / (self.tf - 64.848))) * 100000

                # Compute for Tf via ideal gas equation
                self.tf = pf / (R * cf)

                if abs(self.tf - tf_1) > delta:
                    continue
                else:
                    break

    def _jside(self):

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        try:
            # Solve for hvside and hlb
            hvside = PropsSI("HMOLAR", "T", self.tbar, "Q", 1, "Water")  # J/mol
            self.hlb = PropsSI("HMOLAR", "T", self.tb, "Q", 0, "Water")  # J/mol

            self.jside = self.qside / (hvside - self.hlb)  # mol/sm2

        except ValueError:
            print("Warning: Value Error for hvside and/hlb. Switching to default heat of vaporization")
            # Solve for jside
            self.jside = self.qside / 43988

    def _qout(self):

        try:
            self.qout = ((self.a ** 2) * (self.qcond + self.hvf * self.jevap) - 4 * self.a * self.b * self.qside - self.hlb * ((self.a ** 2) * self.jevap - 4 * self.a * self.b * self.jside)) / (self.a ** 2)
        except AttributeError:
            print("Warning: hvf and/or hlb were not found. Switching to default heat of vaporization")
            # If there was an attribute error in previous code ie. self.hvf or self.hlb not defined,
            # Set hlb = 0 and
            # Define hvf such that its value is hlb + enthalpy of vaporization
            self.hlb = 0
            self.hvf = self.hlb + 43988
            self.qout = ((self.a ** 2) * (self.qcond + self.hvf * self.jevap) - 4 * self.a * self.b * self.qside - self.hlb * ((self.a ** 2) * self.jevap - 4 * self.a * self.b * self.jside)) / (self.a ** 2)

    def _ntot(self):
        self.ntot = ((self.a ** 2) * (self.qout - self.qcond) + 4 * self.a * self.b * self.qside) / (
                    (self.a ** 2) * self.qin)

    def solve(self):

        # Guess Tf
        self.tf = (self.tb - self.tinf) + self.tb

        # Get delta
        load_dotenv()
        delta = float(os.environ.get("delta"))

        while True:
            # Initial value of tb
            tf_1 = self.tf

            # Update values
            self._qcond()
            self._Jevap()
            self._tf()

            if abs(tf_1 - self.tf) > delta:
                continue
            else:
                break

        # Update Values
        self._psi()
        self._jside()
        self._qout()
        self._ntot()




class Stage1(Stagei):
    """A class that models the first stage of the desalinator"""

    def __init__(self, qsun, tinf, tb, t, a, b, k, epsilon, dwa298):
        """Initialize attributes of parent class"""
        super().__init__(qsun, tinf, tb, t, a, b, k, epsilon, dwa298)

    def _qrad(self):
        # Computes for q''rad

        # Get sigma
        load_dotenv()
        sigma = float(os.environ.get("sigma"))

        # Output qrad
        self.qrad = self.epsilon * sigma * ((self.tf ** 4) - (self.tinf ** 4))

    def _Jevap(self):
        # Get P
        load_dotenv()
        p = float(os.environ.get("P"))

        try:
            # get hlf and hvf
            hlf = PropsSI("HMOLAR", "T", self.tinf, "P", p, "Water")
            self.hvf = PropsSI("HMOLAR", "T", self.tf, "Q", 1, "Water")

            # Return jevap
            self.jevap = (self.qin - self.qrad - self.qcond) / (self.hvf - hlf)  # mol/sm2

        except ValueError:
            print("Warning: hvf and/or hlf were not found. Switching to default heat of vaporization")
            self.jevap = (self.qin - self.qrad - self.qcond) / 43988

    def solve(self):

        # Guess Tf
        self.tf = (self.tb - self.tinf) + self.tb

        # Get delta
        load_dotenv()
        delta = float(os.environ.get("delta"))

        while True:
            # Initial value of tb
            tf_1 = self.tf

            # Update values
            self._qrad()
            self._qcond()
            self._Jevap()
            self._tf()

            if abs(tf_1 - self.tf) > delta:
                continue
            else:
                break

        # Update Values
        self._psi()
        self._jside()
        self._qout()
        self._ntot()