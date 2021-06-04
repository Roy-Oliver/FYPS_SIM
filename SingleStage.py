from CoolProp.CoolProp import PropsSI, HAPropsSI
from dotenv import load_dotenv
import os
import math


class Stagei:
    """A class that models the ith stage of the desalinator"""

    def __init__(self, qin, tinf, tf, t, a, b, k, epsilon, dwa298):
        # Initializes the stage object

        self.qin = qin
        self.tinf = tinf
        self.tf = tf
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

        # Get P and hvap
        load_dotenv()
        p = float(os.environ.get("P"))
        hvap = float(os.environ.get("hvap"))

        try:
            # get hlf and hvf
            self.hlf = PropsSI("HMOLAR", "T", self.tinf, "P", p, "Water")
            self.hvf = PropsSI("HMOLAR", "T", self.tf, "Q", 1, "Water")

            # Return jevap
            self.jevap = (self.qin - self.qcond) / (self.hvf - self.hlf)  # mol/sm2
        except ValueError:
            print("Warning: Value Error for hlf and/or hvf. Switching to default heat of vaporization")
            self.jevap = (self.qin - self.qcond) / hvap

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
            print("Warning: Value Error for ka. Switching to ka = 0.027")
            ka = 0.027

        # Solve for qcond
        self.qcond = ka * (self.tf - self.tb) / self.b

    def _tb(self):

        # Solve for tbar
        self.tbar = (self.tf + self.tb) / 2

        try:
            # Solve for cf
            cf = PropsSI("DMOLAR", "T", self.tf, "Q", 1, "Water")  # mol/m3
        except ValueError:
            print("Warning: Value Error for cb. Switching to Antoine Equation and Ideal Gas Equation")
            # Solve for vapor pressure
            pf = (10 ** (4.6543 - 1435.264 / (self.tf - 64.848))) * 100000 # pascals
            cf = pf / (8.314 * self.tf) # ideal gas law

        # Solve for dwa
        dwa = self.dwa298 * ((self.tbar / 298.15) ** 1.75)

        # solve for cb:
        cb = -(self.jevap * self.b / dwa - cf)  # mol/m3


        # Solve for Tb
        try:
            self.tb = PropsSI("T", "DMOLAR", cb, "Q", 1, "Water")
        except ValueError:
            print("Warning: Value Error for tb. Switching to Antoine Equation and Ideal Gas Law")
            # Get tolerance
            load_dotenv()
            delta = float(os.environ.get("delta"))
            R = float(os.environ.get("R"))

            # The upper value for solution to exist for tf is 59259 mol/m3.
            if cb > 59259:
                print("No solution for Tb. Exiting")
                exit()

            # Iteratively solve for Tb based on Ideal gas law and Antoine Equation
            self.tb = self.tf # Guess an initial value of tb

            while True:

                tb_1 = self.tb # store an initial value of tb

                # COmpute for back wall vapor pressure via ideal gas equation
                pb = (cb * R * self.tb) / (10 ** 5) # vapor pressure in bar

                # Compute for new back wall temperature using antoine equation
                self.tb = 1435.264 / (4.6543 - math.log(pb, 10)) + 64.848

                # check for convergence
                if abs(self.tb - tb_1) > delta:
                    continue
                else:
                    break


    def _qout(self):

        # Load hvap
        load_dotenv()
        hvap = float(os.environ.get("hvap"))

        try:
            self.hlb = PropsSI("HMOLAR", "T", self.tb, "Q", 0, "Water")  # J/mol
            self.qout = ((self.a ** 2) * (self.qcond + self.hvf * self.jevap - self.hlb * self.jevap) - 4 * self.a * self.b * self.qside) / (self.a ** 2)
        except AttributeError:
            print("Warning: hvf and/or hlb were not found. Switching to default heat of vaporization")
            # If there was an attribute error in previous code ie. self.hvf or self.hlb not defined,
            # Set hlb = 0 and
            # Define hvf such that its value is hlb + enthalpy of vaporization
            self.hlb = 0
            self.hvf = self.hlb + hvap
            self.qout = ((self.a ** 2) * (self.qcond + self.hvf * self.jevap - self.hlb * self.jevap) - 4 * self.a * self.b * self.qside) / (self.a ** 2)

    def _qevap(self):

        # Solves for qevap using energy balance
        self.qevap = ((self.a ** 2) * (self.qout - self.qcond) + 4 * self.a * self.b * self.qside) / (self.a ** 2)

    def solve(self):

        # Guess Tb
        self.tb = self.tinf

        # Get delta
        load_dotenv()
        delta = float(os.environ.get("delta"))

        while True:
            # Initial value of tb
            tb_1 = self.tb

            # Update values
            self._qcond()
            self._Jevap()
            self._tb()

            if abs(tb_1 - self.tb) > delta:
                continue
            else:
                break

        # Update Values
        self._psi()
        self._qout()
        self._qevap()


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
        # Get P and hvap
        load_dotenv()
        p = float(os.environ.get("P"))
        hvap = float(os.environ.get("hvap"))

        try:
            # get hlf and hvf
            self.hlf = PropsSI("HMOLAR", "T", self.tinf, "P", p, "Water")
            self.hvf = PropsSI("HMOLAR", "T", self.tf, "Q", 1, "Water")

            # Return jevap
            self.jevap = (self.qin - self.qrad - self.qcond) / (self.hvf - self.hlf)  # mol/sm2

        except ValueError:
            print("Warning: hvf and/or hlf were not found. Switching to default heat of vaporization")
            self.jevap = (self.qin - self.qrad - self.qcond) / hvap

    def solve(self):

        # Guess Tf
        self.tb = self.tinf

        # Get delta
        load_dotenv()
        delta = float(os.environ.get("delta"))

        while True:
            # Initial value of tb
            tb_1 = self.tb

            # Update values
            self._qrad()
            self._qcond()
            self._Jevap()
            self._tb()

            if abs(tb_1 - self.tb) > delta:
                continue
            else:
                break

        # Update Values
        self._psi()
        self._qout()
        self._qevap()