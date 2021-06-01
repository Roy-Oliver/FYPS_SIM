from CoolProp.CoolProp import PropsSI
from dotenv import load_dotenv
import os
import math


class NoHeatSink:
    """Models a Flat Plate (No Heat Sink)"""

    def __init__(self, a, tbn, tinf):
        self.a = a
        self.tbn = tbn
        self.tinf = tinf

    def solve(self):
        # Solve for Tfilms
        tfilms = (self.tbn + self.tinf) / 2

        # Get P and g
        load_dotenv()
        p = int(os.environ.get("P"))
        g = float(os.environ.get("g"))

        # Obtain properties of air
        myus = PropsSI("V", "T", tfilms, "P", p, "Air")
        rhos = PropsSI("DMASS", "T", tfilms, "P", p, "Air")
        cps = PropsSI("CPMASS", "T", tfilms, "P", p, "Air")
        kas = PropsSI("CONDUCTIVITY", "T", tfilms, "P", p, "Air")
        betas = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", tfilms, "P", p, "Air")

        # Solve for Grs
        grs = (((self.a ** 2) / (4 * self.a)) ** 3) * (rhos ** 2) * g * betas * (self.tbn - self.tinf) / (myus ** 2)

        # Solve for Prs
        prs = (cps * myus) / kas

        # Solve for NumThins
        numthins = 0.527 * ((grs * prs) ** (1 / 5)) / ((1 + ((1.9 / prs) ** (9 / 10))) ** (2 / 9))

        # Solve for has
        has = kas * 2.5 / ((self.a ** 2) * math.log(1 + 2.5 / numthins) / (4 * self.a))

        # Solve for qs
        self.qs = has * (self.a ** 2) * (self.tbn - self.tinf)

class RectangularFin:
    """Models a Multiple Rectangular Fin Array Heat Sink"""

    def __init__(self, a, tf, l, n, tbn, tinf, ks):
        self.a = a
        self.tf = tf  # thickness of fin
        self.l = l
        self.n = n
        self.tbn = tbn
        self.tinf = tinf
        self.ks = ks

    def solve(self):
        # Solve for Tfilms
        tfilms = (self.tbn + self.tinf) / 2

        # Get P and g
        load_dotenv()
        p = int(os.environ.get("P"))
        g = float(os.environ.get("g"))

        # Solve for properties of air
        myus = PropsSI("V", "T", tfilms, "P", p, "Air")
        rhos = PropsSI("DMASS", "T", tfilms, "P", p, "Air")
        cps = PropsSI("CPMASS", "T", tfilms, "P", p, "Air")
        kas = PropsSI("CONDUCTIVITY", "T", tfilms, "P", p, "Air")
        betas = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", tfilms, "P", p, "Air")

        # Solve for z
        z = self.a / self.n - self.tf

        # Solve for El
        el = g * betas * (self.tbn - self.tinf) * (z ** 4) / ((myus / rhos) * (kas * self.l / (rhos * cps)))

        # Solve for hf
        hf = (kas / z) * (((576 / (el ** 2)) + 2.873 / (el ** (1/2))) ** (-1/2))

        # Solve for prs
        prs = cps * myus / kas

        # Solve for grs
        grs = (((self.a / self.n - self.tf) * self.a / (2 * (self.a / self.n - self.tf) + 2 * self.a)) ** 3) * (rhos ** 2) * g * betas * (self.tbn - self.tinf) / (myus ** 2)

        # Solve for numthins
        numthins = 0.527 * ((grs * prs) ** (1/5)) / ((1 + ((1.9 / prs) ** (9/10))) ** (2/9))

        # Solve for hbs
        antf = (self.a / self.n - self.tf) # Helper variable for (a / n -tf)
        hbs = kas * 2.5 / ((antf * self.a * math.log(1 + 2.5 / numthins)) / (2 * antf + 2 * self.a))

        # Solve for qs
        fac1 = hf * self.ks * self.tf # Helper variable for hf * ks * tf
        fac2 = hf * self.l / (self.ks * self.tf) # Helper variable for hf * l / ks * tf
        self.qs = self.n * self.a * (((2 * fac1) ** (1/2)) * (self.tbn - self.tinf) * math.tanh(2 * fac2) + hbs * (self.tbn - self.tinf) * z)

class PinFin:
    """Models a multiple rectangular pin fin heat sink"""

    def __init__(self, a, d, l, n, tbn, tinf, ks):
        self.a = a
        self.d = d
        self.l = l
        self.n = n
        self.tbn = tbn
        self.tinf = tinf
        self.ks = ks

    def _hbs(self):
        # Solve for Tfilms
        tfilms = (self.tbn + self.tinf) / 2

        # Get P and g
        load_dotenv()
        p = int(os.environ.get("P"))
        g = float(os.environ.get("g"))

        # Solve for properties of air
        myus = PropsSI("V", "T", tfilms, "P", p, "Air")
        rhos = PropsSI("DMASS", "T", tfilms, "P", p, "Air")
        cps = PropsSI("CPMASS", "T", tfilms, "P", p, "Air")
        kas = PropsSI("CONDUCTIVITY", "T", tfilms, "P", p, "Air")
        betas = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", tfilms, "P", p, "Air")

        # Solve for prs
        prs = cps * myus / kas

        # Solve for grs1
        fac1 = ((self.a / self.n - self.d) * self.a / (2 * (self.a / self.n - self.d) + 2 * self.a)) ** 3  # Helper variable
        grs1 = fac1 * (rhos ** 2) * g * betas * (self.tbn - self.tinf) / (myus ** 2)

        # Solve for numthins1
        fac2 = (1 + (1.9 / prs) ** (9/10)) ** (2/9)  # Helper variable
        numthins1 = 0.527 * ((grs1 * prs) ** (1/5)) / fac2

        # Solve for hbs
        fac3 = self.a / self.n - self.d  # Helper variable
        fac4 = 2 * fac3 + 2 * self.a
        hbs = kas * 2.5 / (fac3 * math.log(1 + 2.5 / numthins1) / fac4)
        return hbs

    def _hf(self):
        # Solve for Tfilms
        tfilms = (self.tbn + self.tinf) / 2

        # Get P and g
        load_dotenv()
        p = int(os.environ.get("P"))
        g = float(os.environ.get("g"))

        # Solve for properties of air
        myus = PropsSI("V", "T", tfilms, "P", p, "Air")
        rhos = PropsSI("DMASS", "T", tfilms, "P", p, "Air")
        cps = PropsSI("CPMASS", "T", tfilms, "P", p, "Air")
        kas = PropsSI("CONDUCTIVITY", "T", tfilms, "P", p, "Air")
        betas = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", tfilms, "P", p, "Air")

        # Solve for prs
        prs = cps * myus / kas

        # Solve for grs2
        grs2 = (self.l ** 3) * (rhos ** 2) * g * betas * (self.tbn - self.tinf) / (myus ** 2)

        # Solve for cl
        fac1 = (1 + (0.492 / prs) ** (9/16)) ** (4/9)
        cl = 0.671 / fac1

        # Solve for numthins2
        numthins2 = cl * ((grs2 * prs) ** (1/4))

        # Solve for zeta
        zeta = (1.8 * self.l / self.d) / numthins2

        # Solve for numlams
        numlams = zeta * 2 / (math.log(1 + zeta) * math.log(1 + 2 / numthins2))

        # Solve for numturbs
        fac2 = (1 + 0.61 * (prs ** 0.22)) ** 0.42
        fac3 = 1 + 1400000000 / grs2
        numturbs = 0.13 * (prs ** 0.22) * ((grs2 * prs) ** (1/3)) / (fac2 * fac3)

        # Solve for hf
        hf = (kas / self.l) * (((numlams ** 6) + (numturbs ** 6)) ** (1/6))
        return hf

    def solve(self):
        # Get pi
        load_dotenv()
        pi = float(os.environ.get("pi"))

        # Solve for hf and hbs
        hf = self._hf()
        hbs = self._hbs()

        # Solve for qs
        fac1 = math.tanh(((4 * hf / (self.ks * self.d)) ** (1/2)) * self.l)
        fac2 = (4 * hf / (self.ks * self.d)) ** (1/2)
        fac3 = (self.a ** 2) - (self.n ** 2) * pi * (self.d ** 2) / 4
        self.qs = (self.n ** 2) * (fac1 / fac2) * (pi * self.d * hf) * (self.tbn - self.tinf) + hbs * (self.tbn - self.tinf) * fac3
