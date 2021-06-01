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
        self.tf = tf # thickness of fin
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
        hf = (kas / z) * ((((576 / (el ** 2)) + 2.873 / (el ** 1/2))) ** (-1/2))

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






