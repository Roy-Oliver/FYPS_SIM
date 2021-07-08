from SingleStage import Stage1, Stagei
from HeatSinks import NoHeatSink, RectangularFin, PinFin
from dotenv import load_dotenv
import os

class Distiller:
    """Models a Distiller"""
    def _MultiStageAssemble(self):
        # Assembles the multiple stages

        # Create an empty list for the stages
        self.stages = []

        for stagenum in range(self.N):
            if stagenum == 0:
                # The stage is the first stage
                self.stages.append(Stage1(self.qsun, self.tinf, self.tf, self.t, self.a, self.b, self.k, self.epsilon, self.dwa298, self.c))

            elif stagenum == self.N - 1:
                # The last stage
                self.stages.append(Stagei(0, self.tinf, 0, self.t, self.a, self.b, self.k, self.epsilon, self.dwa298, self.c))

                # qin and Tb is unknown here

            else:
                # Middle stages
                self.stages.append(Stagei(0, self.tinf, 0, self.t, self.a, self.b, self.k, self.epsilon, self.dwa298, self.c))

                # qin and Tb is unknown here

    def _HeatSinkAssemble(self):
        if self.htsnk == 1:
            # No Heat Sink Type

            # Parameters
            a = self.a
            tbn = 0 # Leave blank
            tinf = self.tinf
            c = self.c

            # Create Object
            self.heatsink = NoHeatSink(a, tbn, tinf, c)

        elif self.htsnk == 2:
            # Rectangular Fin

            # Parameters
            a = self.a
            tf = self.param[0]
            l = self.param[1]
            n = self.param[2]
            tbn = 0 # Leave blank
            tinf = self.tinf
            ks = self.param[3]
            c = self.c

            # Create object
            self.heatsink = RectangularFin(a, tf, l, n, tbn, tinf, ks, c)

        else:
            # Pin Fin

            # Parameters
            a = self.a
            d = self.param[0]
            l = self.param[1]
            n = self.param[2]
            tbn = 0 # Leave blank
            tinf = self.tinf
            m = self.param[3]
            z = self.param[4]
            ks = self.param[5]
            c = self.c


            # Create Object
            self.heatsink = PinFin(a, d, l, n, tbn, tinf, ks, c, m, z)

    def __init__(self, qsun, tinf, t, a, b, k, epsilon, dwa298, N, htsnk, param, c):
        self.qsun = qsun
        self.tinf = tinf
        self.tf = 200 # A dummy
        self.t = t
        self.a = a
        self.b = b
        self.k = k
        self.epsilon = epsilon
        self.dwa298 = dwa298
        self.N = N # Number of stages
        self.htsnk = htsnk # Type of Heat Sink
        self.param = param # Heat Sink parameters
        self.c = c

        self._MultiStageAssemble()
        self._HeatSinkAssemble()

    def _solvestages(self):
        # Solves the stages from top to bottom

        # Solve multistage
        for stage_num, stage in enumerate(self.stages):
            if stage_num == len(self.stages) - 1:
                # Solve last stage
                stage.solve()


            else:
                # Solve stages before the last stage
                stage.solve()

                # Transfer qout -> qin and tb -> tf
                self.stages[stage_num + 1].qin = stage.qout
                self.stages[stage_num + 1].tf = stage.tb

        # Solve heatsink
        # Transfer last stage temp to last heatsink
        self.heatsink.tbn = self.stages[-1].tb
        # Solve the heat sink
        self.heatsink.solve()

    def solve(self, tf_guess):

        # Input Guess to Initial Stage
        self.stages[0].tf = tf_guess

        # Solve the Stages
        self._solvestages()

        # Return error
        err = self.stages[-1].qout - (self.heatsink.qs / (self.a ** 2))
        return err



