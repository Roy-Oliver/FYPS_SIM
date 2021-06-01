from SingleStage import Stagei, Stage1

class MultiStage:
    """Models a Multi-Stage Distiller without the heat sink"""

    def _assemble(self):
        # Assembles the multiple stages

        # Create an empty list for the stages
        self.stages = []

        for stagenum in range(self.N):
            if stagenum == 0:
                # The stage is the first stage
                self.stages.append(Stage1(self.qsun, self.tinf, 0, self.t, self.a, self.b, self.k, self.epsilon, self.dwa298))
                # Tb of the first stage is unknown here

            elif stagenum == self.N - 1:
                # The last stage
                self.stages.append(Stagei(0, self.tinf, self.tb, self.t, self.a, self.b, self.k, self.epsilon, self.dwa298))

                # qin is unknown here

            else:
                # Middle stages
                self.stages.append(Stagei(0, self.tinf, 0, self.t, self.a, self.b, self.k, self.epsilon, self.dwa298))

                # qin and Tb is unknown here

    def __init__(self, qsun, tinf, tb, t, a, b, k, epsilon, dwa298, N):

        self.qsun = qsun
        self.tinf = tinf
        self.tb = tb
        self.t = t
        self.a = a
        self.b = b
        self.k = k
        self.epsilon = epsilon
        self.dwa298 = dwa298
        self.N = N # Number of stages

        self._assemble()

