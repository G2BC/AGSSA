class Model:
    def __init__(self, name="", rep=2, nmbofclasses=0, nmbofclusters=0):

        self.name = name
        self.ListOfVarSites = []  # CBUC
        self.NmbOfClasses = nmbofclasses  # CBUC
        self.CodeOfClass = []  # CBUC
        self.GroundTruth = []  # ANNOTATION
        self.NmbOfClusters = nmbofclusters  # CLOPE
        self.TheClusterOfClass = {}  # CLOPE
        self.TheSpeciesOfCluster = {}  # CLOPE + ANNOTATION
        self.repulsion = rep  # CLOPE
        self.MinClusterResolution = 0  # MODEL
        self.Resolution = 0  # MODEL
        self.DistributionMatrix = []  # MODEL
        self.Accession = []  # DATA - GVF
        self.CLOPE_repulsion = 0

    @staticmethod
    def print_model(model, n):
        print("name:", model.name)
        print("ListOfVarSites:\n", model.ListOfVarSites)
        print("NmbOfClasses:", model.NmbOfClasses)
        print("CodeOfClass[:", n, "]:\n", model.CodeOfClass[:n])
        print("GroundTruth[:", n, "]:\n", model.GroundTruth[:n])
        print("Accession[:", n, "]:\n", model.Accession[:n])
        print("repulsion:", model.repulsion)
        print("NmbOfClusters:", model.NmbOfClusters)
        print("TheSpeciesOfCluster:\n", model.TheSpeciesOfCluster)
        print("MinClusterAcc:", model.MinClusterResolution)
        print("Resolution:", model.Resolution)
        print("CLOPE repulsion parameter:", model.CLOPE_repulsion)
        print("DistributionMatrix:\n", model.DistributionMatrix)
        print("TheClusterOfClass:\n", model.TheClusterOfClass)