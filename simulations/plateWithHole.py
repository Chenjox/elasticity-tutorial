import math

import os 
import inspect
from pathlib import Path
from os.path import abspath


class ParameterSet(dict):
    def __init__(self, *args, **kwargs):
        super(ParameterSet, self).__init__(*args, **kwargs)
        self.__dict__ = self

parameterSet = ParameterSet()

#############################################
#  File Names
#############################################
# Get path of executing module
modulePath = Path(inspect.getfile(inspect.currentframe()))
# Grab directory
thisModuleDir = modulePath.parent
# ASSUMPTION: Every File and Dir is relative to this file!

# Grid File Name
gridFileName = "plateWithHole.msh"
gridFileNamePath = Path(thisModuleDir, Path(gridFileName) )

#############################################
#  Grid parameters
#############################################

# Path specification
#relative = Path("diskWithHole.msh")
#absolute = relative.absolute()  # absolute is a Path object

parameterSet.gridFile = abspath(gridFileNamePath.absolute()) #"diskWithHole.msh" #abspath(absolute) # "plateWithHole.msh"
parameterSet.refinementLevels = 0

# Print
print(parameterSet.gridFile)
#print(self.__file__)

#############################################

youngsModulus = 210000.0 # N/mm^2
poissonsRatio = 0.3

parameterSet.firstLameParameter = poissonsRatio/((1.0 - 2.0*poissonsRatio)*(1.0 + poissonsRatio))*youngsModulus # N/mm^2
parameterSet.secondLameParameter = 0.5*1.0/(1.0 + poissonsRatio) * youngsModulus # N/mm^2

#############################################

# A indicator function for dirichlet values. 
# It gets a position and should return a list of BOOL
def DirichletIndicatorFunction(x):
    return [x[0] < 1e-10, x[1] < 1e-10]

class DirichletBoundary:
    def __init__(self, homotopyParameter):
        self.homotopyParameter = homotopyParameter

    # Deformation of 3d classical materials
    def dirichletValues(self, x):
        # Clamp the L-shape in its reference configuration
        return [0.0, 0.0]

class NeumannBoundary:
    def __init__(self, homotopyParameter):
        self.homotopyParameter = homotopyParameter

    # Current Position and Normal of the boundary.
    def neumannValues(self, x):
        return [x[0],x[1]]

def bodyforce(x):
    return [0.0, -9.81]