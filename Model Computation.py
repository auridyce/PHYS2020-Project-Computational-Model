###PRESSURE-BASED ROCKET EQUATION-DERIVED COMPUTATIONAL MODEL###
###FOR WATER BOTTLE ROCKET PROJECT IN THERMAL AND STATISTICAL PHYSICS###
###BY U7470329 CASSIDY O'BRIEN###

#IMPORTS
import numpy as np #Numpy arrays for easier handling
from uncertainties import ufloat #Uncertainty value class
from uncertainties.umath import * #Values with uncertainties
#import matplotlib as plot #Graphical plots for video
#import pandas #Output into csvs, comparison against models


##CONSTANTS##
#Physical Constants - Fundamental regardless of environment, no measurements
gasConst = ufloat(8.314, 0) #J(k*mol)^-1
waterDensity = ufloat(10**5, 0) #kgm^-3
gravityLocal = ufloat(9.796, 0) #nkg^-1

#Experimental Constants - Unchanging with given experimental setup, measured
#ejectionRatio = ufloat(30, 8) #Dimensionless
bottleArea = ufloat(6.08*(10**-3),5*(10**-4)) #m^2
nozzleArea = ufloat(2.01*(10**-4),5*(10**-5)) #m^2
airMolPerPump = ufloat(0.1,0.01) #mol/pumps
temperature = ufloat(283,2) #kelvin
bottleVolume = ufloat(1.5*(10**-3),0) #m^3
bottleMass = ufloat(0.032, 0.01) #kg
airDensity = ufloat(1.23, 0.1) #kgm^-3
dragCoefficient = ufloat(0.8, 0.05) #Dimensionless, estimated from water bottle shape. (Mostly blunt cylinder)

#Notes on other constants:
#Many constants (I.E air pressure) are obscured in other ratios or cancelled
#See scans of text in PDF for detailed derivations. Video covers in simple detail.

##MODEL EQUATIONS##
#General Model Equations - Used for both instant and variable models

#WaterMass - Calculates mass of water in bottle from total mass of bottle
def waterMass(totalLaunchMass):
    return (totalLaunchMass - bottleMass)

#Differential pressure at launch - Calculates the differential pressure induced by the bike pump
def dPressureLaunch(pumps, waterMassLaunch):
    numerator = (airMolPerPump*pumps*gasConst*temperature)
    denominator = (bottleVolume- waterMassLaunch/waterDensity)
    return (numerator/denominator)

#Exhaust velocity - Calculates the exhaust velocity of the water from the rocket
def exhaustVelocity(pumps, waterMassLaunch):
    numerator = dPressureLaunch(pumps, waterMassLaunch)*nozzleArea
    return (numerator/waterMassLaunch)

#Air Resistance Term, general calculation
def resistanceTerm(velocity, area):
    return 0.5*(velocity**2)*area*dragCoefficient

#Instant Model Equations - Assume instant action of launch, and water release
def equivalentVelocityInstantModel(pumps, totalLaunchMass):
    waterMassLaunch = waterMass(totalLaunchMass)
    numerator = nozzleArea*dPressureLaunch(pumps, waterMassLaunch)+exhaustVelocity(pumps, waterMassLaunch)
    #numerator = nozzleArea*dPressureLaunch(pumps, waterMassLaunch)
    massEject = waterMassLaunch
    return (numerator/massEject)

def velocityInstantFinal(pumps, totalLaunchMass):
    return ((equivalentVelocityInstantModel(pumps, totalLaunchMass) * log(totalLaunchMass/bottleMass))-gravityLocal)

#Instant with air resistance


def instantVelocityAirResist(pumps, totalLaunchMass):
    return (velocityInstantFinal(pumps, totalLaunchMass)-(resistanceTerm(velocityInstantFinal(pumps,totalLaunchMass),bottleArea)/bottleMass))

#Computational time-based model - Assumes initial pressure and water transfers are instant
#But water release and acceleration is computed over time
timeStep = 1/60

def instantaneousMassEjection(pumps, waterMassLaunch):
    return (exhaustVelocity(pumps, waterMassLaunch) * nozzleArea * waterDensity)

def variableMass(pumps, totalLaunchMass, t):
    waterMassLaunch = waterMass(totalLaunchMass)
    massFactor = t*instantenousMassEjection(pumps, waterMassLaunch)
    if massFactor <= waterMass:
        return (totalLaunchMass - massFactor)
    else:
        return (totalLaunchMass - waterMass)



