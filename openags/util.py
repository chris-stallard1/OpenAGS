from abc import ABC, abstractmethod
import math
import itertools

import numpy as np
from numpy import NaN, inf, arange, asarray
import xylib
from sigfig import round

class KnownPeak:
    def __init__(self, elementName="None", center=0.0, unit="mg", sensitivity=None, mass=None, divisor=None, output=None, halfLife=None, decayConstant = None, decayUnit="min"):
        """Represents a peak generated by isotope "elementName" with a predetermined energy centroid and possbly other parameters.
        
        The peak can have a known sensitivity or mass attached to it. This will be the divisor which the area is divided by at the end of the analysis.
        If and only if the analysis is NAA, the peak will also have a half-life or decay constant, representing that of the isotope that is producing it.
        These are used for time correction in the get_tcf() function.
        """
        self.elementName = elementName
        self.center = center
        self.delayed = True

        if divisor and output:
            self.divisor = divisor
            self.output = output
        else:
            if sensitivity != None and mass != None:
                raise TypeError("Please provide no more than 1 of the following: mass, sensitivity")
            if mass != None:
                self.divisor = mass
                self.output = "Sensitivity ("+unit+")"
            elif sensitivity != None:
                self.divisor = sensitivity
                self.output = "Mass ("+unit+")"
            else:
                self.divisor = None
                self.output = "Peak Area (cps)"
        
        if halfLife != None and decayConstant != None:
            raise TypeError("Please provide no more than 1 of the following: half-life, decay constant")
        elif halfLife != None:
            self.decayConstant = math.log(2) / halfLife
        elif decayConstant != None:
            self.decayConstant = decayConstant
        else:
            self.delayed = False

        if self.delayed:
            #these seem reversed because decay constant unit is 1/time
            if decayUnit in ("s", "sec"):
                self.decayConstant *= 60
            if decayUnit in ("h", "hr"):
                self.decayConstant /= 60
    
    #Input Functions: load data into object
    @staticmethod
    def load_from_dict(d):
        """Takes in a dictionary represenation of this object (generated by export_to_dict()), and sets parameters accordingly."""
        elementName = d["ele"]
        center = d["ctr"]

        if "divisor" in d:
            divisor = d["divisor"]
            output = d["output"]
        else:
            divisor = None
            output = None
        
        if "decayConstant" in d:
            decayConstant = d["decayConstant"]
        else:
            decayConstant = None
        return KnownPeak(elementName, center, divisor=divisor, output=output, decayConstant=decayConstant)

    #Interesting functions: Do more than get or set but aren't inputs/outputs

    def get_tcf(self):
        """Returns the time correction factor to divide by: 1 for PGAA, and dependent on irradiation, wait, and count time for NAA"""
        if self.delayed:
            return (1 - math.e ** (-1 * self.decayConstant * self.irrTime)) * math.e ** (-1 * self.decayConstant * self.waitTime) * (1 - math.e ** (-1 * self.decayConstant * self.countTime))
        else:
            return 1
    

    #Getters and setters

    def set_delay_times(self, irr, wait, count):
        self.irrTime = irr
        self.waitTime = wait
        self.countTime = count

    def get_ctr(self):
        return self.center

    def get_ele(self):
        return self.elementName

    def get_output(self):
        return self.output
    
    def set_NAA_params(self, halfLife=None, decayConstant=None, unit="min"):
        """Set certain parameters for NAA. 
        
        This mirrors some functionality in __init__ but helps to avoid a certain 
        tangled mess of if cases that used to be in parsers.py
        """
        if halfLife != None and decayConstant != None:
            raise TypeError("Please provide no more than 1 of the following: half-life, decay constant")
        elif halfLife != None:
            self.decayConstant = math.log(2) / halfLife
        elif decayConstant != None:
            self.decayConstant = decayConstant
        else:
            raise TypeError("Either half-life or decay constant should be provided")

        if unit in ("s", "sec"):
            self.decayConstant *= 60
        elif unit in ("h", "hr"):
            self.decayConstant /= 60

    #Output Methods: Transform the data in the object into something

    def export_to_dict(self):
        """Exports all information about this peak to a dictionary, from which the peak can later be loaded."""
        outDict = {
                "ele" : self.elementName,
                "ctr" : self.center
                }
        
        if self.divisor != None:
            outDict["divisor"] = self.divisor
            outDict["output"] = self.output

        if self.delayed:
            outDict["decayConstant"] = self.decayConstant
        
        return outDict

    def __str__(self):
        """Returns a string representation of this peak as its element name concatenated with its pek centroid."""
        return self.elementName + " : " + str(round(float(self.center), decimals=1))
    
    def get_results(self, area, areaStdev):
        """Combines this peak (with known sensitivity or mass value) with a peak fround in the data to get the mass/sensitivity of that peak.

        The peak in the data is fitted and found to have area "area" and area std. dev. "areaStdev". It is then matched to this known peak,
        and this method returns the desired output for the original peak given that match.
        """
        if self.divisor == None:
            return [area/self.get_tcf(), areaStdev/self.get_tcf()]
        return [area/self.divisor/self.get_tcf(), areaStdev/self.divisor/self.get_tcf()]

def multiple_peak_and_background(peaks, background, x, params):
    """Function used for curve_fit, splits the parameters between the component peaks and background.
    
    Uses get_num_params() function of peaks/backgrounds in order to determine how many paramaters to allocate them.
    Computes y values for a given set of parameters and returns them.
    """
    y = np.zeros_like(x)
    leftCounter = 0
    rightCounter = background.get_num_params()
    y += background.get_ydata_with_params(x,params[leftCounter:rightCounter])
    
    for peak in peaks:
        leftCounter = rightCounter
        rightCounter = leftCounter + peak.get_num_params()
        y += peak.get_ydata_with_params(x,params[leftCounter:rightCounter])
        
    return y

def set_all_params(peaks, background, params, variances, reanalyze):
    """Sets parameters of background/peaks to the result of the fitter

    If the data is being fit for the first time (not in a batch), original params and variances are set. 
    These are the parameters stored in the state file, so that the reloaded fits match the file the user adjusts the fits on.
    If this step wasn't taken, the program would apply the fit from the last file to the data from the first file.
    """
    leftCounter = 0
    rightCounter = background.get_num_params()

    noVars = (variances[0] == inf) or (variances[0] == -1 * inf) #sometimes, variances cannot be determined

    background.set_params(params[leftCounter:rightCounter])
    if not noVars:
        background.set_variances(variances[leftCounter:rightCounter])
    if not reanalyze:
        background.set_original_params(params[leftCounter:rightCounter])
        if not noVars:
            background.set_original_variances(variances[leftCounter:rightCounter])
    for peak in peaks:
        leftCounter = rightCounter
        rightCounter = leftCounter + peak.get_num_params()
        peak.set_params(params[leftCounter:rightCounter])
        if not noVars:
            peak.set_variances(variances[leftCounter:rightCounter])
        if not reanalyze:
            peak.set_original_params(params[leftCounter:rightCounter])
            if not noVars:
                peak.set_original_variances(variances[leftCounter:rightCounter])

def get_curve(peaks, background, x):
    """Sums peaks and background counts for a given set of X values"""
    y = np.zeros_like(x)
    y += background.get_ydata(x)
    for peak in peaks:
        y += peak.get_ydata(x)
    return list(y)

def binary_search_find_nearest(l, e):
    """Find the index of the element in list l which is closest to element e."""
    upperBound = len(l)
    lowerBound = 0
    guess = (upperBound + lowerBound)//2

    #edge cases, where len(l) <= 2
    if guess == 0:
        return guess
    if guess == 1:
        if abs(e - l[0]) < abs(e - l[1]):
            return 0
        return 1

    #Binary Search    
    while not (e < l[guess+1] and e > l[guess-1]):
        if e > l[guess]:
            lowerBound = guess + 1
        else:
            upperBound = guess - 1
        guess = (upperBound + lowerBound)//2
        if guess <= 2 or guess >= len(l)-2:
            break
    if e > l[guess]:
        guess += 1
    return guess

def ivw_combine(meas, stdev = None, variance = None):
    """Combines measurements through Inverse Variance Weighting.

    Arguments: meas -- The measurements to combine
    Keyword Arguments: stdev -- Standard deviations of the measurements
                       variance -- variances of the measurements
    (only stdev XOR variance should be provided)
    """
    var = None
    res = None
    if variance != None:
        var = 1/sum([1/v for v in variance])
        res = var * sum([m/v for m,v in zip(meas, variance)])
    else:
        var = 1/sum([1/s**2 for s in stdev])
        res = var * sum([m/s**2 for m,s in zip(meas, stdev)])
    return [res, var**.5]
