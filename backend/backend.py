import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
import itertools
from util import multiple_peak_and_background, get_curve, binary_search_find_nearest, set_all_params, KnownPeak
from constants import default_prefs, som
from parsers import SpectrumParser, StandardsFileParser, CSVWriter, ExcelWriter
from models import GaussianPeak, LinearBackground
from sigfig import round
import os

class ActivationAnalysis:
    """An application class representing the whole backend. 
    
    This is the ONLY backend class that the frontend should interact with.

    ----Variables----
    userPrefs : User Preferences
    fileData : Extracted contents and metadata from project files
    fileList : List of spectrum filenames, data from each correlates to fileData
    knownPeaks : dictionary mapping peak location to KnownPeak onjects
    ROIs : Regions of Interest, see ROI class below
    isotopes: Isotopes being analyzed
    title : Project Title
    ROIsFitted: Whether all regions of interest for this project have been fitted
    resultsGenerated : Whether evaluators have been run and results generated for this project
    delayed : Whether or not this is a Delayed Gamma Analysis
    """
    def __init__(self, userPrefs = default_prefs, title = ""):
        self.userPrefs = userPrefs
        self.fileData = []
        self.fileList = []
        self.knownPeaks = {}
        self.ROIs = []
        self.isotopes = []
        self.title = title
        self.ROIsFitted = False
        self.resultsGenerated = False
        self.delayed = False

    def load_from_dict(self, stored_data):
        """Sets variables for an analysis object based on a dictionary exported by the export_to_dict() function."""

        if "userPrefs" in stored_data.keys(): #otherwise keep default prefs
            self.userPrefs = stored_data["userPrefs"]
        
        self.title = stored_data["title"]
        self.add_files(stored_data["files"])
        self.load_known_peaks(stored_data["standardsFilename"])
        energies = self.fileData[0]["energies"]
        cps = self.fileData[0]["cps"]
        self.ROIsFitted = stored_data["ROIsFitted"]
        for ROIData in stored_data["ROIs"]:
            lowerIndex = ROIData["indicies"][0]
            upperIndex = ROIData["indicies"][1]
            r = ROI(energies[lowerIndex:upperIndex], cps[lowerIndex:upperIndex], [lowerIndex, upperIndex])
            r.load_from_dict(ROIData)
            self.ROIs.append(r)
        self.isotopes = list(set(itertools.chain(*[r.get_isotopes() for r in self.ROIs])))
        self.resultsGenerated = stored_data["resultsGenerated"]
        if self.resultsGenerated:
            for i in range(len(stored_data["results"])):
                self.fileData[i]["results"] = stored_data["results"][i]
                self.fileData[i]["resultHeadings"] = stored_data["resultHeadings"]
                self.fileData[i]["evaluatorNames"] = stored_data["evaluatorNames"]
        self.delayed = stored_data["delayed"]
        if self.delayed and "NAATimes" in stored_data.keys():
            for i in range(len(stored_data["NAATimes"])):
                self.fileData[i]["NAATimes"] = stored_data["NAATimes"][i]


    def export_to_dict(self):
        """Exports the current project state to a dictionary"""
        exportROIs = [r.export_to_dict() for r in self.ROIs]
        outDict =  {
            "userPrefs" : self.userPrefs,
            "title" : self.title,
            "files" : self.fileList,
            "standardsFilename" : self.standardsFilename,
            "ROIsFitted" : self.ROIsFitted,
            "ROIs" : exportROIs,
            "resultsGenerated" : self.resultsGenerated,
            "delayed" : self.delayed
        }

        if self.resultsGenerated:
            outDict["results"] = [fd["results"] for fd in self.fileData]
            outDict["resultHeadings"] = self.fileData[0]["resultHeadings"]
            outDict["evaluatorNames"] = self.fileData[0]["evaluatorNames"]

        if self.delayed and "NAATimes" in self.fileData[0].keys():
            outDict["NAATimes"] = [fd["NAATimes"] for fd in self.fileData]
            
        return outDict

    def add_files(self, files):
        """Parse and add the spectrum files specified in the files list"""
        self.fileList = files
        self.fileData = [SpectrumParser(f).getValues() for f in files]

    def load_known_peaks(self, standardsFilename):
        """Parse and add known peaks from a standards file"""
        self.standardsFilename = standardsFilename
        l = StandardsFileParser(standardsFilename).extract_peaks(self.delayed)
        for p in l:
            c = p.get_ctr() # avoid collisions by changing the center by .01 eV in this dictionary, without affecting the actual object
            while c in self.knownPeaks.keys():
                c += 0.00001
            self.knownPeaks[c] = p

    def update_ROIs(self, addedIsotopes, removedIsotopes = []):
        """Update our ROIs, adding some isotopes and potentially removing some.
        
        This function can be used to create ROIs by calling it with only 1 argument.
        """
        rmvLst = []
        editList = []
        alreadyAdded = []

        if addedIsotopes == [] and removedIsotopes == []:
            #If we are told to do nothing, do nothing
            return None
        else:
            #If we are doing anything, the new ROIs we create might not be fitted.
            self.ROIsFitted = False

        #ensure no duplicate additions
        for iso in addedIsotopes:
            if iso in self.isotopes:
                alreadyAdded.append(iso)
            else:
                self.isotopes.append(iso)

        for iso in alreadyAdded:
            addedIsotopes.remove(iso)

        #remove isotopes that the user wanst to remove
        for iso in removedIsotopes:
            try:
                self.isotopes.remove(iso)
            except:
                pass
        
        #either remove ROIs completely or add to an "edit list" if some, but not all, isotopes in roi have been removed
        for r in self.ROIs:
            isotopes = r.get_isotopes()
            removed = len([x for x in isotopes if x in removedIsotopes])
            if removed == len(isotopes):
                rmvLst.append(r)
            elif removed > 0:
                editList += [kp.get_ctr() for kp in r.get_known_peaks()]
        
        for r in rmvLst:
            self.ROIs.remove(r)

        regions = []
        peaks = []
        otherPeaks = []

        #create new ROIs
        sortedKeys = sorted(self.knownPeaks.keys())

        for k in sortedKeys:
            p = self.knownPeaks[k]
            if p.get_ele() in addedIsotopes or p.get_ctr() in editList:
                if p.get_ele() == "B-11": #special case for boron
                    lowerBound = max(p.get_ctr() - self.userPrefs["B_roi_width"], 0)
                    upperBound = min(p.get_ctr() + self.userPrefs["B_roi_width"], self.fileData[0]["energies"][-1])
                else:
                    lowerBound = max(p.get_ctr() - self.userPrefs["roi_width"], 0)
                    upperBound = min(p.get_ctr() + self.userPrefs["roi_width"], self.fileData[0]["energies"][-1])
                
                regions.append(lowerBound)
                regions.append(upperBound)
                
                peaks.append([p])
                otherPeaks.append([self.knownPeaks[e] for e in sortedKeys[binary_search_find_nearest(sortedKeys, lowerBound):binary_search_find_nearest(sortedKeys, upperBound)]])

        if self.userPrefs["overlap_rois"]:
            i=0
            while i < len(regions) - 1:
                if regions[i] > regions[i+1]: #if there is an overlap, delete both points that overlap, leaving a single, larger region
                    del regions[i]
                    del regions[i]
                    peaks[i//2] += peaks[i//2+1]
                    del peaks[i//2+1]
                    otherPeaks[i//2] += otherPeaks[i//2+1]
                    del otherPeaks[i//2+1]
                else:
                    i += 1

        energies = self.fileData[0]["energies"]
        cps = self.fileData[0]["cps"]
        for i in range(0,len(regions),2):
            lowerIndex = binary_search_find_nearest(energies, regions[i])
            upperIndex = binary_search_find_nearest(energies, regions[i+1])
            r = ROI(energies[lowerIndex:upperIndex],cps[lowerIndex:upperIndex], [lowerIndex, upperIndex], "B-11" in [p.get_ele() for p in peaks[i//2]] and not self.delayed, self.userPrefs)
            r.set_known_peaks(peaks[i//2], otherPeaks[i//2])
            self.ROIs.append(r)
        self.ROIs = sorted(self.ROIs, key=lambda x:x.get_range()[0])

    def get_fitted_ROIs(self):
        """Fits all ROIs that aren't fitted: convenience function that calls several functions on each unfitted ROI."""
        for ROI in self.ROIs:
            if not ROI.fitted:
                ROI.add_peaks()
                ROI.add_bg()
                ROI.fit()
        self.ROIsFitted = True
        return self.ROIs

    def get_entry_repr(self, model, name, ROIIndex, params):
        """Get Entry Represenaation for provided entry, given peak and background type."""
        if model == "peaks":
            testObj = som[model][name]()
            testObj.handle_entry(params, bounds=self.ROIs[ROIIndex].get_range())
            return testObj.to_string(), testObj.get_params()
        elif model == "backgrounds":
            tmpObj = som[model][name].guess_params(self.ROIs[ROIIndex].get_energies(), self.ROIs[ROIIndex].get_cps())
            return tmpObj.to_string(), tmpObj.get_params()

    def set_ROI_range(self, ROIIndex, newRange):
        """Set the range (of energy values) of the ROI at index ROIIndex to the values in values"""
        energies = self.fileData[0]["energies"]
        cps = self.fileData[0]["cps"]
        lowerIndex = binary_search_find_nearest(energies, newRange[0])
        upperIndex = binary_search_find_nearest(energies, newRange[1])
        self.ROIs[ROIIndex].set_data([energies[lowerIndex], energies[upperIndex]], energies[lowerIndex:upperIndex], cps[lowerIndex:upperIndex], [lowerIndex, upperIndex])

    def run_evaluators(self, evaluators, e_args):
        """Run a list of evaluators on our ROIs, with arguments specified in the list e_args"""
        ROIsToEval = [r for r in self.ROIs if r.fitted]
        for i in range(len(self.fileData)):
            if i != 0:
                energies = self.fileData[i]["energies"]
                cps = self.fileData[i]["cps"]
                for r in ROIsToEval:
                    if self.delayed:
                        for kp in r.get_known_peaks():
                            kp.set_delay_times(*self.fileData[i]["NAATimes"], self.fileData[i]["realtime"]/60)
                    bounds = r.get_range()
                    lowerIndex = binary_search_find_nearest(energies, bounds[0])
                    upperIndex = binary_search_find_nearest(energies, bounds[1])
                    r.reanalyze(energies[lowerIndex:upperIndex], cps[lowerIndex:upperIndex])
            self.fileData[i]["results"] = [e(ROIsToEval).get_results(*args) for e, args in zip(evaluators, e_args)]
            self.fileData[i]["resultHeadings"] = [e.get_headings(ROIsToEval[0]) for e in evaluators]
            self.fileData[i]["evaluatorNames"] = [e.get_name() for e in evaluators]
        self.resultsGenerated = True
        return None
    def write_results_file(self, projectID, filename):
        """Writes a results file, format/spec depends on the filename of the request.txt

        Implements ExcelWriter to write results and CSVWriter to write results or spectrum file data.
        """
        if filename.split(".")[-1] == "xlsx":
            headings = [fd["resultHeadings"] for fd in self.fileData]
            data = [fd["results"] for fd in self.fileData]
            ew = ExcelWriter(projectID, self.get_title(), self.fileList, headings, data)
            ew.write()
        elif filename[-21:] == "_Analysis_Results.csv":
            origFilename = filename.replace("_Analysis_Results.csv","")
            for i in range(len(self.fileList)):
                if os.path.split(self.fileList[i])[1].split('.')[0] == origFilename:
                    cw = CSVWriter(projectID, filename, self.fileData[i]["resultHeadings"][0], self.fileData[i]["results"])
                    cw.write()
                    break
        else:
            origFilename = filename.replace("_xy.csv","")
            for i in range(len(self.fileList)):
                if os.path.split(self.fileList[i])[1].split('.')[0] == origFilename:
                    cw = CSVWriter(projectID, filename, ["Energy (keV)", "Counts Per Second"], zip(self.fileData[i]["energies"], self.fileData[i]["cps"]))
                    cw.write()
                    break
    #Getters and Setters
    def set_delayed_times(self, i, irr, wait, count):
        self.fileData[i]["NAATimes"] = [irr, wait, count]
        
    def get_all_isotopes(self):
        return set([self.knownPeaks[key].get_ele() for key in self.knownPeaks.keys()])
    
    def get_known_annots(self):
        return [[[kp.get_ctr(), kp.get_ele()] for kp in r.peaksInRegion] for r in self.ROIs]
    
    def get_naa_times(self):
        return [fd["NAATimes"] for fd in self.fileData]
    
    def get_unfitted_ROIs(self):
        return [i for i in range(len(self.ROIs)) if not self.ROIs[i].fitted]
    
    def set_user_prefs(self, newPrefs):
        for k in newPrefs.keys():
            self.userPrefs[k] = newPrefs[k]

    def get_known_peaks(self):
        return self.knownPeaks
    
    def get_title(self):
        return self.title
    
    def set_title(self, newTitle):
        self.title = newTitle

    def get_isotopes(self):
        return self.isotopes
        
    def get_filename_list(self):
        return [os.path.split(f)[1] for f in self.fileList]

    def get_all_entry_fields(self):
        return {
            "peaks" : {k : som["peaks"][k].get_entry_fields() for k in som["peaks"].keys()},
            "backgrounds" : {k : som["backgrounds"][k].get_entry_fields() for k in som["backgrounds"].keys()}
        }


class ROI:
    def __init__(self, energies, cps, indicies, boronROI = False, userPrefs = default_prefs):
        self.energies = energies
        self.range = (energies[0], energies[-1])
        self.cps = cps
        self.knownPeaks = []
        self.peaksInRegion = []
        self.userPrefs = userPrefs
        self.indicies = indicies
        self.peaks = []
        self.peakPairs = None
        self.fitted = False
        self.boronROI = boronROI

    def load_from_dict(self, stored_data):
        """Sets variables for an ROI object based on a dictionary exported by the export_to_dict() function."""
        if "peaks" in stored_data.keys():
            self.peaks = [som["peaks"][p["type"]](*p["params"], variances=p["variances"]) for p in stored_data["peaks"]]
            self.bg = som["backgrounds"][stored_data["background"]["type"]](*stored_data["background"]["params"],variances=stored_data["background"]["variances"])
            self.fitted = (self.peaks[0].get_variances()[0] != None)
        else:
            self.fitted = False
        for kp in stored_data["knownPeaks"]:
            knownPeakObj = KnownPeak()
            knownPeakObj.load_from_dict(kp)
            self.knownPeaks.append(knownPeakObj)
    def export_to_dict(self):
        """Exports the current ROI state to a dictionary"""
        PIR = [p.export_to_dict() for p in self.peaksInRegion]
        try:
            exportPeaks = [{"type" : p.get_type(), "params" : p.get_original_params(), "variances": p.get_original_variances()} for p in self.peaks]
            exportBackground = {"type" : self.bg.get_type(), "params" : self.bg.get_original_params(), "variances": self.bg.get_original_variances()}
            exportKnownPeaks = [kp.export_to_dict() for kp in self.knownPeaks]
            return {
                "indicies" : self.indicies,
                "peaks" : exportPeaks,
                "background" : exportBackground,
                "knownPeaks" : exportKnownPeaks,
                "peaksInRegion" : PIR
            }
        except:
            exportKnownPeaks = [kp.export_to_dict() for kp in self.knownPeaks]
            return {
                "indicies" : self.indicies,
                "knownPeaks" : exportKnownPeaks,
                "peaksInRegion" : PIR
            }
    def add_peaks(self):
        """Find and add peaks to own model (guesss params)"""
        if self.boronROI:
            BPeak = som["peaks"][self.userPrefs["boron_peak_type"]]
            self.peaks = [BPeak.guess_params(self.energies, self.cps)]
            self.peaks += som["peaks"][self.userPrefs["peak_type"]].guess_params(self.energies, BPeak.remove_from_data(self.energies, self.cps))
        else:
            self.peaks = som["peaks"][self.userPrefs["peak_type"]].guess_params(self.energies, self.cps)
   
    def add_bg(self):
        """Find and add background to own model (guesss params)"""
        self.bg = som["backgrounds"][self.userPrefs["background_type"]].guess_params(self.energies, self.cps)

    def fit(self, reanalyze = False):
        """Fit our model to the data within the ROI, using the guessed params as initial ones"""
        f = lambda x,*params: multiple_peak_and_background(self.peaks, self.bg, x, params)
        p0 = np.array(self.bg.get_params() + list(itertools.chain.from_iterable([p.get_params() for p in self.peaks])))
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            try:
                params, cov = curve_fit(f, self.energies, self.cps, p0=p0)
                variances = np.diag(cov)
                set_all_params(self.peaks, self.bg, params, variances, reanalyze)
                self.fitted = True
            except Exception as e:
                print(e)
                self.fitted = False
                pass
    
    def get_fitted_curve(self, xdata = None):
        """Get the output of our fit (ydata) given x values"""
        if xdata == None:
            xdata = np.arange(self.range[0], self.range[-1], .01)
        return [list(xdata), get_curve(self.peaks, self.bg, xdata), list(self.bg.get_ydata(xdata))]

    def get_closest_peak(self, kp):
        """Find the closest peak in our data to a given known peak. Works well for auto matching."""
        c = kp.get_ctr()
        minSep = 99999
        minPeak = None
        for p in self.peaks:
            if(abs(p.get_ctr() - c) < minSep):
                minSep = abs(p.get_ctr() - c)
                minPeak = p
        return minPeak
    
    def set_original_peak_pairs(self, energyPairs):
        """Original peak pairs are set so that they don't change when the ROI is reanalyzed on new data and can be exported/imported easily."""
        pairs = []
        for ep in energyPairs:
            pair = []
            for p in self.peaks:
                if p.get_ctr() == ep[0]:
                    pair.append(p)
                    break
            for kp in self.knownPeaks:
                if kp.get_ctr() == ep[1]:
                    pair.append(kp)
                    break
            pairs.append(pair)
        self.peakPairs = pairs
        self.originalPeakPairs = pairs

    def reanalyze(self, newEnergies, newCPS):
        """Re-runs the fit on a new set of energies and cps from another spectrum file, and re-match peaks"""
        if self.peakPairs == None:
            raise RuntimeError("Reanalyze called before peak pairs created!")
        self.energies = newEnergies
        self.cps = newCPS
        self.fit(True)
        outputs = []
        for p in self.originalPeakPairs:
            closestMatch = self.peaks[min(range(len(self.peaks)), key = lambda i: abs(self.peaks[i].get_ctr()-p[0].get_ctr()))]
            outputs.append([closestMatch, p[1]])  
        self.peakPairs = outputs      
        return outputs
    
    #Getters and Setters
    def set_peaks(self, peaks):
        self.peaks = peaks
    
    def get_peaks(self):
        return self.peaks

    def get_isotopes(self):
        return [kp.get_ele() for kp in self.knownPeaks]

    def get_peak_ctrs(self):
        return [p.get_ctr() for p in self.peaks]

    def get_known_peaks(self):
        return self.knownPeaks
    
    def get_range(self):
        return list(self.range)

    def get_formatted_range(self):
        return [str(round(float(self.range[0]), decimals=1)), str(round(float(self.range[1]), decimals=1))]

    def set_range(self, newRange):
        self.range = newRange
        self.energies = np.arange(newRange[0], newRange[1], .01)

    def get_energies(self):
        return list(self.energies)

    def get_cps(self):
        return list(self.cps)

    def set_data(self, newRange, energies, cps, indicies):
        self.range = newRange
        self.energies = energies
        self.cps = cps
        self.indicies = indicies

    def set_known_peaks(self, peaks, otherPeaks):
        self.knownPeaks = peaks
        self.peaksInRegion = otherPeaks

    def set_background(self, bg):
        self.bg = bg

    def get_background(self):
        return self.bg

    def get_peak_pairs(self):
        return self.peakPairs

    def get_indicies(self):
        return self.indicies





