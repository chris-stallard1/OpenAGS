import os
import warnings
import itertools

import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
from sigfig import round

from openags.util import multiple_peak_and_background, get_curve, binary_search_find_nearest, set_all_params, KnownPeak
from openags.constants import default_prefs, som
from openags.parsers import SpectrumParser, StandardsFileParser, CSVWriter, ExcelWriter
from openags.models import GaussianPeak, LinearBackground



class ActivationAnalysis:
    """An application class representing the whole backend. 
    
    This is the ONLY backend class that the frontend should interact with.

    ----Variables----
    #User Preferences
    userPrefs: Dict[str, any] = None
    #Extracted contents and metadata from project files
    fileData: List[Dict[str, any]] = None
    #List of spectrum filenames, data from each correlates to fileData
    fileList: List[str] = None
    #dictionary mapping peak location to KnownPeak onjects
    knownPeaks: Dict[float, KnownPeak] = None
    #Regions of Interest, see ROI class below
    ROIs: List[ROI] = None
    #Isotopes being analyzed
    isotopes: List[str] = None
    #Project Title
    title: str
    #Whether all regions of interest for this project have been fitted
    ROIsFitted: bool
    #Whether evaluators have been run and results generated for this project
    resultsGenerated: bool
    #Whether or not this is a Delayed Gamma Analysis
    delayed: bool"""
    
    def __init__(self, userPrefs = default_prefs, title = ""):
        self.userPrefs = userPrefs
        self.fileData = None
        self.fileList = None
        self.knownPeaks = None
        self.ROIs = None
        self.isotopes = None
        self.title = title
        self.ROIsFitted = False
        self.resultsGenerated = False
        self.delayed = False

    def load_from_dict(self, stored_data):
        """Sets variables for an analysis object based on a dictionary exported by the export_to_dict() function."""

        if "userPrefs" in stored_data: #otherwise keep default prefs
            self.userPrefs = stored_data["userPrefs"]
        
        self.title = stored_data["title"]
        self._load_spectra(stored_data["files"])
        self._load_known_peaks(stored_data["standardsFilename"])
        energies = self.fileData[0]["energies"]
        cps = self.fileData[0]["cps"]

        self.ROIsFitted = stored_data["ROIsFitted"]
        self.ROIs = []
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
        if self.delayed and "NAATimes" in stored_data:
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

        if self.delayed and "NAATimes" in self.fileData[0]:
            outDict["NAATimes"] = [fd["NAATimes"] for fd in self.fileData]
            
        return outDict

    def _load_spectra(self, files):
        """Parse and add the spectrum files specified in the files list"""
        self.fileList = files
        self.fileData = [SpectrumParser(f).getValues() for f in files]

    def _load_known_peaks(self, standardsFilename):
        """Parse and add known peaks from a standards file"""
        self.standardsFilename = standardsFilename
        self.knownPeaks = {}
        peakList = StandardsFileParser(standardsFilename).extract_peaks(self.delayed)
        for p in peakList:
            c = p.get_ctr() # avoid collisions by changing the center by .01 eV in this dictionary, without affecting the actual object
            while c in self.knownPeaks:
                c += 0.00001
                
            self.knownPeaks[c] = p

    def update_ROIs(self, addedIsotopes, removedIsotopes = []):
        """Update our ROIs, adding some isotopes and potentially removing some.
        
        This function can be used to create ROIs by calling it with only 1 argument.
        """

        if not addedIsotopes and not removedIsotopes:
            #If we are told to do nothing, do nothing
            return

        #If we are doing anything, the new ROIs we create might not be fitted.
        self.ROIsFitted = False

        #ensure no duplicate additions
        addedIsotopes = [iso for iso in addedIsotopes if iso not in self.isotopes and iso in self.get_all_isotopes()]
        self.isotopes.extend(addedIsotopes)

        #remove isotopes that the user wants to remove
        for iso in removedIsotopes:
            try:
                self.isotopes.remove(iso)
            except ValueError: #if there is a duplicate, move on
                pass
        
        #either remove ROIs completely or add to an "edit list" if some, but not all, isotopes in roi have been removed
        editList, ROIs = [], []
        for r in self.ROIs:
            isotopes = r.get_isotopes()
            filtered = [iso for iso in isotopes if iso not in removedIsotopes]
            if filtered:
                # some isotopes remaining after remove, so keep the ROI.
                ROIs.append(r)
                if len(filtered) < len(isotopes):
                    editList.extend(kp.get_ctr() for kp in r.get_known_peaks)
        self.ROIs = ROIs

        regions = []
        peaks = []
        otherPeaks = []

        #create new ROIs
        sortedKeys = sorted(self.knownPeaks.keys())

        for k in sortedKeys:
            p = self.knownPeaks[k]
            if p.get_ele() in addedIsotopes or p.get_ctr() in editList:
                if p.get_ele() == "B-11" and p.get_ctr() < 480 and p.get_ctr() > 470: #special case for boron
                    lowerBound = max(p.get_ctr() - self.userPrefs["B_roi_width"], 0)
                    upperBound = min(p.get_ctr() + self.userPrefs["B_roi_width"], self.fileData[0]["energies"][-1])
                else:
                    lowerBound = max(p.get_ctr() - self.userPrefs["roi_width"], 0)
                    upperBound = min(p.get_ctr() + self.userPrefs["roi_width"], self.fileData[0]["energies"][-1])
                
                regions.append(lowerBound)
                regions.append(upperBound)
                
                peaks.append([p])
                lowerIndex = binary_search_find_nearest(sortedKeys, lowerBound)
                upperIndex = binary_search_find_nearest(sortedKeys, upperBound)
                otherPeaks.append([self.knownPeaks[e] for e in sortedKeys[lowerIndex:upperIndex]])

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
            boronRegion = "B-11" in [p.get_ele() for p in peaks[i//2]] and not self.delayed and regions[i] < 477.6 and regions[i+1] > 477.6
            r = ROI(energies[lowerIndex:upperIndex],cps[lowerIndex:upperIndex], [lowerIndex, upperIndex], boronRegion, self.userPrefs)
            r.set_known_peaks(peaks[i//2], otherPeaks[i//2])
            self.ROIs.append(r)
        self.ROIs = sorted(self.ROIs, key=lambda x:x.get_range()[0])

    def fit_ROIs(self):
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
            return str(testObj), testObj.get_params()
        elif model == "backgrounds":
            tmpObj = som[model][name].guess_params(self.ROIs[ROIIndex].get_energies(), self.ROIs[ROIIndex].get_cps())
            return str(tmpObj), tmpObj.get_params()

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

    def write_results_file(self, projectID, filename):
        """Writes a results file, format/spec depends on the filename of the request.

        Implements ExcelWriter to write results and CSVWriter to write results or spectrum file data.
        """
        if filename.split(".")[-1] == "xlsx":
            headings = [fd["resultHeadings"] for fd in self.fileData]
            data = [fd["results"] for fd in self.fileData]
            ew = ExcelWriter(projectID, self.get_title(), self.fileList, headings, data)
            ew.write()
        elif filename[-21:] == "_Analysis_Results.csv":
            origFilename = filename.replace("_Analysis_Results.csv","")
            for filename, fileData in zip(self.fileList, self.fileData):
                if os.path.split(filename)[1].split('.')[0] == origFilename:
                    cw = CSVWriter(projectID, filename, fileData["resultHeadings"][0], fileData["results"])
                    cw.write()
                    break
        elif filename[-7:] == "_xy.csv":
            origFilename = filename.replace("_xy.csv","")
            for filename, fileData in zip(self.fileList, self.fileData):
                if os.path.split(filename)[1].split('.')[0] == origFilename:
                    cw = CSVWriter(projectID, filename, ["Energy (keV)", "Counts Per Second"], zip(fileData["energies"], fileData["cps"]))
                    cw.write()
                    break
    #Getters and Setters
    def set_delayed_times(self, i, irr, wait, count):
        self.fileData[i]["NAATimes"] = [irr, wait, count]
        
    def get_all_isotopes(self):
        return set(v.get_ele() for k, v in self.knownPeaks.items())
    
    def get_known_annots(self):
        return [[[kp.get_ctr(), kp.get_ele()] for kp in r.peaksInRegion] for r in self.ROIs]
    
    def get_naa_times(self):
        return [fd["NAATimes"] for fd in self.fileData]
    
    def get_unfitted_ROIs(self):
        return [i for i, ROI in enumerate(self.ROIs) if not ROI.fitted]
    
    def set_user_prefs(self, newPrefs):
        self.userPrefs.update(newPrefs)

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
            "peaks" : {k : v.get_entry_fields() for k, v in som["peaks"].items()},
            "backgrounds" : {k : v.get_entry_fields() for k, v in som["backgrounds"].items()}
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
        self.bg = None
        self.peakPairs = None
        self.fitted = False
        self.boronROI = boronROI

    def load_from_dict(self, stored_data):
        """Sets variables for an ROI object based on a dictionary exported by the export_to_dict() function."""
        if "peaks" in stored_data:
            self.peaks = [som["peaks"][p["type"]](*p["params"], variances=p["variances"]) for p in stored_data["peaks"]]
            self.bg = som["backgrounds"][stored_data["background"]["type"]](*stored_data["background"]["params"],variances=stored_data["background"]["variances"])
            self.fitted = (self.peaks[0].get_variances()[0] != None)
        else:
            self.fitted = False
        
        self.knownPeaks = [KnownPeak().load_from_dict(kp) for kp in stored_data["knownPeaks"]] 
        if "peakPairs" in stored_data:
            self.peakPairs = self.originalPeakPairs = [(self.peaks[i], self.knownPeaks[j]) for i, j in stored_data["peakPairs"]]
    def export_to_dict(self):
        """Exports the current ROI state to a dictionary"""
        PIR = [p.export_to_dict() for p in self.peaksInRegion]
        exportKnownPeaks = [kp.export_to_dict() for kp in self.knownPeaks]
        outDict = {
                "indicies" : self.indicies,
                "knownPeaks" : exportKnownPeaks,
                "peaksInRegion" : PIR
        }

        if self.bg != None:
            outDict["peaks"] = [{"type" : p.get_type(), "params" : p.get_original_params(), "variances": p.get_original_variances()} for p in self.peaks]
            outDict["background"] = {"type" : self.bg.get_type(), "params" : self.bg.get_original_params(), "variances": self.bg.get_original_variances()}
        if self.peakPairs != None:
            outDict["peakPairs"] = [(self.peaks.index(p), self.knownPeaks.index(kp)) for p, kp in self.originalPeakPairs]
        return outDict

    def add_peaks(self):
        """Find and add peaks to own model (guesss params)"""
        if self.boronROI:
            BPeak = som["peaks"][self.userPrefs["boron_peak_type"]]
            self.peaks = [BPeak.guess_params(self.energies, self.cps)]
            scrubbedCPS = BPeak.remove_from_data(self.energies, self.cps)
            self.peaks += som["peaks"][self.userPrefs["peak_type"]].guess_params(self.energies, scrubbedCPS)
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

    def get_closest_peak(self, peak):
        if len(self.peaks) == 0:
            return None
        target = peak.get_ctr()
        distance = [abs(p.get_ctr() - target) for p in self.peaks]
        index = np.argmin(distance)
        return self.peaks[index]
    
    def set_original_peak_pairs(self, energyPairs):
        """Original peak pairs are set so that they don't change when the ROI is reanalyzed on new data and can be exported/imported easily."""
        pairs = []
        peakCtr = np.array([p.get_ctr() for p in self.peaks])
        knownCtr = np.array([p.get_ctr() for p in self.knownPeaks])
        nearest_peak = lambda target: self.peaks[np.argmin(abs(peakCtr-target))]
        nearest_known = lambda target: self.knownPeaks[np.argmin(abs(knownCtr-target))]
        pairs = [(nearest_peak(p), nearest_known(kp)) for p, kp in energyPairs]
        self.peakPairs = self.originalPeakPairs = pairs

    def reanalyze(self, newEnergies, newCPS):
        """Re-runs the fit on a new set of energies and cps from another spectrum file, and re-match peaks"""
        if self.peakPairs == None:
            raise RuntimeError("Reanalyze called before peak pairs created!")
        self.energies = newEnergies
        self.cps = newCPS
        self.fit(True)
        outputs = []
        peakCtrs = np.array([p.get_ctr() for p in self.peaks])
        for peak, knownPeak in self.originalPeakPairs:
            target = peak.get_ctr()
            closestMatch = self.peaks[np.argmin(abs(peakCtrs - target))]
            outputs.append([closestMatch, knownPeak])  
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





