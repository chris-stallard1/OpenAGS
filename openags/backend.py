import os
import warnings
import itertools
from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
from sigfig import round
from collections.abc import Iterable

from openags.util import (
    multiple_peak_and_background,
    get_curve,
    binary_search_find_nearest,
    set_all_params,
    KnownPeak,
)
from openags.constants import default_prefs, som
from openags.parsers import SpectrumParser, StandardsFileParser, CSVWriter, ExcelWriter

from openags.baseClasses import Background


@dataclass
class ActivationAnalysis:
    """An application class representing the whole backend.

    This is the ONLY backend class that the frontend should interact with.

    ----Variables----
    #User Preferences
    user_prefs: Dict[str, any] = None
    #Extracted contents and metadata from project files
    file_data: List[Dict[str, any]] = None
    #List of spectrum filenames, data from each correlates to file_data
    file_list: List[str] = None
    #dictionary mapping peak location to KnownPeak onjects
    known_peaks: Dict[float, KnownPeak] = None
    #Regions of Interest, see ROI class below
    ROIs: List[ROI] = None
    #Isotopes being analyzed
    isotopes: List[str] = None
    #Project Title
    title: str
    #Whether all regions of interest for this project have been fitted
    ROIs_fitted: bool
    #Whether evaluators have been run and results generated for this project
    results_generated: bool
    #Whether or not this is a Delayed Gamma Analysis
    delayed: bool"""

    file_data: list
    file_list: list
    known_peaks: dict
    batches: list
    primary_files: list
    region_data: list
    isotopes: list
    title: str
    k0_fname: str
    eff_fname: str
    flux: float
    user_prefs: dict
    delayed: bool = False

    @staticmethod
    def load_from_dict(stored_data):
        """Sets variables for an analysis object based on a dictionary exported by the export_to_dict() function."""

        if "user_prefs" in stored_data:  # otherwise keep default prefs
            user_prefs = stored_data["user_prefs"]
        else:
            user_prefs = default_prefs

        title = stored_data["title"]
        file_list = stored_data["files"]
        file_data = [SpectrumParser(f).getValues() for f in file_list]
        delayed = stored_data["delayed"]
        known_peaks = ActivationAnalysis.load_known_peaks(stored_data["k0_fname"], stored_data["eff_fname"], delayed, stored_data["flux"])
        region_data = [ROIData.load_from_dict(d) for d in stored_data["region_data"]]
        isotopes = list(
            set(itertools.chain(*[r.get_isotopes() for r in region_data]))
        )

        if delayed and "NAATimes" in stored_data:
            for i in range(len(stored_data["NAATimes"])):
                file_data[i]["NAATimes"] = stored_data["NAATimes"][i]

        a = ActivationAnalysis(file_data, file_list, known_peaks, [], [], region_data, isotopes, title,
                               stored_data["k0_fname"], stored_data["eff_fname"], stored_data["flux"], user_prefs, delayed)
        a.batches = [BatchAnalysis.load_from_dict(a, d) for d in stored_data["batches"]]
        a.primary_files = [b.pf_number for b in a.batches]
        return a

    def export_to_dict(self):
        """Exports the current project state to a dictionary"""
        export_batches = [b.export_to_dict() for b in self.batches]
        export_regions = [r.export_to_dict() for r in self.region_data]
        outDict = {
            "user_prefs": self.user_prefs,
            "title": self.title,
            "files": self.file_list,
            "k0_fname": self.k0_fname,
            "eff_fname": self.eff_fname,
            "flux": self.flux,
            "batches": export_batches,
            "region_data": export_regions,
            "delayed": self.delayed,
        }

        if self.delayed and "NAATimes" in self.file_data[0]:
            outDict["NAATimes"] = [fd["NAATimes"] for fd in self.file_data]

        return outDict

    @staticmethod
    def load_known_peaks(k0_fname, eff_fname, delayed, flux=1):
        """Parse and add known peaks from a standards file"""
        known_peaks = {}
        peakList = StandardsFileParser(k0_fname, eff_fname, flux).extract_peaks(delayed)
        for p in peakList:
            c = (
                p.get_ctr()
            )  # avoid collisions by changing the center by .01 eV in this dictionary, without affecting the actual object
            while c in known_peaks:
                c += 0.00001

            known_peaks[c] = p
        return known_peaks

    def update_ROIs(self, added_isotopes, removed_isotopes=()):
        """Update our ROIs, adding some isotopes and potentially removing some.

        This function can be used to create ROIs by calling it with only 1 argument.
        """

        if not added_isotopes and not removed_isotopes:
            # If we are told to do nothing, do nothing
            return

        # If we are doing anything, the new ROIs we create might not be fitted.
        for b in self.batches:
            b.ROIs_fitted = False

        # ensure no duplicate additions
        added_isotopes = [
            iso for iso in added_isotopes
            if iso not in self.isotopes and iso in self.get_all_isotopes()
        ]
        self.isotopes.extend(added_isotopes)

        # remove isotopes that the user wants to remove
        for iso in removed_isotopes:
            try:
                self.isotopes.remove(iso)
            except ValueError:  # if there is a duplicate, move on
                pass

        # either remove ROIs completely or add to an "edit list" if some, but not all, isotopes in roi have been removed
        editList, new_region_data = [], []
        for r in self.region_data:
            isotopes = r.get_isotopes()
            filtered = [iso for iso in isotopes if iso not in removed_isotopes]
            if filtered:
                # some isotopes remaining after remove, so keep the ROI.
                new_region_data.append(r)
                if len(filtered) < len(isotopes):
                    editList.extend(kp.get_ctr() for kp in r.get_known_peaks)
            else:
                for b in self.batches:
                    b.remove_ROI(r)

        self.region_data = new_region_data

        regions = []
        peaks = []
        other_peaks = []

        # create new ROIs
        sortedKeys = sorted(self.known_peaks.keys())

        for k in sortedKeys:
            p = self.known_peaks[k]
            if p.get_ele() in added_isotopes or p.get_ctr() in editList:
                if (
                        p.get_ele() == "B-11" and p.get_ctr() < 480 and p.get_ctr() > 470
                ):  # special case for boron
                    lowerBound = max(
                        p.get_ctr() - self.user_prefs["B-11 ROI Width (keV)"], 0
                    )
                    upperBound = min(
                        p.get_ctr() + self.user_prefs["B-11 ROI Width (keV)"],
                        self.file_data[0]["energies"][-1],
                    )
                else:
                    lowerBound = max(p.get_ctr() - self.user_prefs["ROI Width (keV)"], 0)
                    upperBound = min(
                        p.get_ctr() + self.user_prefs["ROI Width (keV)"],
                        self.file_data[0]["energies"][-1],
                    )

                regions.append(lowerBound)
                regions.append(upperBound)

                peaks.append([p])
                lowerIndex = binary_search_find_nearest(sortedKeys, lowerBound)
                upperIndex = binary_search_find_nearest(sortedKeys, upperBound)
                other_peaks.append(
                    [self.known_peaks[e] for e in sortedKeys[lowerIndex:upperIndex]]
                )

        if self.user_prefs["Overlap ROIs"]:
            i = 0
            while i < len(regions) - 1:
                # if there is an overlap, delete both points that overlap, leaving a single, larger region
                if regions[i] > regions[i + 1]:
                    del regions[i]
                    del regions[i]
                    peaks[i // 2] += peaks[i // 2 + 1]
                    del peaks[i // 2 + 1]
                    other_peaks[i // 2] += other_peaks[i // 2 + 1]
                    del other_peaks[i // 2 + 1]
                else:
                    i += 1

        for i in range(0, len(regions), 2):
            boron_region = (
                    "B-11" in [p.get_ele() for p in peaks[i // 2]]
                    and not self.delayed
                    and regions[i] < 477.6
                    and regions[i + 1] > 477.6
            )
            r = ROIData(
                min_e=regions[i],
                max_e=regions[i+1],
                known_peaks=peaks[i//2],
                region_peaks=other_peaks[i//2],
                boron_region=boron_region
            )
            for b in self.batches:
                b.add_ROI(r)
            self.region_data.append(r)

    def get_entry_repr(self, model, name, batch_num, region_index, params):
        """Get Entry Representation for provided entry, given peak and background type."""
        batch = self.get_batch(batch_num)
        if model == "peaks":
            testObj = som[model][name]()
            testObj.handle_entry(params, bounds=batch.ROIs[region_index].get_range())
            return str(testObj), testObj.get_params()
        elif model == "backgrounds":
            tmpObj = som[model][name].guess_params(
                batch.ROIs[region_index].get_energies(), batch.ROIs[region_index].get_cps()
            )
            return str(tmpObj), tmpObj.get_params()

    def write_results_file(self, projectID, filename):
        """Writes a results file, format/spec depends on the filename of the request.

        Implements ExcelWriter to write results and CSVWriter to write results or spectrum file data.
        """
        if filename.split(".")[-1] == "xlsx":
            headings = [fd["resultHeadings"] for fd in self.file_data]
            data = [fd["results"] for fd in self.file_data]
            ew = ExcelWriter(projectID, self.get_title(), self.file_list, headings, data)
            ew.write()
        elif filename[-21:] == "_Analysis_Results.csv":
            origFilename = filename.replace("_Analysis_Results.csv", "")
            for filename, fileData in zip(self.file_list, self.file_data):
                if os.path.split(filename)[1].split(".")[0] == origFilename:
                    cw = CSVWriter(
                        projectID,
                        filename,
                        fileData["resultHeadings"][0],
                        fileData["results"],
                    )
                    cw.write()
                    break
        elif filename[-7:] == "_xy.csv":
            origFilename = filename.replace("_xy.csv", "")
            for filename, fileData in zip(self.file_list, self.file_data):
                if os.path.split(filename)[1].split(".")[0] == origFilename:
                    cw = CSVWriter(
                        projectID,
                        filename,
                        ["Energy (keV)", "Counts Per Second"],
                        zip(fileData["energies"], fileData["cps"]),
                    )
                    cw.write()
                    break

    # Getters and Setters
    def get_batch(self, file_num):
        for b in self.batches:
            if b.pf_number == file_num:
                return b
        return None

    def get_unanalyzed_batch(self):
        for b in self.batches:
            if not b.results_generated:
                return b.pf_number
        return self.batches[0].pf_number

    def all_analyzed(self):
        return all([b.results_generated for b in self.batches])

    def get_analyzed_batches(self):
        return [b for b in self.batches if b.results_generated]

    def set_delayed_times(self, i, irr, wait, count):
        self.file_data[i]["NAATimes"] = [irr, wait, count]

    def get_all_isotopes(self):
        return set(v.get_ele() for k, v in self.known_peaks.items())

    def get_known_annots(self):
        return [
            [[kp.get_ctr(), kp.get_ele()] for kp in r.region_peaks] for r in self.region_data
        ]

    def get_naa_times(self):
        return [fd["NAATimes"] for fd in self.file_data]

    def set_user_prefs(self, new_prefs):
        self.user_prefs.update(new_prefs)

    def get_known_peaks(self):
        return self.known_peaks

    def get_title(self):
        return self.title

    def set_title(self, new_title):
        self.title = new_title

    def get_isotopes(self):
        return self.isotopes

    def get_filename_list(self):
        return [os.path.split(f)[1] for f in self.file_list]

    def get_all_entry_fields(self):
        return {
            "peaks": {k: v.get_entry_fields() for k, v in som["peaks"].items()},
            "backgrounds": {
                k: v.get_entry_fields() for k, v in som["backgrounds"].items()
            },
        }

@dataclass
class BatchAnalysis:
    parent_analysis: ActivationAnalysis
    pf_number: int
    primary_file: dict
    reanalysis_files: list
    results_generated: bool
    ROIs: list
    ROIs_fitted: bool

    @staticmethod
    def load_from_dict(parent_analysis, stored_data):
        pf_number = stored_data["pf_number"]
        primary_file = parent_analysis.file_data[pf_number]
        reanalysis_files = stored_data["reanalysis_files"]
        ROIs = []
        for i, ROI_data in enumerate(stored_data["ROIs"]):
            lowerIndex = ROI_data["indicies"][0]
            upperIndex = ROI_data["indicies"][1]
            r = ROI.load_from_dict(
                primary_file["energies"][lowerIndex:upperIndex],
                primary_file["cps"][lowerIndex:upperIndex],
                [lowerIndex, upperIndex],
                parent_analysis.region_data[i],
                ROI_data,
                parent_analysis.user_prefs
            )
            ROIs.append(r)
        ROIs_fitted = all([r.fitted for r in ROIs])

        results_generated = stored_data["results_generated"]
        if results_generated:
            fd = parent_analysis.file_data
            rf = [pf_number] + reanalysis_files
            for i in range(len(stored_data["results"])):
                fd[rf[i]]["results"] = stored_data["results"][i]
                fd[rf[i]]["resultHeadings"] = stored_data["resultHeadings"]
                fd[rf[i]]["evaluatorNames"] = stored_data["evaluatorNames"]

        return BatchAnalysis(parent_analysis, pf_number, primary_file, reanalysis_files, results_generated, ROIs,
                             ROIs_fitted)

    def export_to_dict(self):
        out_dict = {"ROIs": [r.export_to_dict() for r in self.ROIs], "pf_number": self.pf_number,
                    "reanalysis_files": self.reanalysis_files, "results_generated": self.results_generated}

        if self.results_generated:
            fd = self.parent_analysis.file_data
            out_dict["results"] = [self.primary_file["results"]] + [fd[i]["results"] for i in self.reanalysis_files]
            out_dict["resultHeadings"] = self.primary_file["resultHeadings"]
            out_dict["evaluatorNames"] = self.primary_file["evaluatorNames"]

        return out_dict

    def run_evaluators(self, evaluators, e_args):
        """Run a list of evaluators on our ROIs, with arguments specified in the list e_args"""
        fitted_ROIs = [r for r in self.ROIs if r.fitted]
        self.primary_file["results"] = [e(fitted_ROIs).get_results(*args) for e, args in zip(evaluators, e_args)]
        self.primary_file["resultHeadings"] = [e.get_headings(fitted_ROIs[0]) for e in evaluators]
        self.primary_file["evaluatorNames"] = [e.get_name() for e in evaluators]

        for i in self.reanalysis_files:
            successfulROIs = []
            energies = self.parent_analysis.file_data[i]["energies"]
            cps = self.parent_analysis.file_data[i]["cps"]
            for r in fitted_ROIs:
                if self.parent_analysis.delayed:
                    for kp in r.get_known_peaks():
                        kp.set_delay_times(
                            *self.parent_analysis.fileData[i]["NAATimes"],
                            self.parent_analysis.fileData[i]["realtime"] / 60
                        )
                bounds = r.get_range()
                lowerIndex = binary_search_find_nearest(energies, bounds[0])
                upperIndex = binary_search_find_nearest(energies, bounds[1])
                try:
                    r.reanalyze(
                        energies[lowerIndex:upperIndex], cps[lowerIndex:upperIndex]
                    )
                    successfulROIs.append(r)
                except Exception:
                    continue
            self.parent_analysis.file_data[i]["results"] = [
                e(successfulROIs).get_results(*args) for e, args in zip(evaluators, e_args)
            ]
            self.parent_analysis.file_data[i]["resultHeadings"] = [
                e.get_headings(fitted_ROIs[0]) for e in evaluators
            ]
            self.parent_analysis.file_data[i]["evaluatorNames"] = [
                e.get_name() for e in evaluators
            ]
        self.results_generated = True

    def add_ROI(self, ROI_data):
        low_index = binary_search_find_nearest(self.primary_file["energies"], ROI_data.min_e)
        high_index = binary_search_find_nearest(self.primary_file["energies"], ROI_data.max_e)
        r = ROI(
            self.primary_file["energies"][low_index:high_index],
            (self.primary_file["energies"][low_index], self.primary_file["energies"][high_index]),
            self.primary_file["cps"][low_index:high_index],
            ROI_data,
            ROI_data.known_peaks,
            ROI_data.region_peaks,
            self.parent_analysis.user_prefs,
            (low_index, high_index),
            (),
            None,
            (),
            (),
            False,
            ROI_data.boron_region
        )
        self.ROIs.append(r)

    def add_all_ROIs(self):
        for ROI in self.parent_analysis.region_data:
            self.add_ROI(ROI)

    def remove_ROI(self, ROI_data):
        to_rmv = [r for r in self.ROIs if r.shared_data == ROI_data]
        self.ROIs.remove(to_rmv)

    def set_ROI_range(self, ROI_index, new_range):
        """Set the range (of energy values) of the ROI at index ROI_index to the values in values"""
        energies = self.primary_file["energies"]
        cps = self.primary_file["cps"]
        lowerIndex = binary_search_find_nearest(energies, new_range[0])
        upperIndex = binary_search_find_nearest(energies, new_range[1])
        self.ROIs[ROI_index].set_data(
            [energies[lowerIndex], energies[upperIndex]],
            energies[lowerIndex:upperIndex],
            cps[lowerIndex:upperIndex],
            [lowerIndex, upperIndex],
        )

    def fit_ROIs(self):
        """Fits all ROIs that aren't fitted: convenience function that calls several functions on each unfitted ROI."""
        for ROI in self.ROIs:
            if not ROI.fitted:
                ROI.add_peaks()
                ROI.add_bg()
                ROI.fit()
        self.ROIs_fitted = True
        return self.ROIs

    def get_unfitted_ROIs(self):
        return [i for i, ROI in enumerate(self.ROIs) if not ROI.fitted]

    # Basic Python Overrides
    def __contains__(self, item):
        return item == self.pf_number or item in self.reanalysis_files

    def __getitem__(self, index):
        if index == 0:
            return self.pf_number
        return self.reanalysis_files[index - 1]

@dataclass
class ROIData:
    min_e: float
    max_e: float
    known_peaks: list
    region_peaks: list
    boron_region: bool

    def get_formatted_range(self):
        return [
            str(round(float(self.min_e), decimals=1)),
            str(round(float(self.max_e), decimals=1)),
        ]

    def get_isotopes(self):
        return [kp.get_ele() for kp in self.known_peaks]

    def export_to_dict(self):
        return {
            "min_e": self.min_e,
            "max_e": self.max_e,
            "known_peaks": [kp.export_to_dict() for kp in self.known_peaks],
            "region_peaks": [rp.export_to_dict() for rp in self.region_peaks],
            "boron_region": self.boron_region
        }

    @staticmethod
    def load_from_dict(data_dict):
        kp = [KnownPeak.load_from_dict(p) for p in data_dict["known_peaks"]]
        rp = [KnownPeak.load_from_dict(p) for p in data_dict["region_peaks"]]
        return ROIData(data_dict["min_e"], data_dict["max_e"], kp, rp, data_dict["boron_region"])


@dataclass
class ROI:
    energies: np.array
    range: tuple
    cps: np.array
    shared_data: ROIData
    known_peaks: list
    region_peaks: list
    user_prefs: dict
    indicies: tuple
    peaks: Iterable
    bg: Background
    peakPairs: Iterable
    originalPeakPairs: Iterable
    fitted: bool
    boron_ROI: bool

    @staticmethod
    def load_from_dict(energies, cps, indicies, shared_data, stored_data, user_prefs):
        """Sets variables for an ROI object based on a dictionary exported by the export_to_dict() function."""
        if "peaks" in stored_data:
            peaks = [
                som["peaks"][p["type"]](*p["params"], variances=p["variances"])
                for p in stored_data["peaks"]
            ]
            bg = som["backgrounds"][stored_data["background"]["type"]](
                *stored_data["background"]["params"],
                variances=stored_data["background"]["variances"]
            )
            fitted = peaks[0].get_variances()[0] is not None
        else:
            peaks = []
            bg = None
            fitted = False
        known_peaks = shared_data.known_peaks
        isotopes = [kp.get_ele for kp in known_peaks]
        boron_ROI = "B-11" in isotopes and energies[0] < 477 < energies[-1]
        if "peakPairs" in stored_data:
            peakPairs = originalPeakPairs = [
                (peaks[i], known_peaks[j]) for i, j in stored_data["peakPairs"]
            ]
        else:
            peakPairs = originalPeakPairs = None
        return ROI(energies, (energies[0], energies[-1]), cps, shared_data, shared_data.known_peaks,
                   shared_data.region_peaks, user_prefs, indicies, peaks, bg, peakPairs, originalPeakPairs, fitted, boron_ROI)

    def export_to_dict(self):
        """Exports the current ROI state to a dictionary"""
        PIR = [p.export_to_dict() for p in self.region_peaks]
        exportKnownPeaks = [kp.export_to_dict() for kp in self.known_peaks]
        outDict = {
            "indicies": self.indicies,
            "known_peaks": exportKnownPeaks,
            "region_peaks": PIR,
        }

        if self.bg != None:
            outDict["peaks"] = [
                {
                    "type": p.get_type(),
                    "params": p.get_original_params(),
                    "variances": p.get_original_variances(),
                }
                for p in self.peaks
            ]
            outDict["background"] = {
                "type": self.bg.get_type(),
                "params": self.bg.get_original_params(),
                "variances": self.bg.get_original_variances(),
            }
        if self.peakPairs != None:
            outDict["peakPairs"] = [
                (self.peaks.index(p), self.known_peaks.index(kp))
                for p, kp in self.originalPeakPairs
            ]
        return outDict

    def add_peaks(self):
        """Find and add peaks to own model (guesss params)"""
        if self.boron_ROI:
            BPeak = som["peaks"][self.user_prefs["Boron Peak Type"]]
            self.peaks = [BPeak.guess_params(self.energies, self.cps)]
            scrubbedCPS = BPeak.remove_from_data(self.energies, self.cps)
            self.peaks += som["peaks"][self.user_prefs["Peak Type"]].guess_params(
                self.energies, scrubbedCPS
            )
        else:
            self.peaks = som["peaks"][self.user_prefs["Peak Type"]].guess_params(
                self.energies, self.cps
            )

    def add_bg(self):
        """Find and add background to own model (guesss params)"""
        self.bg = som["backgrounds"][self.user_prefs["Background Type"]].guess_params(
            self.energies, self.cps
        )

    def fit(self, reanalyze=False):
        """Fit our model to the data within the ROI, using the guessed params as initial ones"""
        f = lambda x, *params: multiple_peak_and_background(
            self.peaks, self.bg, x, params
        )
        p0 = np.array(
            self.bg.get_params()
            + list(itertools.chain.from_iterable([p.get_params() for p in self.peaks]))
        )
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

    def get_fitted_curve(self, xdata=None):
        """Get the output of our fit (ydata) given x values"""
        if xdata == None:
            xdata = np.arange(self.range[0], self.range[-1], 0.01)
        return [
            list(xdata),
            get_curve(self.peaks, self.bg, xdata),
            list(self.bg.get_ydata(xdata)),
        ]

    def get_closest_peak(self, peak):
        if len(self.peaks) == 0:
            return None
        target = peak.get_ctr()
        distance = [abs(p.get_ctr() - target) for p in self.peaks]
        index = np.argmin(distance)
        return self.peaks[index]

    def set_original_peak_pairs(self, energy_pairs):
        """Original peak pairs are set so that they don't change when the ROI is reanalyzed on new data and can be exported/imported easily."""
        pairs = []
        peakCtr = np.array([p.get_ctr() for p in self.peaks])
        knownCtr = np.array([p.get_ctr() for p in self.known_peaks])

        def nearest_peak(target):
            return self.peaks[np.argmin(abs(peakCtr - target))]

        def nearest_known(target):
            return self.known_peaks[np.argmin(abs(knownCtr - target))]
        pairs = [(nearest_peak(p), nearest_known(kp)) for p, kp in energy_pairs]
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

    # Getters and Setters
    def set_peaks(self, peaks):
        self.peaks = peaks

    def get_peaks(self):
        return self.peaks

    def get_isotopes(self):
        return [kp.get_ele() for kp in self.known_peaks]

    def get_peak_ctrs(self):
        return [p.get_ctr() for p in self.peaks]

    def get_known_peaks(self):
        return self.known_peaks

    def get_range(self):
        return list(self.range)

    def get_formatted_range(self):
        return [
            str(round(float(self.range[0]), decimals=1)),
            str(round(float(self.range[1]), decimals=1)),
        ]

    def set_range(self, new_range):
        self.range = new_range
        self.energies = np.arange(new_range[0], new_range[1], 0.01)

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
        self.known_peaks = peaks
        self.region_peaks = otherPeaks

    def set_background(self, bg):
        self.bg = bg

    def get_background(self):
        return self.bg

    def get_peak_pairs(self):
        return self.peakPairs

    def get_indicies(self):
        return self.indicies
