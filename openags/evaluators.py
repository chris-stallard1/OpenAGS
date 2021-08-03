import math
from sigfig import round
from openags.baseClasses import Evaluator

class HBondAnalysis(Evaluator):
    """Experimental Evaluator for looking at the location of the H Peak relative to some other common ones."""
    def __init__(self, ROIs):
        self.HPeak = None
        self.AlPeak = None
        self.FPeak = None
        for r in ROIs:
            if r.get_range()[0] < 1633 and r.get_range()[1] > 1633:
                pairs = r.get_peak_pairs()
                for p in pairs:
                    if p[1].get_ele() == "F-20":
                        self.FPeak = p[0]
                        break
            
            if r.get_range()[0] < 1778 and r.get_range()[1] > 1778:
                pairs = r.get_peak_pairs()
                for p in pairs:
                    if p[1].get_ele() == "Al-28":
                        self.AlPeak = p[0]
                        break
            
            if r.get_range()[0] < 2224 and r.get_range()[1] > 2224:
                pairs = r.get_peak_pairs()
                for p in pairs:
                    if p[1].get_ele() == "H-2":
                        self.HPeak = p[0]
                        break
    @staticmethod
    def get_headings(_):
        return ["Value", "95% CI +/-", "Comparator Used", "Width Ratio (H/Al)"]
    @staticmethod
    def get_name():
        return "Hydrogen Peak Location Results"
    
    def get_results(self):
        HCtr = self.HPeak.get_ctr()
        HWidth = self.HPeak.get_params()[2]
        HVar = self.HPeak.get_variances()[0]

        AlCtr = self.AlPeak.get_ctr()
        AlWidth = self.AlPeak.get_params()[2]
        AlVar = self.AlPeak.get_variances()[0]
        return [((HCtr - AlCtr) - 444.33)*1000, (2*math.sqrt(HVar+AlVar))*1000, "Al-28", HWidth/AlWidth]

class MassSensEval(Evaluator): 
    """Standard Mass/Sensitivity Evaluator"""              
    def __init__(self, ROIs):
        self.ROIs = ROIs
    @staticmethod
    def get_name():
        return "Mass/Sensitivity Results"
    @staticmethod
    def get_headings(ROI):
        output = ROI.get_peak_pairs()[0][1].get_output()
        if output == "Peak Area (cps)":
            return ["Isotope", "Peak Centroid (keV)", "FWHM Width (keV)"] + [output, output + " St. Dev"]
        return ["Isotope", "Peak Centroid (keV)", "Peak Area (cps)", "FWHM Width (keV)"] + [output, output + " St. Dev"]

    def get_results(self):
        results = []
        for r in self.ROIs:
            for p in r.get_peak_pairs():
                peak_results = p[1].get_results(p[0].get_area(), p[0].get_area_stdev())
                output = p[1].get_output()
                if output == "Peak Area (cps)":
                    formatted_results = [round(float(peak_results[0]), uncertainty = 2 * float(peak_results[1]), sep=list)[0], round(float(peak_results[1]), sigfigs = 3)]
                    results.append([p[1].get_ele(), round(float(p[0].get_ctr()), decimals=2), p[0].get_fwhm(), *formatted_results])
                else:
                    formatted_results = [round(float(peak_results[0]), uncertainty = 2 * float(peak_results[1]), sep=list)[0], round(float(peak_results[1]), sigfigs = 3)]
                    results.append([p[1].get_ele(), round(float(p[0].get_ctr()), decimals=2), round(float(p[0].get_area()), decimals=2), p[0].get_fwhm(), *formatted_results])
        return results