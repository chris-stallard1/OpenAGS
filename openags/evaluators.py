import math
from sigfig import round
from openags.baseClasses import Evaluator

class MassSensEval(Evaluator): 
    """Standard Mass/Sensitivity Evaluator"""              
    def __init__(self, ROIs, relative_to=None):
        self.ROIs = ROIs
        self.relative_to = relative_to
    @staticmethod
    def get_name():
        return "Mass/Sensitivity Results"
    @staticmethod
    def get_headings(ROI):
        output = ROI.get_peak_pairs()[0][1].get_output()
        if relative_to is not None:
            output += f"(Relative to {relative_to.get_ele()} @ {round(relative_to.get_ctr(), decimals=2)} keV)"
        if output == "Peak Area (cps)":
            return ["Isotope", "Peak Centroid (keV)", "FWHM Width (keV)"] + [output, output + " St. Dev"]
        return ["Isotope", "Peak Centroid (keV)", "Peak Area (cps)", "FWHM Width (keV)"] + [output, output + " St. Dev"]

    def get_results(self):
        results = []
        if self.relative_to is not None:
            base_results = None
            for r in self.ROIs:
                for p in r.peak_pairs():
                    if p[1].get_ele() == self.relative_to.get_ele() and abs(p[1].get_ctr() - self.relative_to.get_ctr()) < 1:
                        base_results = p[1].get_results(p[0].get_area(), p[0].get_area_stdev())
        else:
            base_results = None
        for r in self.ROIs:
            for p in r.get_peak_pairs():
                peak_results = p[1].get_results(p[0].get_area(), p[0].get_area_stdev())
                output = p[1].get_output()
                try:
                    if output == "Peak Area (cps)":
                        if base_results:
                            total_unc = 2 * (float(peak_results[1]) ** 2 + float(base_results[1]) ** 2) ** .5
                            formatted_results = [round(float(peak_results[0]/base_results[0]), uncertainty = total_unc, sep=list)[0], round(total_unc, sigfigs = 3)]
                        else:
                            formatted_results = [round(float(peak_results[0]), uncertainty = 2 * float(peak_results[1]), sep=list)[0], round(float(peak_results[1]), sigfigs = 3)]
                        results.append([p[1].get_ele(), round(float(p[0].get_ctr()), decimals=2), p[0].get_fwhm(), *formatted_results])
                    else:
                        if base_results:
                            total_unc = 2 * (float(peak_results[1]) ** 2 + float(base_results[1]) ** 2) ** .5
                            formatted_results = [round(float(peak_results[0]), uncertainty = total_unc, sep=list)[0], round(total_unc, sigfigs = 3)]
                        else:
                            formatted_results = [round(float(peak_results[0]), uncertainty = 2 * float(peak_results[1]), sep=list)[0], round(float(peak_results[1]), sigfigs = 3)]
                        results.append([p[1].get_ele(), round(float(p[0].get_ctr()), decimals=2), round(float(p[0].get_area()), decimals=2), p[0].get_fwhm(), *formatted_results])
                except Exception as e:
                    print(f"Warning: {e}")
                    results.append([p[1].get_ele(), round(float(p[0].get_ctr()), decimals=2), "Error", "Error", "Error", "Error"])
        return results