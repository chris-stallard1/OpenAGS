from openags.models import LinearBackground, QuadraticBackground, ArctanBackground, GaussianPeak, KuboSakaiBoronPeak

#Default user preferences, to be changed from settings tab
default_prefs = {
    "Peak Type" : "Simple Gaussian",
    "Boron Peak Type" : "Physical B-11",
    "Background Type" : "linear",
    "Overlap ROIs" : True,
    "ROI Width (keV)" : 15,
    "B-11 ROI Width (keV)" : 20
}

#String-Object Map: used for UI interactions, makes a lot of things easier & allows for peaks/background types to be interchangeable
som = {
    "backgrounds":
    {
        "Linear":LinearBackground,
        "Quadratic":QuadraticBackground,
        "Arctan" : ArctanBackground
    },
    "peaks":
    {
        "Simple Gaussian":GaussianPeak,
        "Physical B-11":KuboSakaiBoronPeak
    }
}