from models import LinearBackground, QuadraticBackground, ArctanBackground, GaussianPeak, KuboSakaiBoronPeak

#Default user preferences, to be changed from settings tab
default_prefs = {
    "peak_type" : "gaussian",
    "boron_peak_type" : "kubo_sakai",
    "background_type" : "linear",
    "overlap_rois" : True,
    "roi_width" : 15,
    "B_roi_width" : 20
}

#String-Object Map: used for UI interactions, makes a lot of things easier & allows for peaks/background types to be interchangeable
som = {
    "backgrounds":
    {
        "linear":LinearBackground,
        "quadratic":QuadraticBackground,
        "arctan" : ArctanBackground
    },
    "peaks":
    {
        "gaussian":GaussianPeak,
        "kubo_sakai":KuboSakaiBoronPeak
    }
}