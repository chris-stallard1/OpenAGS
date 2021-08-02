from abc import ABC, abstractmethod

class Model(ABC):
    """Defines some standard functions for peaks/background models to rpevent code reuse"""
    def __init__(self):
        super().__init__()
    
    @staticmethod
    @abstractmethod
    def get_entry_fields(self):
        """Returns names of fields to be used when manually adding this peak/background model. 
        
        Example: ["Center (keV)", "Amplitude (cps)"] for Gaussian-type peaks.
        """
        pass

    @abstractmethod
    def get_type(self):
        """Returns the key to this peak/background class in the String-Object Map (constants.som)"""
        pass

    @abstractmethod
    def get_num_params(self):
        """Returns the number of parameters this model takes"""
        pass
    
    @abstractmethod
    def get_ydata(self, xdata):
        """Using the current model parameters, model the y values at every point in xdata, and return them"""
        pass
    @abstractmethod
    def get_ydata_with_params(self,xdata,params):
        """Using the provided parameters, model the y values at every point in xdata, and return them (used by the fitter)"""
        pass
    
    @abstractmethod
    def handle_entry(self, entry):
        """Handle user data entry. 
        
        The entry list is in the same order as the one provided by get_entry_fields()
        """
        pass
    @abstractmethod
    def __str__(self):
        pass
    
    #Getters and setters

    @abstractmethod
    def set_params(self, newParams):
        pass
    
    @abstractmethod
    def get_params(self):
        pass

    @abstractmethod
    def set_variances(self, variances):
        pass
    
    @abstractmethod
    def get_variances(self):
        pass

    @abstractmethod 
    def set_original_params(self, params):
        pass

    @abstractmethod
    def get_original_params(self):
        pass

    @abstractmethod
    def set_original_variances(self, vars):
        pass

    @abstractmethod
    def get_original_variances(self):
        pass

class Peak(Model):
    def __init__(self):
        super().__init__()
    
    @abstractmethod
    def get_area(self):
        """Returns the area under this peak"""
        pass
    
    @abstractmethod
    def get_area_stdev(self):
        """Returns the standard deviation of the area under this peak"""
        pass

    @abstractmethod
    def get_fwhm(self):
        """Returns a formatted FWHM width for the peak, or "N/A" if the peak doesn't have a FWHM"""

    #Getters and Setters

    @abstractmethod
    def get_ctr(self):
        pass
    

class StandardPeak(Peak):
    def __init__(self):
        super().__init__()

    @staticmethod
    @abstractmethod
    def guess_params(xdata, ydata):
        """Returns a list of StandardPeak objects representing all found peaks within a set of x-y data represented by xdata and ydata"""
        pass

class BoronPeak(Peak):
    def __init__(self):
        super().__init__()
        
    @staticmethod
    @abstractmethod
    def guess_params(xdata, ydata):
        """Returns a single BoronPeak object representing the closest match for a boron peak within the provided x-y data"""
        pass

    @staticmethod
    @abstractmethod
    def remove_from_data(xdata, ydata):
        """Returns cleaned ydata with the influence fo the boron peak subtracted from it.
        
        This is used because the Boron peak can often overshadow other peaks and make it harder for the peak finding function to work.
        Therefore, we subtract the boron peak from the data and find peaks in the resultant data.
        """
        pass
    

class Background(Model):
    def __init__(self):
        super().__init__()
    
    @staticmethod
    @abstractmethod
    def guess_params(xdata,ydata):
        """Returns a single Background object representing the closest match for a background within the provided x-y data"""
        pass


class Evaluator(ABC):
    """A class which takes ROIs containing peak matches and converts those peak matches to readable results."""
    def __init__(self):
        super().__init__()
    
    @staticmethod
    @abstractmethod
    def get_name():
        "Returns a name for the Evaluator"
        pass

    @staticmethod
    @abstractmethod
    def get_headings(self):
        """Return a list of column headings to be used with the results"""
        pass
    
    @abstractmethod
    def get_results(self):
        """Return a list of results from each ROI"""
        pass