from abc import ABC, abstractmethod

class Peak(ABC):
    def __init__(self):
        super().__init__()
    
    @staticmethod
    @abstractmethod
    def get_ydata_with_params(xdata,params):
        pass

    @staticmethod
    @abstractmethod
    def get_entry_fields(self):
        pass 

    @abstractmethod
    def get_type(self):
        pass
   
    @abstractmethod
    def get_num_params(self):
        pass
    
    @abstractmethod
    def set_params(self, newParams):
        pass
    
    @abstractmethod
    def get_params(self):
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
    @abstractmethod
    def get_ctr(self):
        pass
    
    @abstractmethod
    def set_variances(self, variances):
        pass
    
    @abstractmethod
    def get_area(self):
        pass
    
    @abstractmethod
    def get_area_stdev(self):
        pass
    
    @abstractmethod
    def get_ydata(self, xdata):
        pass

    @abstractmethod
    def handle_entry(self, entry):
        pass

    @abstractmethod
    def to_string(self):
        pass

class StandardPeak(Peak):
    def __init__(self):
        super().__init__()

    @staticmethod
    @abstractmethod
    def guess_params(xdata, ydata):
        pass

class BoronPeak(Peak):
    def __init__(self):
        super().__init__()
        
    @staticmethod
    @abstractmethod
    def guess_params(xdata, ydata):
        pass

    @staticmethod
    @abstractmethod
    def remove_from_data(xdata, ydata):
        pass
    

class Background(ABC):
    def __init__(self):
        super().__init__()
    
    @staticmethod
    @abstractmethod
    def guess_params(xdata,ydata):
        pass
    @staticmethod
    @abstractmethod
    def get_entry_fields(self):
        pass

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def get_num_params(self):
        pass
    
    @abstractmethod
    def set_params(self, newParams):
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
    @abstractmethod
    def get_ydata(self, xdata):
        pass
    @abstractmethod
    def get_ydata_with_params(self,xdata,params):
        pass
    
    @abstractmethod
    def handle_entry(self, entry):
        pass
    @abstractmethod
    def to_string(self):
        pass

class Evaluator(ABC):
    def __init__(self):
        super().__init__()

    @abstractmethod
    def get_headings(self):
        pass
    
    @abstractmethod
    def get_results(self):
        pass