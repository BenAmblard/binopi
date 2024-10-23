import numpy as np
from scipy import constants

mas = constants.arcsecond/1000

class Param():
    def __init__(self, name, min, max, isfree = False, value = None, type = float):
        self.name = name
        self.min = min
        self.max = max
        if value == None:
            if type == float:
                self.value = (self.max + self.min)/2
            if type == list:
                self.value = []
                for i in range(len(self.min)):
                    print(i)
                    self.value.append((self.max[i] + self.min[i])/2)
            if type == int:
                print("Specify a value for type int parameters")
                raise KeyError
        else:
            self.value = value
        self.type = type
        self.isfree = isfree

class Model():
    def __init__(self):
        self.collection = []
        self.params = {}
        self.limits = {}
        self.free = {}
        self.type = {}

    def addParams(self, param: Param):
        self.params[param.name] = param.value
        self.limits[param.name] = (param.min, param.max)
        self.free[param.name] = param.isfree
        self.type[param.name] = param.type
        self.collection.append(param)

    def modifyParams(self, paramname: str, value = None, isfree = None):
        return

def mergeModel(model1: Model, model2: Model):
    "Merges model2 into model1. This modifies model1."
    for param in model2.collection:
        model1.addParams(param)
    return model1

class ModulatedRing(Model):
    def __init__(self, id: int):
        super().__init__()
        self.id = 'vis_modulated_ring;{}:'.format(id)
        self.addParams(Param(self.id+'POSANG', 0, 2*np.pi))
        self.addParams(Param(self.id+'FLUX', 0, 1)) 
        self.addParams(Param(self.id+'OFX', -1*mas, 1*mas))
        self.addParams(Param(self.id+'OFY', -1*mas, 1*mas)) 
        self.addParams(Param(self.id+'RADIUS', 0.0, 5*mas)) 
        self.addParams(Param(self.id+'RINGRADIUS', 0.0, 5*mas)) 
        self.addParams(Param(self.id+'ERATIO', 1, 5)) 
        self.addParams(Param(self.id+'FLOR', 0, 1)) 
        self.addParams(Param(self.id+'MOD', 0, 1, value = 1, type = int)) 
        self.addParams(Param(self.id+'CM',[-1], [1], value = [0.3], type = list)) 
        self.addParams(Param(self.id+'SM', [-1], [1], value = [-0.3], type = list)) 
        self.addParams(Param(self.id+'SPI', -4, 1)) 

class PointSource(Model):
    def __init__(self, id: int):
        super().__init__()
        self.id = 'vis_point_source;{}:'.format(id)
        self.addParams(Param(self.id+'FLUX', 0, 1))
        self.addParams(Param(self.id+'OFX', 0, 1))
        self.addParams(Param(self.id+'OFY', 0, 1))
        self.addParams(Param(self.id+'SPI', -4, 1))

class Ellipsoid(Model):
    def __init__(self, id: int):
        super().__init__()
        self.id = 'vis_ellipsoid;{}:'.format(id)
        self.addParams(Param(self.id+'POSANG', 0, 2*np.pi))
        self.addParams(Param(self.id+'FLUX', 0, 1)) 
        self.addParams(Param(self.id+'OFX', -1*mas, 1*mas))
        self.addParams(Param(self.id+'OFY', -1*mas, 1*mas)) 
        self.addParams(Param(self.id+'RADIUS', 0.0, 5*mas)) 
        self.addParams(Param(self.id+'SPI', -4, 1))
        self.addParams(Param(self.id+'FLOR', 0, 1))

class Gaussian(Model):
    def __init__(self, id: int):
        super().__init__()
        self.id = 'vis_gaussian;{}:'.format(id)
        self.addParams(Param(self.id+'POSANG', 0, 2*np.pi))
        self.addParams(Param(self.id+'FLUX', 0, 1)) 
        self.addParams(Param(self.id+'OFX', -1*mas, 1*mas))
        self.addParams(Param(self.id+'OFY', -1*mas, 1*mas)) 
        self.addParams(Param(self.id+'FWHM', 0.0, 5*mas)) 
        self.addParams(Param(self.id+'SPI', -4, 1))

class UniformRing(Model):
    def __init__(self, id: int):
        super().__init__()
        self.id = 'vis_uniform_ring;{}:'.format(id)
        self.addParams(Param(self.id+'POSANG', 0, 2*np.pi))
        self.addParams(Param(self.id+'FLUX', 0, 1)) 
        self.addParams(Param(self.id+'OFX', -1*mas, 1*mas))
        self.addParams(Param(self.id+'OFY', -1*mas, 1*mas)) 
        self.addParams(Param(self.id+'INDIAM', 0.0, 5*mas)) 
        self.addParams(Param(self.id+'WIDTH', 0.0, 5*mas)) 
        self.addParams(Param(self.id+'ERATIO', 1, 5)) 
        self.addParams(Param(self.id+'SPI', -4, 1)) 
    
