import sys
import pathlib
sys.path += [str(pathlib.Path().resolve().parent.parent)+'/Modules',str(pathlib.Path().resolve().parent.parent)+'/Notebooks']
from oivis.oifits.oifits_read import Oifits
import numpy as np
from oivis.oivisfit import fitfunctions

def noisymodel(model, exdata, noiselevel, bias = 1, refwave = 2.15e-6):
    """ Adds noise to previously generated synthetic data points and reorganises in the Oifits format

    Parameters:
    -model:

    -noiselevel: relative noise level (type: float)
        example: 1e-2
    
    Returns:
    oifits_read.Oifits object with added noise
    """

    ucoord = exdata.v2data['UCOORD']
    vcoord = exdata.v2data['VCOORD']
    v2wavetable = exdata.v2data['EFF_WAVE']
    v2bp = exdata.v2data['BP']

    u1coord = exdata.cpdata['U1COORD']
    v1coord = exdata.cpdata['V1COORD']
    u2coord = exdata.cpdata['U2COORD']
    v2coord = exdata.cpdata['V2COORD']
    cpwavetable = exdata.cpdata['EFF_WAVE']
    cpbp = exdata.cpdata['BPGEOM']

    x = [ucoord, vcoord, v2wavetable, refwave]
    y = [u1coord, v1coord, u2coord, v2coord, cpwavetable, refwave] 

    llvisibility = fitfunctions.combine_functions(x, model.params)
    
    v2 = np.abs(llvisibility['VIS'])**2
    cp = fitfunctions.compute_closure_phase(y, "combine_functions", model.params)

    noisydata = Oifits()
    numv2 = len(v2)
    numcp = len(cp)

    b1 = np.sqrt(u1coord ** 2 + v1coord ** 2)
    b2 = np.sqrt(u2coord ** 2 + v2coord ** 2)
    b3 = np.sqrt((u1coord + u2coord) ** 2 + (v1coord + v2coord) ** 2)

    for i in range(numv2):
        noisypoint = {}
        noisypoint['INSNAME'] = 'SIMULATION'
        noisypoint['MJD'] = '0'
        noisypoint['TYPE'] = 'V2'
        noisypoint['EFF_WAVE'] = v2wavetable[i]
        noisypoint['WAVE_TABLE'] = v2wavetable
        noisypoint['UCOORD'] = ucoord[i]
        noisypoint['VCOORD'] = vcoord[i]
        if (ucoord[i] != 0):
            noisypoint['THETA'] = np.arctan(vcoord[i] / ucoord[i]) * 180 / np.pi
        else: 
            # if some telescope are ignored (VLTI) the oifits will still contain data at 
            #baselines where theres is no noisypointurements but with - values
            noisypoint['THETA'] = 512
            raise("Some U coordinates = 0 . Missing Baselines ? THETA set to 512 degrees")
        noisypoint['BP'] = v2bp[i]
        #noisypoint['VIS2DATA'] = (model['vis2'][i])*(1+noiselevel*np.random.randn())
        noisypoint['VIS2DATA'] = (v2[i]+noiselevel*np.random.randn())
        noisypoint['VIS2ERR'] = 2*np.abs(noiselevel)
        noisydata.oidata.append(noisypoint)
        
    for i in range(numcp):
        noisypoint = {}
        noisypoint['INSNAME'] = 'SIMULATION'
        noisypoint['MJD'] = '0'
        noisypoint['TYPE'] = 'T3'
        noisypoint['EFF_WAVE'] = cpwavetable[i]
        noisypoint['WAVE_TABLE'] = cpwavetable
        noisypoint['U1COORD'] = u1coord[i]
        noisypoint['V1COORD'] = v1coord[i]
        noisypoint['U2COORD'] = u2coord[i]
        noisypoint['V2COORD'] = v2coord[i]
        noisypoint['T3PHI'] = cp[i] + bias*np.random.randn()
        noisypoint['T3PHIERR'] = bias
        noisypoint["BP1"] = b1[i]
        noisypoint["BP2"] = b2[i]
        noisypoint["BP3"] = b3[i]
        noisypoint['BMAX'] = max(b1[i], b2[i], b3[i])
        noisypoint['BGEO'] = (b1[i] * b2[i] * b3[i]) ** (1. / 3)

        noisydata.oidata.append(noisypoint)

    noisydata.get_cpdata()
    noisydata.get_v2data()
    return noisydata