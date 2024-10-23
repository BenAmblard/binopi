import numpy as np
import scipy as sc
import oivis.oivisfit.fitfunctions as fitfunctions
import binopi.opiparams as opiparams
import binopi.opiutils as opiutils
import emcee

#############################
##### Least Squares Fit #####
#############################

def residuals(theta, model: opiparams.Model, keys, fitdata, refwave = 2.15e-6, cp_fit = False):

    for i in range(len(keys)):
        if model.type[keys[i]] == float:
            model.params[keys[i]] = theta[i] 
        elif model.type[keys[i]] == list:
            model.params[keys[i]].append(theta[i])
            model.params[keys[i]].pop(0)
        else:
            print('Wrong type')
            return

    ucoord = fitdata.v2data['UCOORD']
    vcoord = fitdata.v2data['VCOORD']
    v2wavetable = fitdata.v2data['EFF_WAVE']

    x = [ucoord, vcoord, v2wavetable, refwave]

    llvisibility = fitfunctions.combine_functions(x, model.params)
    
    llvis2model = np.abs(llvisibility['VIS'])**2    
    llv2datapts = fitdata.v2data['VIS2DATA']
    llv2residuals = llv2datapts - llvis2model
    resid = list(llv2residuals)

    if cp_fit:

        u1coord = fitdata.cpdata['U1COORD']
        v1coord = fitdata.cpdata['V1COORD']
        u2coord = fitdata.cpdata['U2COORD']
        v2coord = fitdata.cpdata['V2COORD']
        cpwavetable = fitdata.cpdata['EFF_WAVE']

        y = [u1coord, v1coord, u2coord, v2coord, cpwavetable, refwave]

        llcpmodel = fitfunctions.compute_closure_phase(y, "combine_functions", model)
        llcpdatapts = fitdata.cpdata['T3PHI']
        llcpresiduals = llcpdatapts - llcpmodel
        resid = resid + list(llcpresiduals)

    return resid


def fit_leastSquares(model: opiparams.Model, data, refwave = 2.15e-6, cp_fit = False): 

    ubound = []
    lbound = []
    x0 = []
    keys = []
    for key in model.free.keys():

        if model.free[key]:        
            if model.type[key] == float:
                min = model.limits[key][0]
                max = model.limits[key][1]
                lbound.append(min)
                ubound.append(max)
                x0.append(opiutils.randInterval(min,max))
                keys.append(key)
            elif model.type[key] == list:
                x0 += model.params[key]
                lbound += model.limits[key][0]
                ubound += model.limits[key][1]
                keys += [key]*len(model.params[key])
            else:
                print('This parameter is either of the wrong type or cannot be fitted')
            bounds = (lbound, ubound)


    def fun(theta):
        return residuals(theta, model, keys, data, refwave = refwave, cp_fit = cp_fit)
    
    lsresult = sc.optimize.least_squares(fun, x0, bounds=bounds, method='trf')

    x = lsresult['x']
    cost = lsresult['cost']

    return x, cost, fun(x)


def master_leastSquares(model: opiparams.Model, data, n_it, refwave = 2.15e-6, cp_fit = False):
    freeparams = [i for i in model.free.keys() if model.free[i]]
    x, cost, resid = fit_leastSquares(model, data, refwave = refwave, cp_fit = cp_fit)
    for i in range(n_it-1):
        new_x, new_cost, new_resid = fit_leastSquares(model, data, refwave = refwave, cp_fit = cp_fit)
        if new_cost < cost:
            x = new_x
            cost = new_cost
    print("LS done")
    for i in range(len(freeparams)):
        if model.type[freeparams[i]] == float:
            model.params[freeparams[i]] = x[i]
        elif model.type[freeparams[i]] == list:
            model.params[freeparams[i]].append(x[i])
            model.params[freeparams[i]].pop(0)

    if cp_fit:
        dataerr = np.hstack((data.v2data['VIS2ERR'], data.cpdata['T3PHIERR']))
    else:
        dataerr = data.v2data['VIS2ERR']

    chisquared = 0
    for i in range(len(resid)):
        chisquared += (resid[i]/dataerr[i])**2
    redchisquared = chisquared/(len(resid))


    return model, redchisquared

        
####################
##### MCMC Fit #####
####################

def log_prior(theta, model: opiparams.Model, keys):
    #keys = np.unique(keys)
    for i in range(len(keys)):
        #print(keys[i])
        if model.free[keys[i]]:
            if model.type[keys[i]] == float:
                if not(model.limits[keys[i]][0] <= theta[i] <= model.limits[keys[i]][1]):
                    return -np.inf
            elif model.type[keys[i]] == list:
                for j in range(len(model.params[keys[i]])): #???
                    if not(model.limits[keys[i]][0][j] <= theta[i] <= model.limits[keys[i]][1][j]):
                        return -np.inf
            else:
                print('Wrong type')
                return  
    if not( (-10 < theta[-1] < 1) and (-10 < theta[-2] < 1) ):
        return -np.inf
    return 0.0


def log_likelihood(theta, model: opiparams.Model, keys, fitdata, refwave = 2.15e-6, cp_fit = True):

    for i in range(len(keys)):
        if model.type[keys[i]] == float:
            model.params[keys[i]] = theta[i] 
        elif model.type[keys[i]] == list:
            model.params[keys[i]].append(theta[i])
            model.params[keys[i]].pop(0)
        else:
            print('Wrong type')
            return
    log_f,log_g = theta[-2:]

    ucoord = fitdata.v2data['UCOORD']
    vcoord = fitdata.v2data['VCOORD']
    v2wavetable = fitdata.v2data['EFF_WAVE']
    v2bp = fitdata.v2data['BP']

    u1coord = fitdata.cpdata['U1COORD']
    v1coord = fitdata.cpdata['V1COORD']
    u2coord = fitdata.cpdata['U2COORD']
    v2coord = fitdata.cpdata['V2COORD']
    cpwavetable = fitdata.cpdata['EFF_WAVE']
    cpbp = fitdata.cpdata['BPGEOM']


    x = [ucoord, vcoord, v2wavetable, refwave]
    y = [u1coord, v1coord, u2coord, v2coord, cpwavetable, refwave] 

    llvisibility = fitfunctions.combine_functions(x, model.params)
    
    llvis2model = np.abs(llvisibility['VIS'])**2    
    llcpmodel = fitfunctions.compute_closure_phase(y, "combine_functions", model.params)


    llv2dataerr = fitdata.v2data['VIS2ERR']
    llv2datapts = fitdata.v2data['VIS2DATA']

    llcpdataerr = fitdata.cpdata['T3PHIERR']
    llcpdatapts = fitdata.cpdata['T3PHI']
    
    llv2residuals = llv2datapts - llvis2model
    llcpresiduals = llcpdatapts - llcpmodel

    llv2sigma2 = llv2dataerr**2 + llvis2model**2 * np.exp(2 * log_f)
    llcpsigma2 = llcpdataerr**2 + np.exp(2 * log_g)

    llv2 = -0.5 * (np.sum((llv2residuals) ** 2 / llv2sigma2 + np.log(llv2sigma2)))
    llcp = -0.5 * (np.sum((llcpresiduals) ** 2 / llcpsigma2 + np.log(llcpsigma2)))
    return llv2 + cp_fit*llcp


def log_probability(theta, model: opiparams.Model, keys, fitdata, refwave = 2.15e-6, cp_fit = False):
    
    lp = log_prior(theta, model, keys)
    if not np.isfinite(lp):
        return -np.inf
    proba = lp + log_likelihood(theta, model, keys, fitdata, refwave, cp_fit)
    #print(proba)
    return proba


def fit_mcmc(nsteps, nwalkers, model: opiparams.Model, fitdata, refwave = 2.15e-6, cp_fit = False):
    freekeys = [key for key in model.free.keys() if model.free[key]]
    bounds = ([],[])
    x0 = []
    keys = []
    for key in freekeys:
            print(key)
            if model.type[key] == float:
                min = model.limits[key][0]
                max = model.limits[key][1]
                bounds[0].append(min)
                bounds[1].append(max)
                x0.append(opiutils.randInterval(min,max))
                keys.append(key)
            elif model.type[key] == list:  
                x0 += model.params[key]
                keys += [key]*len(model.params[key])
            else:
                print('This parameter is either of the wrong type or cannot be fitted')
                return
    #print(keys)

    def fun_lp(theta):
        return log_probability(theta, model, keys, fitdata, refwave = refwave, cp_fit = cp_fit)
    
    def randpos(model):
        pos = []
        for key in freekeys:
                if model.type[key] == float:
                    pos.append(opiutils.randGauss(model.limits[key][0], model.limits[key][1], model.params[key], spread = 10))
                if model.type[key] == list:
                    pos += [opiutils.randGauss(model.limits[key][0][i], model.limits[key][1][i], model.params[key][i], spread = 10) for i in range(len(model.params[key]))]
        return pos + [opiutils.randInterval(-10, 0), opiutils.randInterval(-10, 0)]

    pos = np.array([tuple(randpos(model))  for _ in range(nwalkers)])

    nwalkers, ndim = pos.shape
    sampler = emcee.EnsembleSampler(nwalkers, ndim, fun_lp)
    mcmc_result = sampler.run_mcmc(pos, nsteps, progress=True, **{'skip_initial_state_check':True}); 

    return sampler, mcmc_result, model, keys