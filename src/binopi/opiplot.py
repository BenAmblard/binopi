import oivis.oifits.oifits_read as pif
import oivis.oifits.oifits_transform as oitrans
import oivis.oivisfit.fitfunctions as fitfunctions
import oivis.oivisfit.fitutils as fitutils
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
import corner
import copy
from fpdf import FPDF
from math import log10, floor
from scipy.special import j0


def find_exp(number) -> int:
    neg = 1
    if abs(number) < 1:
        neg = -1
    base10 = log10(abs(number))
    return neg*abs(floor(base10))

### Global Parameters ### 
max_spfreq = 70 #Mlambda
max_baseline = 200 #m

##########################
### PLOTTING FUNCTIONS ###
##########################
#TODO better control over coloring of the plots (wavelengths, time, etc)

def plot_vis2(oifiles, figVis, scale = 'm', type = 'colors', telescopetype = '', axisScale = ''):
    """ Plots the squared visibility of all imported files.
    
    Parameters:
    -oifiles: OI Data in the Oifits class format

    -figVis: target plot or subplot

    -scale [optional]: Type of scaling for the x axis, possible values are 'm' 
    for baselines in meters or 'freq' for spatial frequency in Mlambda.
    Defaults to 'm'.

    -type [optional]: choose whether to plot with 'colors' based on wavelength or with 'errorbars'.
    Defaults to 'colors'. 
    
    Returns:
    None"""
    
    maxX = max_spfreq
    vis2data = oifiles.v2data
    wavelengths = vis2data['EFF_WAVE']
    baselines = np.sqrt(vis2data['UCOORD']**2  + vis2data['VCOORD']**2)
    if scale == 'freq':
        baselines = baselines/wavelengths/1e6
        figVis.set_xlabel('Spatial Frequency (M'r'$\lambda$'')')
        figVis.set_title('Squared Visibility vs Spatial Frequency')
    elif scale == 'm':
        maxX = max_baseline
        figVis.set_xlabel('Baseline (m)')
        figVis.set_title('Squared Visibility vs Baseline length')
    else:
        print('Not a scale format')

    ## Plot measurements done with different instruments sperately ##

    if len(np.unique(vis2data['INSNAME']))>1:
        cmaps = ['plasma', 'plasma_r']

        for instrument in np.unique(vis2data['INSNAME'])[::-1]:
            index = np.where(np.unique(vis2data['INSNAME'])[::-1] == instrument)[0][0]

            v2 = [vis2data['VIS2DATA'][i] for i in range(len(vis2data['VIS2DATA'])) if vis2data['INSNAME'][i] == instrument]
            v2err = [vis2data['VIS2ERR'][i] for i in range(len(vis2data['VIS2ERR'])) if vis2data['INSNAME'][i] == instrument]
            bls = [baselines[i] for i in range(len(baselines)) if vis2data['INSNAME'][i] == instrument]
            wave = [wavelengths[i] for i in range(len(wavelengths)) if vis2data['INSNAME'][i] == instrument]
            figVis.set_ylabel('Squared Visibility')
            
            if type == 'colors':
                im = figVis.scatter(bls, v2,  c = wave, cmap = cmaps[index], marker = '.', alpha=0.5, label= instrument)
            elif type == 'errorbars':
                im = figVis.scatter(bls, v2, marker = '.', alpha=0.5, label= instrument)
                figVis.errorbar(bls, v2, yerr = v2err, fmt = '.')
            else:
                raise Exception("Invalid plot type. Make sure it is either 'colors' or 'errorbars'.")
            figVis.legend()
    else:

        times = vis2data['MJD']


        v2 = vis2data['VIS2DATA']
        v2err = vis2data['VIS2ERR']

        figVis.set_ylabel('Squared Visibility')

        if type == 'colors':
            im = figVis.scatter(baselines, v2,  c = wavelengths, cmap = 'rainbow', marker = '.', alpha=0.5)
        elif type == 'errorbars':
            im = figVis.scatter(baselines, v2, marker = '.', alpha=0.5, label = telescopetype)
            figVis.errorbar(baselines, v2, yerr = v2err, fmt = '.', alpha=0.5,)
        else:
            raise Exception("Invalid plot type. Make sure it is either 'colors' or 'errorbars'.")

    if axisScale == 'fixed':
        figVis.set_xlim(0, maxX)
        figVis.set_ylim(0,1)
    figVis.grid(True)

    return im
    
def plot_t3phi(oifiles, figT3, scale = 'm', type = 'colors', telescopetype = '', axisScale = ''):
    """ Plots the closure phases of all imported files.
    
    Parameters:
    -oifiles: OI Data in the Oifits class format

    -figT3: target plot or subplot

    -scale [optional]: Type of scaling for the x axis, possible values are 'm' 
    for baselines in meters or 'freq' for spatial frequency in Mlambda.
    Defaults to 'm'.

    -type [optional]: choose whether to plot with 'colors' based on wavelength or with 'errorbars'.
    Defaults to 'colors'. 
    
    Returns:
    None
    """
    maxX = max_spfreq
    t3data = oifiles.cpdata
    wavelengths = t3data['EFF_WAVE'] * 1e6
    baselines = t3data['BPGEOM']
    if scale == 'freq':
        baselines = baselines/wavelengths
        figT3.set_xlabel('Spatial Frequency (M'r'$\lambda$'')')
        figT3.set_title('Closure Phase vs Spatial Frequency')
    elif scale == 'm':
        maxX = max_baseline
        figT3.set_xlabel('Baseline (m)')
        figT3.set_title('Closure Phase vs Baseline length')
    else:
        print('Not a scale format')

    t3 = t3data['T3PHI']
    t3err = t3data['T3PHIERR']

    ## Plot measurements done with different instruments sperately ##

    if len(np.unique(t3data['INSNAME']))>1:
        cmaps = ['plasma', 'plasma_r']
        for instrument in np.unique(t3data['INSNAME'])[::-1]:
            index = np.where(np.unique(t3data['INSNAME'])[::-1] == instrument)[0][0]
            t3 = [t3data['T3PHI'][i] for i in range(len(t3data['T3PHI'])) if t3data['INSNAME'][i] == instrument]
            t3err = [t3data['T3PHIERR'][i] for i in range(len(t3data['T3PHIERR'])) if t3data['INSNAME'][i] == instrument]
            bls = [baselines[i] for i in range(len(baselines)) if t3data['INSNAME'][i] == instrument]
            wave = [wavelengths[i] for i in range(len(wavelengths)) if t3data['INSNAME'][i] == instrument]
            figT3.set_ylabel('Closure Phase (°)')

            if type == 'colors':
                im = figT3.scatter(bls, t3, c = wave, cmap = cmaps[index], marker = '.', alpha=0.5, label= instrument)
            elif type == 'errorbars':
                im = figT3.scatter(bls, t3, marker = '.', alpha=0.5, label= instrument)
                figT3.errorbar(bls, t3, yerr = t3err, fmt = '.', alpha = 0.5)
            else:
                raise Exception("Invalid plot type. Make sure it is either 'colors' or 'errorbars'.")
            figT3.legend()

    else:

        times = t3data['MJD']

        t3 = t3data['T3PHI']
        t3err = t3data['T3PHIERR']

        figT3.set_ylabel('Closure Phase (°)')

        if type == 'colors':
            im = figT3.scatter(baselines, t3, c = wavelengths, cmap = 'rainbow', marker = '.', alpha=0.5)
        elif type == 'errorbars':
            im = figT3.scatter(baselines, t3, marker = '.', alpha=0.5,label = telescopetype)
            figT3.errorbar(baselines, t3, yerr = t3err, fmt = '.')
        else:
            raise Exception("Invalid plot type. Make sure it is either 'colors' or 'errorbars'.")
    
    if axisScale == 'fixed':
        figT3.set_xlim(0, maxX)
        figT3.set_ylim(-180, 180)

    figT3.grid(True)

    return im
    
def plot_spectra(oifiles, figSpec, type = 'colors'):
    """ Plots the spectra of all imported files 
    
    Parameters:
    -oifiles: OI Data in the Oifits class format

    -figSpec: target plot or subplot
     
    -type [optional]: choose whether to plot with 'colors' based on wavelength or with 'errorbars'. 
    Defaults to 'colors'.  
    
    Returns:
    None"""


    v2data = oifiles.v2data
    wavelengths = v2data['EFF_WAVE']*1e6
    v2 = v2data['VIS2DATA']
    v2err = v2data['VIS2ERR']

    if len(np.unique(v2data['INSNAME']))>1:
        cmaps = ['plasma', 'plasma_r']
        for instrument in np.unique(v2data['INSNAME'])[::-1]:
            index = np.where(np.unique(v2data['INSNAME'])[::-1] == instrument)[0][0]
            v2 = [v2data['VIS2DATA'][i] for i in range(len(v2data['VIS2DATA'])) if v2data['INSNAME'][i] == instrument]
            v2err = [v2data['VIS2ERR'][i] for i in range(len(v2data['VIS2ERR'])) if v2data['INSNAME'][i] == instrument]
            wave = [wavelengths[i] for i in range(len(wavelengths)) if v2data['INSNAME'][i] == instrument]
            figSpec.set_ylabel('Closure Phase (°)')

            if type == 'colors':
                figSpec.scatter(wave, v2, c = wave, cmap = cmaps[index], marker = '.', alpha=0.5,label= instrument)
            elif type == 'errorbars':
                figSpec.scatter(wave, v2, marker = '.', alpha=0.5,label= instrument)
                figSpec.errorbar(wave, v2, yerr = v2err, fmt = '.') 
            else:
                raise Exception("Invalid plot type. Make sure it is either 'colors' or 'errorbars'.")
            figSpec.legend()
    else:

        wavelengths = v2data['EFF_WAVE']*1e6
        v2 = v2data['VIS2DATA']
        v2err = v2data['VIS2ERR']

        figSpec.set_ylabel('Closure Phase (°)')
        if type == 'colors':
            figSpec.scatter(wavelengths, v2, c = wave, cmap = 'plasma', marker = '.', alpha=0.5)
        elif type == 'errorbars':
            figSpec.scatter(wavelengths, v2, marker = '.', alpha=0.5)
            figSpec.errorbar(wavelengths, v2, yerr = v2err, fmt = '.')
        else:
            raise Exception("Invalid plot type. Make sure it is either 'colors' or 'errorbars'.")

    figSpec.set_ylabel('Squared Visibility')
    figSpec.set_xlabel('Wavelength ($\mu m$)')
    #figSpec.scatter(wavelengths, v2, marker = '.')
    #figSpec.errorbar(wavelengths, v2, yerr = v2err, fmt = '.')

    figSpec.set_ylim(0,1)
    
    figSpec.set_title('Squared Visibility vs Wavelength')
    figSpec.grid(True)

def plot_uvcov(oifiles, axuv, scale = 'freq', cmap = 'plasma', label = '', setup = True):
    """ Plots the UV coverage of all imported files 
    
    Parameters:
    -oifiles: OI Data in the Oifits class format

    -axuv: target plot or subplot

    -scale [optional]: Type of scaling for the x axis, possible values are 'm' 
    for baselines in meters or 'freq' for spatial frequency in Mlambda.
    Defaults to 'freq'. 

    -cmap [optional]: Color map. If None, plots in solid colors.
    Defaults to 'plasma'
    
    Returns:
    None"""

    v2m = list(filter(lambda x: x['TYPE'] == 'V2', oifiles.oidata))
    based = {'m': 'm', 'freq': '$M\lambda$'}
    ucoord = np.array([x['UCOORD'] for x in v2m])
    vcoord = np.array([x['VCOORD'] for x in v2m])
    wavel = np.array([x['EFF_WAVE'] for x in v2m])
    if (scale == 'freq'):
        maxUV = max_spfreq
        norm = wavel * 1e6
        axuv.set_title('UV Coverage (Spatial Frequency)')
    if (scale == 'm'):
        maxUV = max_baseline
        norm = wavel / wavel
        axuv.set_title('UV Coverage (Baseline length)')
    if str(type(axuv)) == "<class 'matplotlib.projections.polar.PolarAxes'>":
        
        axuv.set_thetagrids(np.array(range(0, 360, 15))[::-1])
        axuv.set_theta_zero_location('N')
        axuv.set_theta_direction(-1)
        axuv.set_theta_offset(np.pi / 2)
        theta = np.arctan2(ucoord,vcoord)
        r = np.sqrt(ucoord ** 2 + vcoord ** 2)
        the = np.hstack((theta, theta + np.pi))
        rad = np.hstack((r / norm, r / norm))
        if cmap == None:
            axuv.scatter(the, rad, marker = '.', alpha=0.5, label = label)
        else:
            c = np.hstack((wavel,wavel))
            axuv.scatter(the, rad, c = c, cmap = cmap, marker = '.', lw=3, alpha=0.5)
        if setup == True:
            axuv.text(np.pi / 2.1, max(rad) * 0.8, 'U (' + based[scale] + ')')
            axuv.text(-np.pi / 20, max(rad) * 0.9, 'V (' + based[scale] + ')')
    else:
        u = np.hstack((ucoord / norm, -ucoord / norm))
        v = np.hstack((vcoord / norm, -vcoord / norm))
        wavel = np.hstack((wavel,wavel))
        if cmap == None:
            axuv.scatter(u, v, marker = '.', alpha=0.5, label = label)
        else:
            axuv.scatter(u, v, c = wavel, cmap = cmap, marker = '.', lw=3, alpha=0.5)
        axuv.set_xlabel('U (' + based[scale] + ')')
        axuv.set_ylabel('V (' + based[scale] + ')')
        axuv.set_ylim(-maxUV, maxUV)
        axuv.set_xlim(-maxUV, maxUV)
        axuv.set_aspect(1.0)
    axuv.grid(True)

def fitReport(savepath, data, mcmc_sampler, model, keys, refwave = 2.15e-6):
    
    def cleanup(l):
    
        n = len(l)
        ave = np.mean(l)
        sigma = np.std(l)
        filterfun = lambda x: ave - 3*sigma <= x <= ave + 3*sigma

        a = [i for i in filter(filterfun, l)]
        if len(a) == len(l):
            return a
        else:
            return cleanup(a)
        
    def get_fit_params(savepath, mcmc_sampler, model, keys):
        ndim = len(keys)
        flat_samples = mcmc_sampler.get_chain(discard=200, thin=15, flat=True)
        fit_params = {}
        with open(savepath+'/params.txt', 'w') as file:
            for i in range(ndim):
                cleansamples = cleanup(flat_samples[:,i])
                mcmc = np.median(cleansamples), np.std(cleansamples), np.std(cleansamples)
                exp = find_exp(mcmc[0])
                if exp != 0:
                    q = mcmc[0]/(10**exp), mcmc[1]/(10**exp), mcmc[2]/(10**exp)
                else:
                    q = mcmc
                txt = r"{3}$= {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}x10^{{{4}}}$"
                txt = txt.format(q[0], q[1], q[2], keys[i], exp)
                file.write(txt+'\n')
                if model.free[keys[i]] == float:
                    fit_params[keys[i]] = q[0], q[1], q[2]
                    model.params[keys[i]] = q[0]
                elif model.free[keys[i]] == list:
                    try:
                        fit_params[keys[i]].append((q[0], q[1], q[2]))
                    except:
                        fit_params[keys[i]] = [(q[0], q[1], q[2])] 
                    model.params[keys[i]].append(q[0])
                    model.params[keys[i]].pop(0)
        return model, fit_params

    def dataplot(data, model, savepath, keys):

        ucoord = data.v2data['UCOORD']
        vcoord = data.v2data['VCOORD']
        v2wavetable = data.v2data['EFF_WAVE']
        v2bp = data.v2data['BP']

        u1coord = data.cpdata['U1COORD']
        v1coord = data.cpdata['V1COORD']
        u2coord = data.cpdata['U2COORD']
        v2coord = data.cpdata['V2COORD']
        cpwavetable = data.cpdata['EFF_WAVE']
        cpbp = data.cpdata['BPGEOM']

        x = [ucoord, vcoord, v2wavetable, refwave]
        y = [u1coord, v1coord, u2coord, v2coord, cpwavetable, refwave] 


        visibility = fitfunctions.combine_functions(x, model.params)

        vis2model = np.abs(visibility['VIS'])**2    
        cpmodel = fitfunctions.compute_closure_phase(y, "combine_functions", model.params)

        dataerr = np.hstack((data.v2data['VIS2ERR'], data.cpdata['T3PHIERR']))
        datapts = np.hstack((data.v2data['VIS2DATA'], data.cpdata['T3PHI']))
        modelpts = np.hstack((vis2model, cpmodel))
        residuals = datapts - modelpts

        chisquared = 0
        for i in range(len(residuals)):
            chisquared += (residuals[i]/dataerr[i])**2
        redchisquared = chisquared/(len(residuals) - len(keys))

        scale = 'freq'
        type = 'errorbars'

        fig_plot = plt.figure()

        axvis2 = plt.subplot(2,1,1)
        axuvcov = plt.subplot(2,2,3, projection = 'polar')
        axt3phi = plt.subplot(2,2,4)

        fig_plot.set_size_inches(16,9)

        plot_uvcov(data, axuvcov, scale)
        plot_vis2(data, axvis2, scale, type, telescopetype = 'Data')
        plot_t3phi(data, axt3phi, scale, type, telescopetype = 'Data')

        axvis2.plot(v2bp / v2wavetable * 1e-6, vis2model, '.r', alpha = 0.5, label = 'Model')
        axvis2.plot(v2bp / v2wavetable * 1e-6, np.abs(data.v2data['VIS2DATA']-vis2model), '.k', alpha = 0.2, label = 'Residuals')
        axt3phi.plot(cpbp / cpwavetable * 1e-6, cpmodel, '.r', alpha=0.5, label = 'Model')
        axt3phi.plot(cpbp / cpwavetable * 1e-6, np.abs(data.cpdata['T3PHI']-cpmodel), '.k', alpha = 0.2, label = 'Residuals')

        fig_plot.suptitle('Data fit results', fontsize = 14)

        axvis2.axhline(linewidth=2, color='gray')

        axvis2.legend(title = "Reduced Chi^2 = {0:.2f}".format(redchisquared))
        axvis2.grid(True)
        axvis2.set_xlim(left = 0)
        axt3phi.legend()
        axt3phi.grid(True)
        plt.tight_layout()
        plt.savefig(savepath+'/mcmcfit.png')

        bmax = 1000 #echantillonnage
        gridim = bmax//20

        pxscale = refwave/(2*bmax)
        maxMas = pxscale * gridim/2

        tmp = {}        
        for key in model.params.keys():
            if key.split(':')[1] == 'FLUX' and key.split(';')[0] == 'vis_point_source':
                tmp[key] = copy.deepcopy(model.params[key])
                model.params[key] = 0
        vgrid = fitutils.vis_grid(bmax, gridim, refwave, model.params, fitfunctions.combine_functions)
        for key in tmp.keys():
            model.params[key] = tmp[key]

        fig1 = plt.figure()

        ax1 = plt.subplot(111)
        ax1.set_xlabel('$\delta$ RA (mas)')
        ax1.set_ylabel('$\delta$ Dec (mas)')
        ax1.set_title('Direct image representation of the best fit')
        ax1.scatter([0], [0], marker='+', color = 'k')
        #ax1.scatter([model.param[key.split(':')[0]+'OFX'] for key in tmp.keys()], [model.param[key.split(':')[0]+'OFY'] for key in tmp.keys()], marker='+', color = 'r')

        im = ax1.imshow((vgrid['IMAGE'])/np.max(vgrid['IMAGE']), cmap = 'jet', aspect='equal', interpolation="bicubic", extent=([maxMas, -maxMas, -maxMas, maxMas]))
        fig1.colorbar(im)
        plt.savefig(savepath+'/image.png')

    def walkersplot(mcmc_sampler, savepath, keys):
        def find_exp(number) -> int:
            neg = 1
            if abs(number) < 1:
                neg = -1
            base10 = log10(abs(number))
            return neg*abs(floor(base10))

        ndim = len(keys)
        flat_samples = mcmc_sampler.get_chain(discard=200, thin=15, flat=True)
        fig_walkers, axes = plt.subplots(ndim, figsize=(10,16), sharex=True)
        samples = mcmc_sampler.get_chain()

        for i in range(ndim):
            cleansamples = cleanup(flat_samples[:,i])
            mcmc = np.median(cleansamples), np.std(cleansamples), np.std(cleansamples)
            exp = find_exp(mcmc[0])
            print(mcmc[0], exp)
            if exp != 0:
                q = mcmc[0]/(10**exp), mcmc[1]/(10**exp), mcmc[2]/(10**exp)
            else:
                q = mcmc
            txt = r"{3}$= {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}x10^{{{4}}}$"
            txt = txt.format(q[0], q[1], q[2], keys[i], exp)
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.plot([mcmc[0] for _ in range(len(samples[:, :, i]))], 'r', label = txt)

            ax.set_xlim(0, len(samples))
            ax.set_ylabel(keys[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
            ax.legend()

        fig_walkers.suptitle('MCMC Fit results', y = 1.0)

        plt.tight_layout()
        plt.savefig(savepath+'/walkers.png')

    def cornerplot(mcmc_sampler, labels, fig): #WIP
        flat_samples = mcmc_sampler.get_chain(discard=1000, thin=15, flat=True)
        corner.corner(
            flat_samples, labels=labels, fig=fig
        );

    class PDF(FPDF):
        def imagex(self, images, x, y, w):
            self.set_xy(x, y)
            self.image(images, w=w )

        def titles(self, title):
            x,y = 6.0,200.0
            self.set_font('Arial', 'B', 8)
            self.cutcell(title, x, y)

        def cutcell(self, txt, x, y, sep = '\n'):
            texts = txt.split('\n')
            for text in texts:
                ntext = text.replace('$','')
                self.set_xy(x,y)
                self.cell(w=100.0, h=5.0, txt=ntext, border=0)
                y += 5.0

    def create_page(pdf, savepath):
        title = "Parameter Values"

        pdf.add_page()
        os.path.split(savepath)
        text = title + '\n'
        with open(savepath+'/params.txt') as file:
            for line in file:
                text += line
        text = r'{}'.format(text)
        pdf.titles(text)
        pdf.imagex(savepath+'/mcmcfit.png', 6, 6, w = 190)
        pdf.imagex(savepath+'/image.png', 6, 126, w = 80)
        pdf.imagex(savepath+'/walkers.png', 106, 126, w = 90)

    model, fit_params = get_fit_params(savepath, mcmc_sampler, model, keys)   

    dataplot(data, model, savepath, keys)
    
    walkersplot(mcmc_sampler, savepath, keys)
    
    #cornerplot(mcmc_sampler, keys)

    pdf = PDF()
    create_page(pdf, savepath)
    pdf.output('ajustements.pdf','F')
    

def abaqPlot(data):
    fig = plt.figure(figsize=(16,9))
    axvis2 = plt.subplot(2,1,1)
    axuvcov = plt.subplot(2,2,3, projection = 'polar')
    axt3phi = plt.subplot(2,2,4)
    plot_vis2(data, axvis2, scale = 'freq', type = 'errorbars')
    plot_t3phi(data, axt3phi, scale = 'freq', type = 'errorbars')
    plot_uvcov(data, axuvcov, scale = 'freq')

    sizes = np.arange(1,10,1)
    kernels = np.array([0.2,0.5, 1, 2])*constants.arcsec/1000
    xabaq = np.arange(0,axvis2.get_xlim()[1], 0.1)
    for size in sizes:
        yabaq = j0(2*np.pi*(xabaq*1e6)*(size*constants.arcsec/1000))#*np.exp(-2 * np.pi * kernels[2] / np.sqrt(3) * (xabaq*1e6))
        yabaq = yabaq**2
        axvis2.plot(xabaq, yabaq, label = 'Ring radius = {} mas'.format(size))


    plt.tight_layout()
    axvis2.legend()
    plt.show()