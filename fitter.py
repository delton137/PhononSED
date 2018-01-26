'''
Code to fit Lorentzian functions to spectra and find the lifetimes.
For use with the PhononSED code.

Daniel C. Elton, 2017

License: MIT
'''

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# -------- User-specified inputs ---------------------------------------------
# -----------------------------------------------------------------------------
#header = 'MgOtest_super'

#header = 'RDXtest'
header = 'silicon_test_20modes'

num_modes_plot = 12  # number of modes to plot per plot window
start_plot = 0         # mode to start the plotting at
num_plot_windows_to_do = 1 #int(np.ceil((num_modes-start_plot)/num_modes_plot))

k_list = [1]

sw = 50      #search width on each side for fitting, in 1/cm
pw = 100     #plot's width on each side in 1/cm
npts = 250   #npts for fit curve in plotting

num_k = len(k_list)

# --------------- functions --------------------------------------------------
def Lorentzian(w, params):
    '''
        The Lorentzian function
        arguments:
            params : a list of parametrs with three parameters: [A, w_0, Gamma]
            w : the frequency to evaluate at
        returns:
            the value of the function
    '''
    A = params[0]
    w_0 = params[1]
    Gamma = params[2]
    D = params[3]
    return (A*Gamma/np.pi)/((w_0 - w)**2 + Gamma**2) + D

# -------------------------------------------
def fit_function(dataX, dataY, fit_fn, params, bounds, differential_evolution=False, TNC=True, SLSQP=True, verbose=False):
    '''
    General purpose function for fitting {X, Y} data with a model.
        arguments:
            dataX : Numpy array, X data to fit
            dataY : Numpy array, Y data to fit
            model_fn : the function to fit which is of the form f(x, params)
            params : list of parameters for function
            bounds : list of bounds for the parameters
        returns:
            params : a list of fitted parameters
    '''

    def costfun(params):
        """Wrapper function needed for the optimization method
            Args:
                params: a list of parameters for the model
            Returns:
                The cost (real scalar)
        """
        #diff = (dataY - fit_fn(dataX, params))/dataY
        diff = np.log10(dataY) - np.log10(fit_fn(dataX, params))
        #diff = dataY - fit_fn(dataX, params)

        return np.dot(diff, diff)

    if (differential_evolution == True):
        resultobject = optimize.differential_evolution(costfun, bounds=bounds, maxiter=20000)
        params = resultobject.x
        if (verbose == True): print("diff. evolv. number of iterations = ", resultobject.nit)

    if (TNC == True):
        resultobject = optimize.minimize(costfun, x0=params, bounds=bounds, method='TNC')
        params = resultobject.x
        if (verbose == True): print("TNC number of iterations = ", resultobject.nit)

    if (SLSQP == True):
        resultobject = optimize.minimize(costfun, x0=params, bounds=bounds, method='SLSQP')
        params = resultobject.x
        if (verbose == True): print("SLSQP number of iterations = ", resultobject.nit)

    return params

# -------------------------------------------
def fit_k(freqs_data, mode_data, num_modes):

    allparams = np.zeros([4, num_modes])
    lifetimes = np.zeros([num_modes])

    freq_step = freqs[5]-freqs[4]
    iw = int(np.floor(sw/freq_step)) #indexwidth

    for m in range(0, num_modes):
        max_height = max(mode_data[:,m])
        idx_peak = list(mode_data[:,m]).index(max_height)
        freq_max = idx_peak*freq_step + freq_step


        if (abs(freq_max - peak_freqs[m]) > sw):
            print("WARNING : for mode ", m, " the location of maximum height is not near GULP value!!")
            print(freq_max, "vs", peak_freqs[m])
            idx_peak = peak_freqs[m]/freq_step - 1

        if (idx_peak > len(mode_data[:,1])-1):
            idx_peak = len(mode_data[:,1])-1
            print("WARNING: according to GULP, peak is at higher freq than avail in file")

        if ((idx_peak - iw) < 0):
            idx_peak = iw +  1

        freqs_2_fit = freqs[idx_peak-iw:idx_peak+iw]

        Y_2_fit = mode_data[idx_peak-iw:idx_peak+iw, m]

        w0 = freqs[idx_peak]

        params = [max_height, w0, 1, 0]

        #this is mostly for handling acoustic modes (ie. when w0 ~ 0.0 )
        if (w0 < sw):
            w0 = sw + 1

        bounds = [(max_height/10, 10*max_height), (w0 - sw, w0 + sw), (.001, 10 ), (0, 0)]

        params = fit_function(freqs_2_fit, Y_2_fit, Lorentzian, params, bounds, verbose=False)

        allparams[:, m] = params
        lifetimes[m] = (1/(params[2]*2.99*1e10))/(1e-9)  #lifetimes in ps

    all_fit_peak_freqs[k,:] = allparams[1,:]
    all_lifetimes[k,:] = lifetimes


# --------------- main loop over k values (one file per k) --------------------
for (i, k) in enumerate(k_list):
    gulp_peak_freqs = np.loadtxt(header+'_'+str(k)+'_frequencies.dat')
    data = np.loadtxt(header+'_'+str(k)+'_SED.dat')

    num_modes = data.shape[1]-1 #number of modes, dropping first column since it is the time data
    num_freqs = data.shape[0]
    print("for k=", k, "read in", num_modes, " modes at (including any acoustic) ", num_freqs, "frequency points")

    freqs_data = data[:,0]
    mode_data = data[:, 1:]

    if (i == 0):
        global all_fit_peak_freqs
        global all_lifetimes
        global all_gulp_peak_freqs
        all_gulp_peak_freqs = np.zeros([num_k, num_modes])
        all_fit_peak_freqs = np.zeros([num_k, num_modes])
        all_lifetimes = np.zeros([num_k, num_modes])

    all_gulp_peak_freqs[k,:] = gulp_peak_freqs

    fit_k(freqs_data, mode_data, num_modes)


#-----------------------------------------
for p in range(num_plot_windows_to_do):

    subplot_index = 1

    for m in range(start_plot + p*num_modes_plot, start_plot + (p+1)*num_modes_plot):

        #max_height = max(mode_data[:,m])
        #idx_peak = list(mode_data[:,m]).index(max_height)

        idx_peak = int(peak_freqs[m]/freq_step - 1) #center on GULP frequencies

        if ((idx_peak - iw) < 1):
            idx_peak = iw + 1

        if (idx_peak > len(freqs) - 1):
            idx_peak = len(freqs)-1-pw

        freqs_2_fit = freqs[idx_peak-iw:idx_peak+iw]


        xmin = freqs[idx_peak] - pw #for plotting
        xmax = freqs[idx_peak] + pw
        modelX = np.linspace(xmin, xmax, npts)
        modelXfit = np.linspace(min(freqs_2_fit), max(freqs_2_fit), npts)

        modelY = Lorentzian(modelX, allparams[:, m] )
        modelYfit = Lorentzian(modelXfit, allparams[:, m] )

        Y = mode_data[:, m]

        ax = plt.subplot(np.ceil(float(num_modes_plot)/3.0), 3, subplot_index)
        subplot_index += 1

        plt.plot(freqs, Y, "g", modelX, modelY,"b-", modelXfit, modelYfit,"y-")
        plt.axvline(x=peak_freqs[m], color='k', linestyle='--')
        plt.xlim([xmin, xmax])
        plt.xlabel(r"$\omega$ (cm$^{-1}$)")
        plt.ylabel(r"")
        plt.yscale('log')
        plt.ylim([.1,max([max(Y),max(modelYfit)])])
        ps_label = ("%6.5f" % lifetimes[m])
        plt.text(.55,.8, ps_label+" ps", fontsize = 10, transform=ax.transAxes)

    plt.show(block=True)

#%%------------ fitting and plotting lifetimes vs frequency -------------------
plt.figure(2)
plt.clf()

for k in k_list:
    def scaling_fn(w, A=10e7):
        return A*1./(w**2)

    def scaling_fn_arb(w, A=10e7, B=2.0):
        return A*1./(w**B)

    fit_peak_freqs = all_fit_peak_freqs[k, 3:num_modes]
    lifetimes = all_fit_lifetimes[k, 3:num_modes]

    A_fit = optimize.curve_fit(scaling_fn, fit_peak_freqs, lifetimes) #p0=
    x_fit = np.linspace(min(fit_peak_freqs), max(fit_peak_freqs),100)
    y_fit = scaling_fn(x_fit, A=A_fit[0])

    plt.plot(fit_peak_freqs, lifetimes[3:num_modes], '*', label ='')
    plt.plot(x_fit, y_fit, '-', label=r'$\omega^{-2}$ fit')

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$\omega$ (cm$^{-1}$)")
plt.ylabel(r"lifetime (ps)")
plt.savefig('lifetimes_'+str(k)+'.png')
plt.show(block=True)
