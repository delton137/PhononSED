'''
Code to fit Lorentzian functions and find the lifetimes.

Daniel C. Elton, 2017

License: MIT
'''
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

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

# -----------------------------------------------------------------------------
# --------------- main code --------------------------------------------------
# -----------------------------------------------------------------------------

#peak_freqs = np.loadtxt('MgOtest1x1x1_frequencies.dat')
#data = np.loadtxt('MgOtest1x1x1_NVT_1_SED.dat')
header = 'MgOtest_super'
k = 2
peak_freqs = np.loadtxt(header+'_'+str(k)+'_frequencies.dat')
data = np.loadtxt(header+'_'+str(k)+'_SED.dat')

num_modes = data.shape[1]-1 #number of modes, dropping one column for time data
num_freqs = data.shape[0]
print("read in", num_modes, " modes at (including any acoustic) ", num_freqs, "frequency points")

freqs = data[:,0]
mode_data = data[:, 1:]

#print("GULP freqs from file:")
#print(peak_freqs)

allparams = np.zeros([4, num_modes])
lifetimes = np.zeros([num_modes])

sw = 5 #search width on each side for fitting, in 1/cm
freq_step = freqs[5]-freqs[4]
iw = int(np.floor(sw/freq_step)) #indexwidth

pw = 50 #plottings width on each side for fitting, in 1/cm
npts = 500 #npts for plotting fit

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

    #this is mostly for handling acoustic modes (w0 ~ 0.0 )
    if (w0 < sw):
        w0 = sw + 1

    bounds = [(max_height/10, 10*max_height), (w0 - sw, w0 + sw), (.001, 10 ), (0, 0)]

    params = fit_function(freqs_2_fit, Y_2_fit, Lorentzian, params, bounds, verbose=False)

    allparams[:, m] = params
    lifetimes[m] = (1/(params[2]*2.99*1e10))/(1e-9)  #lifetimes in ps

fit_peak_freqs = allparams[1,:]


#%%----- plotting ------------------------------------------------------------
num_modes_plot = 12
start_plot = 0
num_plot_windows_to_do = 4 #int(np.ceil((num_modes-start_plot)/num_modes_plot))

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
from scipy.optimize import curve_fit

def scaling_fn(w, A=10e7):
    return A*1./(w**2)

def scaling_fn_arb(w, A=10e7, B=2.0):
    return A*1./(w**B)

A_fit = curve_fit(scaling_fn, peak_freqs, lifetimes) #p0=
x_fit = np.linspace(min(peak_freqs), max(peak_freqs),100)
y_fit = scaling_fn(x_fit, A=A_fit[0])

plt.figure(2)
plt.plot(fit_peak_freqs, lifetimes, '*', label ='')
plt.plot(x_fit, y_fit, '-', label=r'$\omega^{-2}$ fit')
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$\omega$ (cm$^{-1}$)")
plt.ylabel(r"lifetime (ps)")
plt.show(block=True)
