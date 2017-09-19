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
        diff = (dataY - fit_fn(dataX, params))/dataY

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


# --------------- main code --------------------------------------------------

peak_freqs = np.loadtxt('MgOtest_frequencies.dat')
data = np.loadtxt('MgOtest_1_SED.dat')

num_modes = 5 #data.shape[1]-3 #number of modes, dropping the 3 acoustic modes
num_freqs = data.shape[0]
print("read in", num_modes, " non-acoustic modes at ", num_freqs, "frequency points")

freqs = data[:,0]

mode_data = data[:, 4:]
peak_freqs = peak_freqs[3:]

allparams = np.zeros([4, num_modes])
lifetimes = np.zeros([num_modes, 1])

sw = 5 #search width on each side for fitting, in 1/cm
freq_step = freqs[4]-freqs[3]
iw = np.floor(sw/freq_step) #indexwidth

pw = 10 #plottings width on each side for fitting, in 1/cm
npts = 500

for m in range(0, num_modes):
    max_height = max(mode_data[:,m])
    idx_max = list(mode_data[:,m]).index(max_height)
    
    if ((idx_max - iw) < 0): 
        iw = idx_max - 1

    freqs_2_fit = freqs[idx_max-iw:idx_max+iw]
    
    print(min(freqs_2_fit), max(freqs_2_fit))
    
    Y_2_fit = mode_data[idx_max-iw:idx_max+iw, m]

    w0 = freqs[idx_max]
    
    params = [max_height, w0, 1, 0]
    if (w0 < sw):
        sw = w0 - 1
    bounds = [(max_height/10, 10*max_height), (w0 - sw, w0 + sw), (.1, 10 ), (0, 0)]
    params = fit_function(freqs_2_fit, Y_2_fit, Lorentzian, params, bounds, verbose=False)
    
    print(params)
    allparams[:, m] = params
    lifetimes[m] = 1/(params[2]*2.99*1e10)
    
print(lifetimes/1e-9)

#%%----- plotting ------------------------------------------------------------
plt.clf()

for m in range(0, num_modes):

    max_height = max(mode_data[:,m])
    idx_max = list(mode_data[:,m]).index(max_height)
    
    if ((idx_max - iw) < 0): 
        iw = idx_max - 1

    freqs_2_fit = freqs[idx_max-iw:idx_max+iw]

    xmin = freqs[idx_max] - pw #for plotting
    xmax = freqs[idx_max] + pw
    modelX = np.linspace(xmin, xmax, npts)
    modelXfit = np.linspace(min(freqs_2_fit), max(freqs_2_fit), npts)
    
    modelY = Lorentzian(modelX, params)
    modelYfit = Lorentzian(modelXfit, params)


    Y = mode_data[:, m]
    
    plt.subplot(np.ceil(float(num_modes)/3.0), 3, m+1)
    plt.plot(freqs, Y, "g", modelX, modelY,"b-", modelXfit, modelYfit,"y-")
    plt.axvline(x=peak_freqs[m], color='k', linestyle='--')
    plt.xlim([xmin, xmax])
    plt.xlabel(r"$\omega$ (cm$^{-1}$)")
    plt.ylabel(r"")
    plt.yscale('log')

plt.show(block=True)
