import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

def process_GR(filename, window_size = 7, control_condition = '2u:C', normalise=True, starting_row = 6):

    raw_data = pd.read_csv(filename)

    #extract the actual data and the factors from raw_data in the order of the growth data table
    data = np.array(raw_data.iloc[starting_row:,3:]).astype(float)
    factors = list(raw_data.iloc[starting_row:,0])

    #generate matching x-values based on samples being taken every 30 minutes starting from 0
    x= np.linspace(0,(data.shape[1]-1)/2, num=data.shape[1])

    #define the curve fitting function
    def func(x, a, b):
        return a * np.exp(b * x)

    #process windows of size window_size along the x-axis. For each window, fit data to the curve
    # and record the second parameter
    exponents = []

    #process all samples (rows)
    for row in range(data.shape[0]):
        this_parameter_array = np.zeros([1,data.shape[1]+10], dtype=float)
        window_start = 1
        while window_start < (data.shape[1] - (window_size + 1)):
            xdata = x[window_start:(window_start+window_size)]
            ydata = data[row,window_start:(window_start+window_size)]
            popt, pcov = curve_fit(func, xdata, ydata)
            this_parameter_array[:,window_start] = popt[1]
            window_start = window_start + 1

        exponents.append(this_parameter_array.max())

    #normalise exponents to control
    if normalise == True:
        exponents = np.array(exponents)
        controls = [i for i, j in enumerate(factors) if j == control_condition]
        exponents = exponents / np.mean(exponents[controls])

    return pd.DataFrame({'Plasmid':factors,'Rel_GR':exponents})
