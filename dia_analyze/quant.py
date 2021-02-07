import os.path
import numpy as np
import pandas as pd
import statistics
from scipy.optimize import curve_fit

from crozes_base.helpers import helper_functions
from crozes_base.io import boiler_plates


"""
From: Clin Biochem Rev. 2008 Aug; 29(Suppl 1): S49–S52.   (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2556583/)

Limit of Blank (LoB), Limit of Detection (LoD), and Limit of Quantitation (LoQ) are terms used to describe the smallest 
concentration of a measurand that can be reliably measured by an analytical procedure.

LoB is the highest apparent analyte concentration expected to be found when replicates of a blank sample containing no 
analyte are tested.

LoB = meanblank + 1.645(SDblank)

LoD is the lowest analyte concentration likely to be reliably distinguished from the LoB and at which detection is 
feasible. LoD is determined by utilising both the measured LoB and test replicates of a sample known to contain a low 
concentration of analyte.

LoD = LoB + 1.645(SD low concentration sample)

LoQ is the lowest concentration at which the analyte can not only be reliably detected but at which some predefined 
goals for bias and imprecision are met. The LoQ may be equivalent to the LoD or it could be at a much higher 
concentration.
"""


def analyze_dilution_series(rel_concs_list, intensities_list, **kwargs):

    verbose = kwargs.get('verbose', 0)
    method = kwargs.get('method', 'linear from top')
    max_percent_rms_ldr = kwargs.get('max_percent_rms_ldr', 20.0)
    show_plots = kwargs.get('show_plots', False)
    acc_rms_fixed_dil = kwargs.get('acc_rms_fixed_dil', False)  # ex pass 0.1 to quote these stats at the 1/10 dil point

    fit_ydata = None
    fit_xdata = None
    linear_ydata = None

    # initialize output
    ldr = np.nan
    ave_acc = np.nan
    ave_rms = np.nan


    # make copies of input since will alter objects here
    intensities = list(intensities_list)
    rel_concs = list(rel_concs_list)
    top_conc = max(rel_concs)

    # input check - must have no nan
    #if np.nan in intensities:  # this didn't work because were getting non-numpy nan values somehow still float type
    # go through in reverse order so we can delete as we go
    for i in range(len(rel_concs) - 1, -1, -1):
        if np.isnan(intensities[i]):
            del intensities[i]
            del rel_concs[i]


    unique_concs = list(set(rel_concs))

    # if we have not points left from the top conc because it was all nan and removed above, then fail (return)
    if top_conc not in unique_concs:
        return

    # sort descending
    unique_concs.sort(reverse=True)

    # get replicate counts at each conc point
    counts = \
        np.array([sum([1 for i in range(len(intensities)) if rel_concs[i] == x]) for x in unique_concs])

    # get averages at each conc point
    ave_ints = \
        np.array(
            [np.mean([intensities[i] for i in range(len(intensities)) if rel_concs[i] == x]) for x in unique_concs])

    # get averages at each conc point
    std_ints = \
        np.array([np.std([intensities[i] for i in range(len(intensities)) if rel_concs[i] == x]) for x in unique_concs])

    # mask out points w n=1
    std_ints = np.ma.masked_where(std_ints == 0, std_ints)

    cv_ints = std_ints / ave_ints
    cv_ints = np.ma.filled(cv_ints, fill_value=np.nan)

    ave_cv = np.ma.masked_invalid(cv_ints).mean()
    ave_cv = np.ma.filled(ave_cv, fill_value=np.nan)



    # slopes of proximal cons - not used now
    """
    slopes = []
    for i in range(1, len(unique_concs)):
        slope = (ave_ints[i - 1] - ave_ints[i]) / (unique_concs[i - 1] - unique_concs[i])
        slopes.append(slope)

    if verbose > 2:
        print('slps', slopes)
    """

    if show_plots is True:
        # create an linear scale x series that is regular in log space for plotting the resulting curve if showing plots
        non_zero_concs = [x for x in rel_concs if x > 0]
        non_zero_min = min(non_zero_concs)
        extend_orders = 2
        step = (np.log10(max(rel_concs)) - (np.log10(non_zero_min) - extend_orders)) / 1000
        regular_log_x = list(np.arange(np.log10(non_zero_min) - extend_orders, np.log10(max(rel_concs)) + step, step))
        fit_xdata = [10 ** x for x in regular_log_x]

        # Get a linear behavior indicator line that extends from the highest point
        lin_slope = ave_ints[0] / unique_concs[0]   # assumes these are descending sorted
        linear_ydata = [lin_slope * x for x in fit_xdata]



    def interpolate_conc_at_critical_rms(yt, x1, x2, y1, y2):

        #print(yt)
        #print(x1, x2)
        #print(y1, y2)

        if (y1 - y2) == 0:
            xt = np.nan
        else:
            xt = x1 - (y1 - yt) * (x1 - x2) / (y1 - y2)

        return xt



    if method == 'linear from top':

        # convert to relative conc as fraction of max
        unique_rel_concs = [x / max(unique_concs) for x in unique_concs]

        max_obs_ave_int = ave_ints[0]
        dil_expected_int = np.array(unique_rel_concs) * max_obs_ave_int

        # mask out zero conc
        masked_dil_expected_int = np.ma.masked_where(dil_expected_int == 0, dil_expected_int)
        percent_error_vs_dil_expected = 100 * np.abs((ave_ints - masked_dil_expected_int)) / masked_dil_expected_int


        # compute RMS error - don't do this in general as assumes a model for true value
        rms = ((cv_ints * 100)**2 + percent_error_vs_dil_expected**2)**0.5
        rms = np.ma.filled(rms, fill_value=np.nan)

        if acc_rms_fixed_dil is False:

            # get the mean accuracy error and rms, excluding both the top point, and invalid points (ex. 0 conc)
            ave_acc = np.ma.masked_invalid(percent_error_vs_dil_expected[1:]).mean()
            ave_rms = np.ma.masked_invalid(rms[1:]).mean()

            ave_acc = np.ma.filled(ave_acc, fill_value=np.nan)
            ave_rms = np.ma.filled(ave_rms, fill_value=np.nan)

        else:

            # get the index of the specified dilution
            if acc_rms_fixed_dil in unique_rel_concs:

                fixed_dil_ind = unique_rel_concs.index(acc_rms_fixed_dil)

                ave_acc = percent_error_vs_dil_expected[fixed_dil_ind]
                ave_rms = rms[fixed_dil_ind]

            else:
                ave_acc = np.nan
                ave_rms = np.nan

        # ----------------------------------
        # Linear dynamic range - SIMPLE FOR NOW - find lowest conc below threshold RMS error

        indices_passing_rel_concs = [i for i in range(len(unique_rel_concs)) if rms[i] < max_percent_rms_ldr]

        # if the first point's percent cv is worse than the max percent rms, then no LDR, not matter what passes
        # this prevents negative interpolated conc and other weird stuff (ex if top point high cv and lower better)
        if cv_ints[0] * 100 > max_percent_rms_ldr:
            rel_conc_at_thresh = 1.0

        elif len(indices_passing_rel_concs) > 0:

            # get the index of the min conc passing the max rms threshold (the max index w descending conc sort)
            i_max_pass = max(indices_passing_rel_concs)

            """
            print('max i', i_max_pass)
            print('len  ', len(unique_rel_concs), type(len(unique_rel_concs)))
            print('sum  ', (i_max_pass + 1), type((i_max_pass + 1)))
            print('conc ', unique_rel_concs)
            print('rms  ', rms)
            """

            # if passing conc at max index is the end of the list, use that point
            if (i_max_pass + 1) == len(unique_rel_concs):
                rel_conc_at_thresh = unique_rel_concs[i_max_pass]

            # if have passing but no real next point, then can't interpolate past, so just use last passing point
            elif not np.isreal(rms[(i_max_pass + 1)]):
                rel_conc_at_thresh = unique_rel_concs[i_max_pass]

            # if have a real non-passing rms error at the next lower conc, try to interpolate
            else:
                rel_conc_at_thresh = interpolate_conc_at_critical_rms(max_percent_rms_ldr,
                                                                      unique_rel_concs[i_max_pass],
                                                                      unique_rel_concs[i_max_pass + 1],
                                                                      rms[i_max_pass],
                                                                      rms[i_max_pass + 1])

        # if no passing values
        elif np.isreal(rms[1]) and np.isreal(cv_ints[0]):
            # interpolate between top and second point (top point always x=1.0, y=the cv for top)
            rel_conc_at_thresh = interpolate_conc_at_critical_rms(max_percent_rms_ldr,
                                                                  unique_rel_concs[0],
                                                                  unique_rel_concs[1],
                                                                  cv_ints[0] * 100, # this is same as rms[0]
                                                                  rms[1])

        else:
            rel_conc_at_thresh = 1.0

        # compute the actual linear dynamic range as ratio of top rel conc to rel conc a threshold
        ldr = 1.0 / rel_conc_at_thresh

        if verbose > 2 or ldr < 1: # == 1.0:
            print()
            print()
            print('----- ANALYZE DIL SERIES --------------')
            print(intensities)
            print(rel_concs)
            print()
            print('counts:         ', counts)
            #print('unique concs:   ', unique_concs)
            print('aves            ', ave_ints)
            print('stds            ', std_ints)
            print('cvs             ', cv_ints)
            #print('ave cv          ', ave_cv)
            print('uniq rel concs:    ', unique_rel_concs)
            print('max ave int:       ', max_obs_ave_int)
            print('dil expected       ', dil_expected_int)
            print('% error v exptd    ', percent_error_vs_dil_expected)
            print('average % err      ', ave_acc)
            print('% rms error        ', rms)
            print('average rms        ', ave_rms)
            print('threshold %rms:    ', max_percent_rms_ldr)
            print('rel conc at thresh:', rel_conc_at_thresh)
            print('lin dyn range:     ', ldr)


    elif 'param logistic' in method:

        # get initial estimates of parameters using guidance via comments below taken from:
        # https://www.mathworks.com/matlabcentral/fileexchange/38043-five-parameters-logistic-regression-there-and-back-again

        # A is the lower asymptote so guess it with min(y)
        a_init = min(intensities)

        # B is the Hill's slope so guess it with the slope of the line between first and last point.
        b_init = 1.0 #max(slopes) # (max(intensities) - min(intensities)) / (max(rel_concs) - min(rel_concs))

        # C is the inflection point (the concentration of analyte where you have half of the max response) so guess it
        # finding the concentration whose response is nearest to the mid response.
        c_init = statistics.median(rel_concs)

        # D is the upper asymptote so guess it with max(y)
        d_init = max(intensities)

        # E is the asymmetric factor and so guess it with no asymmetry (E=1).
        e_init = 1.0


        if method == '5-param logistic':

            init_param_guesses = (a_init, b_init, c_init, d_init, e_init)
            bounds = ((0, 0, 0, 0, 0), (d_init, np.inf, np.inf, np.inf, np.inf))

            popt, pcov = curve_fit(five_param_logistic, rel_concs, intensities,
                                   p0=init_param_guesses, method='trf', bounds=bounds)

            fit_ydata = five_param_logistic(fit_xdata, *popt)


        elif method == '4-param logistic':
            bounds = ((0, 0, 0, 0), (d_init, np.inf, np.inf, np.inf))
            init_param_guesses = (a_init, b_init, c_init, d_init)

            popt, pcov = curve_fit(four_param_logistic, rel_concs, intensities,
                                   p0=init_param_guesses, method='trf', bounds=bounds)

            fit_ydata = four_param_logistic(fit_xdata, *popt)


        elif method == '3-param logistic':
            # b, the Hill slope becomes fixed at 1.0 in the 3 param model - # todo - this doesn't work w mass spec!!

            bounds = ((0, 0, 0), (d_init, np.inf, np.inf))
            init_param_guesses = (a_init, c_init, np.inf)

            popt, pcov = curve_fit(three_param_logistic, rel_concs, intensities,
                                   p0=init_param_guesses, method='trf', bounds=bounds)

            fit_ydata = three_param_logistic(fit_xdata, *popt)

        if verbose > 2:
            print('Fit results')
            print('a - bottom:     ', popt[0])
            if method != '3-param logistic':
                print('b - Hill slope: ', popt[1])
            print('c - inflection: ', popt[2])
            print('d - top:        ', popt[3])
            if method == '5-param logistic':
                print('e - non-sym:    ', popt[4])

    if show_plots is True:
        fit_plot(rel_concs, intensities, method=method,
                 fit_xdata=fit_xdata, fit_ydata=fit_ydata,
                 linear_ydata=linear_ydata,
                 show_plot=True)

    # make sure empty values are np.nan
    if np.isreal(ldr) is False:
        ldr = np.nan
    if np.isreal(ave_acc) is False:
        ave_acc = np.nan
    if np.isreal(ave_rms) is False:
        ave_rms = np.nan

    return ldr, ave_acc, ave_rms




def five_param_logistic(x, a, b, c, d, e):

    return d + (a - d) / ((1 + ((x / c) ** b)) ** e)




def four_param_logistic(x, a, b, c, d):

    #return d + (a - d) / (1 + ((x / c) ** b))  # this is equivalent to the next line
    return a + (x ** b) * (d - a) / (x ** b + c ** b)




def three_param_logistic(x, a, c, d):

    #return d + (a - d) / (1 + (x / c))
    #return a + x * (d - a) / (c + x)
    return a + x * (d - a) / (x + c)



def fit_plot(xdata, ydata, **kwargs):

    method = kwargs.get('method', 'fit')
    fit_ydata = kwargs.get('fit_ydata', None)
    fit_xdata = kwargs.get('fit_xdata', None)
    linear_ydata = kwargs.get('linear_ydata', None)
    outpath = kwargs.get('outpath', None)
    show_plot = kwargs.get('show_plot', False)
    verbose = kwargs.get('verbose', 0)

    import matplotlib.pyplot as plt

    xlabel = 'Relative Concentration'
    ylabel = 'Measured Intensity'
    title = 'Dose-Response Curve Fitting'


    fig = None
    fig_file = None

    if show_plot is True or outpath is not None:

        fig = plt.figure()
        plt.title(title)
        plt.plot(xdata, ydata, 'b-', color='green', marker='o', linestyle='none', label='raw data')

        if fit_ydata is not None:
            plt.plot(fit_xdata, fit_ydata, 'r-', label=str(method))

        if linear_ydata is not None:
            plt.plot(fit_xdata, linear_ydata, 'g-', label='linear from top')

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(loc=2)

        if show_plot is True:
            plt.show()

        """
        if outpath is not None:
            if '.txt' in title:
                file_base_name = title.replace('.txt', '.png')
            else:
                file_base_name = str(title) + '.png'

            file_base_name = helper_functions.convert_string_to_valid_filename(file_base_name)

            fig_file = os.path.join(outpath, file_base_name)

            fig.savefig(fig_file)
        """

        plt.close()

    return



# ======================================================================================================================
#           TESTS
# ======================================================================================================================


if __name__ == '__main__':

    import os.path

    import pandas as pd
    import numpy as np


    # dilution series factors: E4: 1, E3: 0.5, E2: 0.1, E1: 0.02, E0: 0
    concs = [0, 0, 0, 0.02, 0.02, 0.02, 0.1, 0.1, 0.1, 0.5]

    inta = [np.nan, np.nan, np.nan, 16.01135559, 8.888278971, np.nan, 96.73799184, 95.57018793, 99.19045982, 1222.252661]
    intb = [17.46916991, np.nan, np.nan, 186.039034, 131.8240476, 130.3972418, 1040.829194, 885.3690342, 954.1832112, 5734.029965]  # P00634_GIDALPLTGQYTHYALNK
    intc = [56.28627694, np.nan, 16.48446089, 233.2895956, 180.3766452, 183.4496275, 1478.524412, 1374.22859, 1385.160663, 9346.530875]
    intd = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 333.598966, 362.4935127, 402.7918769, 2776.739788]
    int_trans_b = [20.62874693, np.nan, np.nan, 232.232105, 164.9422788, 171.1232227, 1095.237082, 1039.578169, 976.1289448, 5017.293877]

    intensities_lists = [inta, intb, intc, intd, int_trans_b]
    """
    for ints in intensities_lists:

        #ldr, ave_acc, ave_rms = analyze_dilution_series(concs, ints, method='4-param logistic', show_plots=True, verbose=2)

        ldr, ave_acc, ave_rms = \
            analyze_dilution_series(concs, ints, method='linear from top', show_plots=False, verbose=3)
    """

    """ 
    FAIL CASE REAL DATA - 1
    
    [1047.0, 136.0, 433.0, 639.0, 1542.0, 573.0, 1335.0, 683.0, 1269.0, 1634.0, 1337.0, 1048.0, 2391.0, 1007.0, 3292.0]
    [0.0, 0.0, 0.0, 11.0, 11.0, 11.0, 275.0, 275.0, 275.0, 55.0, 55.0, 55.0, 550.0, 550.0, 550.0]
    
    counts:          [3 3 3 3 3]
    unique concs:    [550.0, 275.0, 55.0, 11.0, 0.0]
    aves             [2230.         1095.66666667 1339.66666667  918.          538.66666667]
    stds             [939.76841119 293.04076319 239.24092924 442.05655747 379.34534252]
    cvs              [0.4214208121918922 0.26745430166087863 0.17858243038348953
     0.4815430909231815 0.7042302150675843]
    ave cv           0.41064617004540527
    uniq rel concs:     [1.0, 0.5, 0.1, 0.02, 0.0]
    dil expected        [2230.  1115.   223.    44.6    0. ]
    % error v exptd     [0.0 1.7339312406576912 500.7473841554559 1958.2959641255604 --]
    average % err       820.2590931738914
    % rms error         [42.14208121918922 26.801577422166993 501.06572381569714
     1958.8879295647866 --]
    average rms         828.9184102675503
    rel conc at thresh: -6.458447096185599
    lin dyn range:      -0.15483598225812
    
    
    FAIL CASE REAL DATA - 2
    
    [nan, nan, nan, 5874.0, nan, nan, 87555.0, 96323.0, 88683.0, 31845.0, 28194.0, 21999.0, 157065.0, 191596.0, 140842.0]
    [0.0, 0.0, 0.0, 11.0, 11.0, 11.0, 275.0, 275.0, 275.0, 55.0, 55.0, 55.0, 550.0, 550.0, 550.0]
    
    counts:          [3 3 3 3 3]
    unique concs:    [550.0, 275.0, 55.0, 11.0, 0.0]
    aves             [163167.66666667  90853.66666667  27346.                     nan
                 nan]
    stds             [21164.81410791  3894.72303965  4064.09128834            nan
                nan]
    cvs              [0.12971205962725127 0.04286808868084714 0.1486173951709496 -- --]
    ave cv           0.10706584782634936
    uniq rel concs:     [1.0, 0.5, 0.1, 0.02, 0.0]
    dil expected        [nan nan nan nan nan]
    % error v exptd     [-- -- -- -- --]
    average % err       --
    % rms error         [-- -- -- -- --]
    average rms         --
    rel conc at thresh: 1.0
    lin dyn range:      1.0
    
    
    FAIL CASE REAL DATA - 3
    
    [nan, nan, nan, 3435.0, nan, nan, 35273.0, 50620.0, 43999.0, 11013.0, 13224.0, 8947.0, 72500.0, 79009.0, 64470.0]
    [0.0, 0.0, 0.0, 11.0, 11.0, 11.0, 275.0, 275.0, 275.0, 55.0, 55.0, 55.0, 550.0, 550.0, 550.0]
    
    counts:          [3 3 3 3 3]
    unique concs:    [550.0, 275.0, 55.0, 11.0, 0.0]
    aves             [71993.         43297.33333333 11061.33333333            nan
                nan]
    stds             [5946.33876151 6285.00089278 1746.41238607           nan           nan]
    cvs              [0.08259606852767012 0.1451590758348048 0.15788443702425917 -- --]
    ave cv           0.12854652712891138
    uniq rel concs:     [1.0, 0.5, 0.1, 0.02, 0.0]
    dil expected        [nan nan nan nan nan]
    % error v exptd     [-- -- -- -- --]
    average % err       --
    % rms error         [-- -- -- -- --]
    average rms         --
    rel conc at thresh: 1.0
    lin dyn range:      1.0


    FAIL CASE REAL DATA - 4

    [nan, 167.0, 20.0, 143.0, 17.0, nan, 285.0, nan, 715.0, 123.0, 172.0, nan, 2411.0, 430.0, nan]
    [0.0, 0.0, 0.0, 11.0, 11.0, 11.0, 275.0, 275.0, 275.0, 55.0, 55.0, 55.0, 550.0, 550.0, 550.0]
    
    counts:          [3 3 3 3 3]
    aves             [nan nan nan nan nan]
    stds             [nan nan nan nan nan]
    cvs              [-- -- -- -- --]
    uniq rel concs:     [1.0, 0.5, 0.1, 0.02, 0.0]
    max ave int:        nan
    dil expected        [nan nan nan nan nan]
    % error v exptd     [-- -- -- -- --]
    average % err       --
    % rms error         [-- -- -- -- --]
    average rms         --
    rel conc at thresh: 1.0
    lin dyn range:      1.0
    
    
    FAIL CASE REAL DATA - 5 - e
    
    
    [4.16740973702218, 6.243683245213836, 4.887248075861549, 96.82794053123334, 95.68512930106034, 89.23513791252893, 27.20394194133066, 20.61482583908347, 21.187153009459443, 194.84466840965504, 201.81330966848688, 206.16448392436567]
    [11.0, 11.0, 11.0, 275.0, 275.0, 275.0, 55.0, 55.0, 55.0, 550.0, 550.0, 550.0]
    
    counts:          [3 3 3 3]
    aves             [200.94082067  93.91606925  23.0019736    5.09944702]
    stds             [4.6622944  3.34263798 2.98041307 0.86081326]
    cvs              [0.02320233 0.03559176 0.12957206 0.16880522]
    uniq rel concs:     [1.0, 0.5, 0.1, 0.02]
    max ave int:        200.94082066750252
    dil expected        [200.94082067 100.47041033  20.09408207   4.01881641]
    % error v exptd     [0.0 6.5236531469357875 14.47138276938745 26.889275221084326]
    average % err       14.47138276938745
    % rms error         [ 2.32023259  7.43140515 19.42447174 31.74878166]
    average rms         19.42447174393767
    rel conc at thresh: -4.027521427342473e-05
    lin dyn range:      -24829.166474722937
    
    """


    concs = [0.0, 0.0, 0.0, 11.0, 11.0, 11.0, 275.0, 275.0, 275.0, 55.0, 55.0, 55.0, 550.0, 550.0, 550.0]

    inta = [1047.0, 136.0, 433.0, 639.0, 1542.0, 573.0, 1335.0, 683.0, 1269.0, 1634.0, 1337.0, 1048.0, 2391.0, 1007.0, 3292.0]
    intb = [np.nan, np.nan, np.nan, 5874.0, np.nan, np.nan, 87555.0, 96323.0, 88683.0, 31845.0, 28194.0, 21999.0, 157065.0, 191596.0, 140842.0]
    intc = [np.nan, np.nan, np.nan, 3435.0, np.nan, np.nan, 35273.0, 50620.0, 43999.0, 11013.0, 13224.0, 8947.0, 72500.0, 79009.0, 64470.0]
    intd = [np.nan, 167.0, 20.0, 143.0, 17.0, np.nan, 285.0, np.nan, 715.0, 123.0, 172.0, np.nan, 2411.0, 430.0, np.nan]
    intb2 = [np.nan, np.nan, np.nan, 5874.0, np.nan, np.nan, 87555.0, 96323.0, 88683.0, np.nan, np.nan, np.nan, 157065.0, 191596.0, 140842.0]
    inte = [np.nan, np.nan, np.nan, 4.167409737, 6.243683245, 4.887248076, 96.82794053, 95.6851293, 89.23513791, 27.20394194, 20.61482584, 21.18715301, 194.8446684, 201.8133097, 206.1644839]
    intf = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 9.577721419, 9.264723934, 10.69751468, np.nan, np.nan, np.nan, 16.99880611, 12.47046258, 18.45213695]

    intensities_lists = [intf] #, intb, intc, intd, intb2, inte, intf]

    for ints in intensities_lists:

        ldr, ave_acc, ave_rms = \
            analyze_dilution_series(concs, ints, method='linear from top',
                                    acc_rms_fixed_dil=0.1,
                                    max_percent_rms_ldr=30.0,
                                    show_plots=False, verbose=3)



