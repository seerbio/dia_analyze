"""
compare.py


Script to compare two DIA data sets.

"""




def main():


    # -----------------------------------------------------------------------------------
    # Get timestamp and start timing

    import time
    import datetime
    from crozes_base.helpers import helper_functions


    # instantiate the timing and RAM use tracker
    Tracker = helper_functions.Time_and_RAM_Tracker()

    processing_instance_time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')




    # ==================================================================================================================
    #   COMMAND LINE CONFIGURATION AND INPUT
    # ==================================================================================================================

    import argparse

    parser = argparse.ArgumentParser(description='')

    # -----------------------------------------------------------------
    # positional (required) arguments - string assumed for all unless specified

    parser.add_argument('comparison_input', help='folder or list DiaData class pkl files semicolon delimited.')
    parser.add_argument('outpath', help='Output path for comparison results')


    # -----------------------------------------------------------------
    # optional arguments - generic

    parser.add_argument('--verbose', '-v',
                        type=int,
                        choices=[0, 1, 2, 3],
                        default=0,
                        help='Set 0, 1, 2, 3 to control level of verbosity. Default is 0.')

    parser.add_argument('--logfile', '-l',
                        default=None,
                        help='Specify a non-default file name and path string.')

    parser.add_argument('--num_threads', '-t',
                        type=int,
                        default=-1,
                        help='Number of threads to use. Omit or set -1 to automatically determine.')

    parser.add_argument('--analysis_to_file', '-f',
                        dest='analysis_to_file',
                        action='store_true',
                        help='Set the flag to capture output to a file instead of console stream.')
    parser.set_defaults(analysis_to_file=False)

    parser.add_argument('--timestamp_folder', '-ts',
                        dest='timestamp_folder',
                        action='store_true',
                        help='Set the flag to create a timestamped subfolder within specified output folder.')
    parser.set_defaults(timestamp_folder=False)

    parser.add_argument('--pyexe', '-p',
                        type=str,
                        default='python',
                        help='Python 3 executable name or full path if something other than the defaul "python"')

    parser.add_argument('--plots_to_file', '-pl_file',
                        dest='plots_to_file',
                        action='store_true',
                        help='Set flag to write plots to files.')
    parser.set_defaults(plots_to_file=False)

    parser.add_argument('--show_plots', '-pl_show',
                        dest='show_plots',
                        action='store_true',
                        help='Set flag to show plots.')
    parser.set_defaults(show_plots=False)

    parser.add_argument('--plot_height_inches', '-pl_h',
                        type=float,
                        default=4.0,
                        help='Set plot height in inches. Default is 4 inches.')

    parser.add_argument('--plot_width_inches', '-pl_w',
                        type=float,
                        default=6.0,
                        help='Set plot width in inches. Default is 6 inches.')


    # -----------------------------------------------------------------
    # optional arguments - domain-related

    parser.add_argument('--mode', '-m',
                        type=str,
                        choices=['normal'],
                        default=None,
                        help='Set mode. Choices: normal')

    parser.add_argument('--run_characterize', '-char',
                        dest='run_characterize',
                        action='store_true',
                        help='Set the flag to run characterization.')
    parser.set_defaults(run_characterize=False)

    parser.add_argument('--run_v_run_plots', '-rvr',
                        dest='run_v_run_plots',
                        action='store_true',
                        help='Set the flag to do same file comparisons - ex. intensities set A vs. set B of same file.')
    parser.set_defaults(run_v_run_plots=False)

    parser.add_argument('--aligned_analyte_correls', '-aac',
                        dest='aligned_analyte_correls',
                        action='store_true',
                        help='Set the flag to plot correlations of analyte in set A vs. set B aligning files.')
    parser.set_defaults(aligned_analyte_correls=False)



    # Python std and installed package imports
    import os.path

    # Get argument list
    args = parser.parse_args()

    # Positional required arguments

    # if input contains semi-colons, treat as a list of pkl file paths
    if ';' in args.comparison_input:
        comparison_input = [os.path.abspath(x) for x in args.comparison_input.split(';')]
    # otherwise assume it's a single path that contains pkl files
    else:
        comparison_input = os.path.abspath(args.comparison_input)    # convert to abspath since passed to other modules

    outpath = os.path.abspath(args.outpath)     # convert to abspath since passed to other modules


    # Optional arguments - generic
    verbose = args.verbose
    num_threads = args.num_threads
    analysis_to_file = args.analysis_to_file
    timestamp_folder = args.timestamp_folder
    plots_to_file = args.plots_to_file
    show_plots = args.show_plots
    plot_height_inches = args.plot_height_inches
    plot_width_inches = args.plot_width_inches

    # Optional arguments - domain-related
    mode = args.mode
    run_characterize = args.run_characterize
    run_v_run_plots = args.run_v_run_plots
    aligned_analyte_correls = args.aligned_analyte_correls


    if args.logfile is not None:
        logfile = args.logfile
    else:
        logfile = os.path.join(outpath, 'LOGFILE.txt')

    # ------------------------------------------------
    # Make subfolder in requested outpath and make sure output exists

    if timestamp_folder is True:
        outpath = os.path.join(outpath, ('comp_' + processing_instance_time_stamp))

    if not os.path.exists(outpath):
        os.makedirs(outpath)


    # ------------------------------------------------
    # Settings not cmd line args yet



    Tracker.capture('argument parsing and parametrization')




    # ==================================================================================================================
    #   INITIATE LOGGING FOR PROCESSING PART OF PROGRAM
    # ==================================================================================================================

    # logging based on accepted solution here: http://stackoverflow.com/questions/1508467/how-to-log-my-traceback-error
    # for all logging levels, see: https://docs.python.org/3/library/logging.html#logging-levels

    # set up logging to file
    import logging

    # change level to any of: DEBUG > INFO > WARNING > ERROR > CRITICAL  (most to least logging)
    if verbose < 2:
        logging_level = logging.WARNING
    else:
        logging_level = logging.DEBUG

    logging.basicConfig(filename=logfile, level=logging_level)

    # run the whole program inside a try/except block for to catch traceback in log file
    try:


        # ==============================================================================================================
        #   IMPORTS
        # ==============================================================================================================

        # Core Python imports
        import sys
        import numpy as np
        import pandas as pd

        # Imports for multithreading
        from multiprocessing.dummy import Pool as ThreadPool
        import subprocess
        import psutil

        # python built in imports
        from contextlib import redirect_stdout

        # crozes imports not already imported
        from crozes_base.io import boiler_plates

        # imports from this project
        from dia_analyze import dia_comparison

        # Less common 3rd party imports
        from tabulate import tabulate


        Tracker.capture('logging setup and imports')




        # ==============================================================================================================
        #   TRANSLATIONS AND INITIALIZATIONS
        # ==============================================================================================================

        # -----------------------------------------------------------------------------------
        # Redirect stdout to an analysis file if one is specified - otherwise goes to console

        if analysis_to_file is True:
            analysis_file = os.path.join(outpath, 'SCRIPT_OUTPUT.txt')
            sys.stdout = open(analysis_file, 'w')

        boiler_plates.section_header(title='Compare DIA Data Sets', level=0)


        # -----------------------------------------------------------------------------------
        # Auto-determine number of cores if needed

        print()
        if num_threads == -1:
            print('Number of threads automatically detected')
            num_threads = psutil.cpu_count()
        else:
            print('Number of threads used as specified by script parameters')
        print('Number of threads used:', num_threads)
        print()

        Tracker.capture('translations and initializations')




        # ==============================================================================================================
        #       PREPARATION: VALIDATIONS AND CHECKING PARAMETERS
        # ==============================================================================================================

        if verbose > 0:
            boiler_plates.section_header(title='Comparison Input Parameters', level=2)
            print('Comparison input:')
            if isinstance(comparison_input, list):
                for x in comparison_input:
                    print('              ', x)
            else:
                print('              ', comparison_input)
            print()
            print('Output path:  ', outpath)
            print()
            print()

        Tracker.capture('validations and checking parameters')




        # ==============================================================================================================
        #       MAIN PROCESSING SEQUENCE
        # ==============================================================================================================


        # -----------------------------------------------------
        # Instantiate DiaCompare

        comp = dia_comparison.DiaComparison(comparison_input,
                                                outpath=outpath,
                                                plots_to_file=plots_to_file, show_plots=show_plots,
                                                plot_height_inches=plot_height_inches,
                                                plot_width_inches=plot_width_inches,
                                                verbose=verbose)



        # -----------------------------------------------------
        # Characterize results

        # these plots have cmd line flags to produce them as they aren't always necessary and are slower to make
        comp.compare_aligned_samples(run_v_run_plots=run_v_run_plots,
                                     aligned_analyte_correls=aligned_analyte_correls,
                                     )

        comp.compare_aligned_analyte_aggregates()

        comp.compare_unaligned_analytes_aggregates()

        comp.summarize()

        #comp.describe()

        #comp.output(output_file)

        #comp.save(output_file)

        comp.profile()

        Tracker.capture('main processing sequence')



    # =================================================================================================================
    # end of try/except block with whole program run inside to catch traceback in log file
    except:
        logging.exception('Exceptions on main handler:')
        raise

    logging.shutdown()
    # If there's nothing in the logfile, don't write it (or remove it after written actually)
    if os.stat(logfile).st_size == 0:
        os.remove(logfile)

    # =================================================================================================================
    # OUTPUT SEGMENT TIMING DATA
    # =================================================================================================================

    # Write out the overall tracking report
    Tracker.report(verbose=verbose)



if __name__ == "__main__":
    main()

