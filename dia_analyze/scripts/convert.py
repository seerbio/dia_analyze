"""
convert.py


Script to convert DIA libraries or results data between various formats.

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

    parser.add_argument('input_file', help='input file txt or xlsx')
    parser.add_argument('outpath', help='output directory')
    parser.add_argument('name_stem', help='Stem used as file name for libraries and subfolder for results processing.')


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

    parser.add_argument('--input_type', '-reslib',
                        type=str,
                        choices=['results', 'library'],
                        default=None,
                        help='Specify input data as "results" or "library" - required.')

    parser.add_argument('--output_format', '-outfmt',
                        type=str,
                        choices=['std', 'skyline'],
                        default=None,
                        help="Set output format. Choices: results only to 'std' and libraries only to 'skyline'")

    parser.add_argument('--expt_type', '-expt',
                        type=str,
                        choices=['biological', 'replicates', 'dilution_series_species_spike'],
                        default=None,
                        help='Set experiment type. Choices: biological, replicates, dilution_series_species_spike')

    parser.add_argument('--species_column', '-sp',
                        type=str,
                        default=None,
                        help='The header of the column containing species info in the input data.')

    parser.add_argument('--const_species', '-c_sp',
                        type=str,
                        default=None,
                        help='For spiked dilution experiments, the recognition string for constant species')

    parser.add_argument('--var_species', '-v_sp',
                        type=str,
                        default=None,
                        help='For spiked dilution experiments, the recognition string for varied species')

    parser.add_argument('--include_species', '-i_sp',
                        type=str,
                        default=None,
                        help='Include only analytes with species string in comma del list - ex. "_ECOLI,_HUMAN"')

    parser.add_argument('--dil_concs', '-concs',
                        type=str,
                        default=None,
                        help="string,float; pairs mapping for dilution conc ex.  E_0 ,0.0; E_11 ,11.0; E_55 ,55.0")

    parser.add_argument('--max_percent_rms_ldr', '-ldr_rms',
                        type=float,
                        default=20.0,
                        help='The max percent RMS used to determine the linear dynamic range limit.')

    parser.add_argument('--acc_rms_fixed_dil', '-acc_dil',
                        type=float,
                        default=None,
                        help='The relative dilution vs top=1.0 where accuracy and RMS should be quoted. Average of all points used otherwise.')

    parser.add_argument('--exclude_sample_strings', '-xsamp',
                        type=str,
                        default=None,
                        help='Exclude samples containing any string in comma delimited list - ex. "blank,wash"')

    parser.add_argument('--min_trans_per_pg', '-min_trans',
                        type=int,
                        default=1,
                        help='Remove peak groups having fewer than this minimum number of transitions.')

    parser.add_argument('--min_pep_per_prot', '-min_pep',
                        type=int,
                        default=1,
                        help='Remove proteins having fewer than this minimum number of peptides.')

    parser.add_argument('--max_qvalue', '-qval',
                        type=float,
                        default=0.01,
                        help='Filter individual peak groups (not whole rows) with greater than this q-value (FDR).')

    parser.add_argument('--min_score', '-min_sc',
                        type=float,
                        default=None,
                        help='Filter individual peak groups (not whole rows) with less than this score.')


    parser.add_argument('--sample_norm', '-norm',
                        type=str,
                        choices=['ratio_to_analyte_ave', 'none'],
                        default=None,
                        help='Sample normalization method. Choices: ratio_to_analyte_ave, none (default).')

    parser.add_argument('--ignore_provided_norm', '-newnorm',
                        dest='ignore_provided_norm',
                        action='store_true',
                        help='Ignore any sample normalization present in data and recalculate here. False by default.')
    parser.set_defaults(ignore_provided_norm=False)

    parser.add_argument('--ignore_provided_rollup', '-newroll',
                        dest='ignore_provided_rollup',
                        action='store_true',
                        help='Ignore any pep and protein rollup in data and recalculate here. False by default.')
    parser.set_defaults(ignore_provided_rollup=False)

    # this is made an optional argument because can figure it out, but it's faster to read if specified
    parser.add_argument('--delimiter', '-sep',
                        default=None,
                        help='delimiter character in input file. reads faster if specified')

    parser.add_argument('--filter_mods', '-fmod',
                        type=str,
                        choices=['high'],
                        default=None,
                        help="Apply mod filter. Options: 'high' to keep only high prob mods.")

    parser.add_argument('--filter_excluded', '-fexcl',
                        dest='filter_excluded',
                        action='store_true',
                        help='Set the flag to filter excluded items.')
    parser.set_defaults(filter_excluded=False)

    parser.add_argument('--run_characterize', '-char',
                        dest='run_characterize',
                        action='store_true',
                        help='Set the flag to run characterization.')
    parser.set_defaults(run_characterize=False)

    parser.add_argument('--min_ms1', '-min_ms1',
                        type=float,
                        default=None,
                        help='Min precursor m/z value.')

    parser.add_argument('--max_ms1', '-max_ms1',
                        type=float,
                        default=None,
                        help='Max precursor m/z value.')

    parser.add_argument('--min_ms2', '-min_ms2',
                        type=float,
                        default=None,
                        help='Min fragment m/z value.')

    parser.add_argument('--max_ms2', '-max_ms2',
                        type=float,
                        default=None,
                        help='Max fragment m/z value.')

    # Python std and installed package imports
    import os.path

    # Get the package bin path
    program_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # Get argument list
    args = parser.parse_args()

    # Positional required arguments
    input_file = os.path.abspath(args.input_file)       # convert to abspath since passed to other modules
    outpath = os.path.abspath(args.outpath)             # convert to abspath since passed to other modules
    name_stem = args.name_stem

    # Optional arguments - generic
    verbose = args.verbose
    num_threads = args.num_threads
    logfile = args.logfile
    analysis_to_file = args.analysis_to_file
    pyexe = args.pyexe
    plots_to_file = args.plots_to_file
    show_plots = args.show_plots
    plot_height_inches = args.plot_height_inches
    plot_width_inches = args.plot_width_inches


    # Optional arguments - domain-related
    mode = args.mode
    reslib = args.input_type                # really required, but set here to use choices validation
    output_format = args.output_format      # required if library conversion but assumed if results conversion
    delimiter = args.delimiter
    include_species = args.include_species
    exclude_sample_strings = args.exclude_sample_strings

    expt_type = args.expt_type
    species_column = args.species_column
    const_species = args.const_species
    var_species = args.var_species
    dil_concs = args.dil_concs
    max_percent_rms_ldr = args.max_percent_rms_ldr
    acc_rms_fixed_dil = args.acc_rms_fixed_dil

    min_trans_per_pg = args.min_trans_per_pg
    min_pep_per_prot = args.min_pep_per_prot

    max_qvalue = args.max_qvalue
    min_score = args.min_score

    min_ms1 = args.min_ms1
    max_ms1 = args.max_ms1
    min_ms2 = args.min_ms2
    max_ms2 = args.max_ms2

    filter_mods = args.filter_mods
    filter_excluded = args.filter_excluded

    sample_norm = args.sample_norm
    ignore_provided_norm = args.ignore_provided_norm

    ignore_provided_rollup = args.ignore_provided_rollup

    run_characterize = args.run_characterize

    # Get the output from the specified output file (needed if logfile dump or additional optional output)
    if reslib == 'library':
        output_file = os.path.join(outpath, str(name_stem) + '_LIB.txt')
        logfile = output_file.replace('_LIB.txt', '_LOGFILE.txt')
        analysis_file = output_file.replace('_LIB.txt', '_SCRIPT_OUTPUT.txt')
    else:
        outpath = os.path.join(outpath, str(name_stem))
        logfile = os.path.join(outpath, 'LOGFILE.txt')
        analysis_file = os.path.join(outpath, str(name_stem) + '_SCRIPT_OUTPUT.txt')


    # ------------------------------------------------
    # Make subfolder in requested outpath and make sure output exists

    if not os.path.exists(outpath):  # need this here for sake of the log file temporarily written here always
        os.makedirs(outpath)


    # ------------------------------------------------
    # Settings not cmd line args yet



    Tracker.capture('argument parsing and parameterization')




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
        from dia_analyze import dia_data

        # Less common 3rd party imports
        from tabulate import tabulate


        Tracker.capture('logging setup and imports')




        # ==============================================================================================================
        #   TRANSLATIONS AND INITIALIZATIONS
        # ==============================================================================================================

        # -----------------------------------------------------------------------------------
        # Redirect stdout to an analysis file if one is specified - otherwise goes to console

        if analysis_to_file is True:
            sys.stdout = open(analysis_file, 'w')

        boiler_plates.section_header(title='Convert format of DIA library or results', level=0)


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


        # -----------------------------------------------------------------------------------
        # Imply output format

        # the only option for output format for results is the standard format
        if reslib == 'results':
            output_format = 'std'

        # -----------------------------------------------------------------------------------
        # Parse dilution concentrations

        if dil_concs is not None:
            input_list = dil_concs.split(';')

            conc_dict = {}

            for x in input_list:
                pattern, conc = x.split(',')
                conc_dict.update({str(pattern): float(conc)})

            dil_concs = conc_dict

        else:
            dil_concs = {}

        # If no fixed dilution provided to quote acc and rms, set to false
        if acc_rms_fixed_dil is None:
            acc_rms_fixed_dil = False


        # -----------------------------------------------------------------------------------
        # Process sample normalization info

        # default when absent is already None. Here just converting string input to proper None
        if sample_norm == 'none':
            sample_norm = None


        # -----------------------------------------------------------------------------------
        # Process include/exclude info

        if exclude_sample_strings is not None:
            exclude_sample_strings = exclude_sample_strings.split(',')
        else:
            exclude_sample_strings = []

        if include_species is not None:
            include_species = include_species.split(',')
        else:
            include_species = []


        # -----------------------------------------------------------------------------------
        # Infer delimiter from filename (if None, file read is way slower when falling back to Python inferred delim)

        # Use explicit delimiter if provided. Otherwise infer from filename here
        if delimiter is None:
            delimiter, format_name = helper_functions.infer_delimiter_or_format(input_file)
        else:
            format_name = 'unknown'

        Tracker.capture('translations and initializations')




        # ==============================================================================================================
        #       PREPARATION: VALIDATIONS AND CHECKING PARAMETERS
        # ==============================================================================================================

        if verbose > 0:
            boiler_plates.section_header(title='Format Conversion Parameters', level=2)
            print('Analysis name stem:                ', name_stem)
            print()
            print('INPUT:')
            print('    input file:                    ', input_file)
            print('    results or library:            ', reslib)
            print('    delimiter/format:              ', format_name)
            print()
            print('OUTPUT')
            if reslib == 'library':
                print('    library file name:             ', output_file)
            else:
                print('    outpath (sub-folder added):    ', outpath)
            print('    output format:                 ', output_format)
            print()
            print('FILTERING')
            print('    precursor m/z range:           ', min_ms1, ' - ', max_ms1)
            print('    fragment m/z range:            ', min_ms2, ' - ', max_ms2)
            print('    max q-value of a peak group:   ', max_qvalue)
            print('    min score of a peak group:     ', min_score)
            print()
            print('    min trans per peak group:      ', min_trans_per_pg)
            print('    min pep per protein:           ', min_pep_per_prot)
            print()
            print('    species column:                ', species_column)
            print('    include analytes w species:    ', include_species)
            print()
            print('    exclude samples with strings:  ', exclude_sample_strings)
            print()
            print('    filter mods:                   ', filter_mods)
            print('    filter excluded items:         ', filter_excluded)
            print()
            print('PROCESSING')
            print('    experiment type:               ', expt_type)
            if expt_type == 'dilution_series_species_spike':
                print('    constant species:              ', const_species)
                print('    variable (spiked) species:     ', var_species)
                print('    dilution concentration map:    ', dil_concs)
                print('    max %RMS for LDR determination:', max_percent_rms_ldr)
                print('    quote acc and RMS at rel dil:  ', acc_rms_fixed_dil)
            print()
            print('    recompute rollup if present:   ', ignore_provided_rollup)
            print()
            print('    ignore norm if already present:', ignore_provided_norm)
            print('    sample normalization:          ', sample_norm)
            print()
            print('SCRIPT SETTINGS')
            print('    verbose:                       ', verbose)
            print('    num_threads:                   ', num_threads)
            print()
            print()

        Tracker.capture('validations and checking parameters')




        # ==============================================================================================================
        #       MAIN PROCESSING SEQUENCE
        # ==============================================================================================================


        # -----------------------------------------------------
        # Instantiate DiaData class (includes assumed need to convert)

        sw = dia_data.DiaData(input_file, name_stem=name_stem, in_delim=delimiter, output_format=output_format,
                                  outpath=outpath,
                                  ignore_provided_norm=ignore_provided_norm,
                                  ignore_provided_rollup=ignore_provided_rollup,
                                  exclude_sample_strings=exclude_sample_strings,
                                  expt_type=expt_type, species_column=species_column,
                                  const_species=const_species, var_species=var_species,
                                  include_species=include_species, dil_concs=dil_concs,
                                  max_percent_rms_ldr=max_percent_rms_ldr, acc_rms_fixed_dil=acc_rms_fixed_dil,
                                  plot_height_inches=plot_height_inches, plot_width_inches=plot_width_inches,
                                  plots_to_file=plots_to_file, show_plots=show_plots,
                                  verbose=verbose)


        # -----------------------------------------------------
        # Characterize results

        sw.characterize(title='Earliest post-conversion characterization before any filtering',
                        full_analysis=False, plots_to_file=False, outpath=None)


        # -----------------------------------------------------
        # Row filtering

        # remove decoys
        # todo - could make this a config flag to keep decoys, counter to default of removing them
        sw.filter_dia_decoys()

        # this can reduce the set of mods from a very large set to only expected workup mods and very common mods
        if filter_mods is not None:
            sw.filter_mods(filter_type=filter_mods)

        # remove any rows that are marked to be excluded in the provided data
        if filter_excluded is True:
            sw.filter_excluded()

        sw.filter_by_scores(min_score=min_score, max_qvalue=max_qvalue)

        # remove any rows outside m/z ranges
        sw.filter_mz(min_ms1=min_ms1, max_ms1=max_ms1, min_ms2=min_ms2, max_ms2=max_ms2)

        # remove low completeness rows (only 0 for now) - Must be after score/qval filter!
        sw.filter_by_completeness(min_count=1)

        # must be after removal of low completeness
        sw.filter_analytes(min_trans_per_pg=min_trans_per_pg, min_pep_per_prot=min_pep_per_prot)

        # recompute rollup if doing rollup here rather than using existing
        # todo - potentially offer better rollup than sum and make this a cmd arg
        sw.rollup(method='sum', ignore_provided_rollup=ignore_provided_rollup)


        # -----------------------------------------------------
        # Characterize results

        sw.describe()
        #sw.characterize(title='Initial before filtering and normalization')


        # -----------------------------------------------------
        # Further processing

        sw.sample_normalization(method=sample_norm,

                                # settings of method='ratio_to_analyte_ave' (there are more not exposed here):
                                apex_method='mode',     # or 'mode_by_ave_of_top'
                                min_f_complete=0.85,

                                # if provided rollup is NOT ignored, then sample norm applied to it
                                # if rollup is ignored and computed in next call below, inherits norm automatically
                                ignore_provided_rollup=ignore_provided_rollup,

                                ignore_provided_norm=ignore_provided_norm,

                                plots_to_file=plots_to_file)


        # -----------------------------------------------------
        # Characterize results and output

        if sw.lib_or_results == 'library':
            sw.characterize(title='Final after filtering',
                            full_analysis=False, plots_to_file=plots_to_file)
        else:
            sw.characterize(title='Final after filtering and normalization',
                            full_analysis=True, plots_to_file=plots_to_file)


        sw.output(outpath)

        # save a pkl version of the DiaData class, deriving the name from the output name (not for libraries)
        if sw.lib_or_results is 'results':
            pkl_file = os.path.join(outpath, str(name_stem) + '.pkl')
            sw.save(pkl_file)

        sw.profile()


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

