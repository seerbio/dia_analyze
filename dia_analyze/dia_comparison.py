import os.path
from pathlib import Path

import numpy as np
import pandas as pd

import functools

from tabulate import tabulate

from crozes_base.helpers import helper_functions
from crozes_base.io import boiler_plates
from crozes_base.stats import plotting, stat_methods

from dia_analyze import dia_data






class DiaComparison(object):
    """
    A DIA data object comparison. Expect this to mainly be used for results comparison, but it should work for
    qualitative comparison of two libraries.

    """


    def __init__(self, comparison_input, **kwargs):
        """
        Create a new DiaComparison class instance.

        :param comparison_input: a list of files or a path pointer to a folder containing DIA files

        """

        self.verbose = kwargs.get('verbose', 0)

        self.outpath = kwargs.get('outpath', None)
        self.show_plots = kwargs.get('show_plots', False)
        self.plots_to_file = kwargs.get('plots_to_file', False)
        self.tables_to_file = kwargs.get('tables_to_file', False)
        self.plot_height_inches = kwargs.get('plot_height_inches', 4)
        self.plot_width_inches = kwargs.get('plot_width_inches', 6)

        # instantiate the timing and RAM use tracker
        self.tracker = helper_functions.Time_and_RAM_Tracker()

        self.file_list = []
        self.name_stem_list = []
        self.dia_data_list = []
        self.levels_present = []

        # results objects
        self.analyte_alignment_dict = {}    # details of alignment of analytes across all data sets
        self.results = {'unaligned_cum_yields': pd.DataFrame()}                   # the master results collector

        # status trackers
        self.samples_aligned = False
        self.analytes_aligned = False

        self.plot_param_dict = \
            {'cv_int':
                 {'axis_range': (-0.08, 0.6), 'log_ratio_axis_range': (-1.5, 1.5),
                  'critical_values': [0.05, 0.20], 'category': 'precision'},

             'f_complete':
                 {'axis_range': (-0.1, 1.02), 'log_ratio_axis_range': (-1.5, 1.5), 'descending': True,
                  'critical_values': [0.5, 0.9], 'category': 'detection', 'suppress_comp_weighting': True},

             'ave_qval':
                 {'axis_range': (-0.0015, 0.01), 'log_ratio_axis_range': (-6.0, 6.0), 'suppress_levels': {'Transition'},
                  'critical_values': [0.002, 0.008], 'category': 'detection'},

             'ldr':
                 {'axis_range': (-10.0, 52), 'log_ratio_axis_range': (-3, 3), 'descending': True,
                  'critical_values': [5.0, 10.0], 'category': 'accuracy', 'suppress_comp_weighting': True},

             'ave_acc':
                 {'axis_range': (-15, 100), 'log_ratio_axis_range': (-3, 3), 'cum_max': 1000.0,
                  'critical_values': [10, 40], 'category': 'accuracy', 'suppress_comp_weighting': True},

             'ave_rms':
                 {'axis_range': (-15, 100), 'log_ratio_axis_range': (-3, 3), 'cum_max': 1000.0,
                  'critical_values': [10, 40], 'category': 'accuracy', 'suppress_comp_weighting': True},

             'ave_delta_rt':
                 {'axis_range': (-0.3, 1.5), 'log_ratio_axis_range': (-3.0, 3.0), 'suppress_levels': {'Transition'},
                  'critical_values': [0.1, 0.5], 'category': 'rt prediction'}
             }


        self.import_data(comparison_input)

        # get all the sample columns (quant data) in the same order in the self.df table in each result
        self.align_samples()

        # get super-set of all analytes across sets, fill in where missing in each set, sort all to common index
        self.align_analytes()

        self.describe_sets()




    def import_data(self, comparison_input):

        if self.verbose > 0:
            boiler_plates.section_header(title='Read input data sets', level=1)

        # -----------------------------------------------------
        # Sort out what kind of input it is

        # if it's a string, it must be a path.
        if type(comparison_input) == str:

            # make sure the path exists and has contents
            if os.path.exists(comparison_input):
                # get all files in path
                _, self.file_list = helper_functions.get_all_files_in_path(comparison_input,
                                                                           ['.pkl', '.txt', '.csv', '.xlsx'])

        elif type(comparison_input) == list:

            # make sure each item in the list exists
            self.file_list = [x for x in comparison_input if os.path.exists(x)]

            if len(self.file_list) < len(comparison_input):
                print('FATAL ERROR!! Not all pkl files in input list exist. Check input!')
                print()
                print('Input file list')
                for x in comparison_input.file_list:
                    print('  ', x)
                print()
                print('Files that exist:')
                for x in self.file_list:
                    print('  ', x)
                quit()



        # -----------------------------------------------------
        # Make sure have at least two DIA data sets

        if len(self.file_list) < 2:
            print('FATAL ERROR!! A minimum of two data sets is required. File pointers found:', len(self.file_list))
            print('File list')
            for x in self.file_list:
                print('  ', x)


        # -----------------------------------------------------
        # Instantiate DiaData class instances

        for file in self.file_list:

            # if pkl in file assume it's a pickled DiaData object
            if '.pkl' in file:
                sw = helper_functions.read_pickled_object(file)

            else:
                sw = dia_data.DiaData(file, verbose=self.verbose)

                # Make sure in std results format
                if sw.lib_or_results == 'results' and sw.dialect != 'std':
                    sw.convert_format('std')

            file_stem = Path(file).stem

            if self.verbose > 0:
                sw.describe(title=('DIA set read in: ' + file_stem), depth='basic')

            # add the data set to collection
            self.dia_data_list.append(sw)
            self.name_stem_list.append(file_stem)

        self.tracker.capture('read input')

        return




    def align_samples(self):

        if self.verbose > 0:
            boiler_plates.section_header(title='Align columns and common sort', level=1)

        ref_headers = self.dia_data_list[0].intensity_columns

        if self.verbose > 0:
            print('Reference data set for column alignment: ', self.dia_data_list[0].name_stem)

        # track synonyms for speed and renaming - {<observed name>: <name in set one>
        synonym_dict = {}

        # track the success or failure of each alignment attempt
        alignment_statuses = []

        # align all sets against the intensity_columns names in the first set
        for s in range(1, len(self.dia_data_list)):

            sw = self.dia_data_list[s]

            for header in sw.intensity_columns:

                matched = False

                # check for exact match
                if header in ref_headers:
                    if self.verbose > 1:
                        print('EXACT MATCH FOUND:', header)
                    matched = True
                    pass

                # check for previously matched
                elif header in synonym_dict:
                    matched = True

                # look for substring match bidirectionally
                else:

                    header_reduced = header.replace('.wiff', '').replace('.raw', '').strip(' Area').strip(' Intensity')

                    for ref_header in ref_headers:

                        ref_reduced = \
                            ref_header.replace('.wiff', '').replace('.raw', '').strip(' Area').strip(' Intensity')

                        if ref_reduced in header_reduced or header_reduced in ref_reduced:
                            if self.verbose > 1:
                                print('MATCH FOUND:', header, '  to  ', ref_header)
                            matched = True

                            # ------------------------------------
                            # capture the match by adding to the synonym dict

                            # check that we don't already have a match, as this would indicate non-unique matching
                            if header in synonym_dict:
                                print('ERROR!! Match already found for header:', header)
                                print('        which means code is not good enough')
                                print('prior match:   ', synonym_dict[header])
                                print('current match: ', ref_header)
                                quit()

                            synonym_dict.update({header: ref_header})

                            break


                if matched is False:
                    print('ERROR!! NO MATCH found for header:', header)
                    print('        which means code is not good enough')
                    print('list matching against from first data set:')
                    for x in ref_headers:
                        print(x)


                alignment_statuses.append(matched)

            if self.verbose > 0:
                if matched is True:
                    print('All sample columns matched in data set:     ', sw.name_stem)
                else:
                    print('FAILED match of sample columns in data set: ', sw.name_stem)

            if matched is True:

                for level in ['trans', 'pep', 'prot']:

                    if level == 'trans':
                        ldf = sw.df
                    elif level == 'pep':
                        ldf = sw.pep_df
                    else:
                        ldf = sw.prot_df

                    if ldf is not None:

                        # rename the headers to first set names
                        ldf.rename(columns=synonym_dict, inplace=True)

                        # sort the columns by the order in set 1 (this also drops unmatched columns)
                        non_sample_cols = [x for x in ldf.columns if x not in ref_headers]
                        new_headers = non_sample_cols + ref_headers
                        ldf = ldf[new_headers]


                # set the intensity_columns attribute in this set (captures renaming and reordering in each set
                sw.intensity_columns = ref_headers

        # Check if all sets were aligned
        if alignment_statuses.count(False) == 0:
            self.samples_aligned = True
            print('YES!! All samples aligned for all sets')

        self.tracker.capture('align columns')

        return




    def align_analytes(self):

        if self.verbose > 0:
            boiler_plates.section_header(title='Align Analytes', level=1)

        self.levels_present = ['Transition']

        # check if pep_df present in all sets (assumed prot_df is also if so)
        count_pep_dfs = sum([1 for x in self.dia_data_list if x.pep_df is not None])

        # keep this tolerance to absence so can compare libraries (results will ALWAYS have higher levels)
        if count_pep_dfs == len(self.dia_data_list):
            self.levels_present += ['Peptide', 'Protein']


        for level in self.levels_present:

            boiler_plates.section_header(title=(level + ' Level'), level=2)

            # todo - for trans level, first check if any have self.mz_indexing True and align w 'trans_ID_mz' if so

            if level == 'Transition':
                list_of_all_index_sets = [set(sw.df.index.values) for sw in self.dia_data_list]
            elif level == 'Peptide':
                list_of_all_index_sets = [set(sw.pep_df.index.values) for sw in self.dia_data_list]
            else:
                list_of_all_index_sets = [set(sw.prot_df.index.values) for sw in self.dia_data_list]

            # get union of all analytes
            union_index = functools.reduce(lambda a, b: a | b, list_of_all_index_sets)


            # build matching summary table at this level
            total_in_all_sets = len(union_index)
            for s in range(len(list_of_all_index_sets)):
                total_in_set = len(list_of_all_index_sets[s])
                percent_intersection = 100 * total_in_set / total_in_all_sets
                set_name = self.dia_data_list[s].name_stem

                # if first one, add level and properties of the union over all sets
                if level not in self.analyte_alignment_dict:
                    self.analyte_alignment_dict.update({level: [{'set_name': 'TOTAL',
                                                                 'total_in_set': total_in_all_sets,
                                                                 'percent_intersection': 100.0}]})

                # add a dictionary for each data set with its metrics
                self.analyte_alignment_dict[level].append({'set_name': set_name,
                                                           'total_in_set': total_in_set,
                                                           'percent_intersection': percent_intersection})


            if self.verbose > 0:
                for s in range(len(list_of_all_index_sets)):
                    print('Unique index values:  ',
                          format(len(list_of_all_index_sets[s]), ','), '  Data set:', self.dia_data_list[s].name_stem)
                print()
                print('Union of all indices: ', format(len(union_index), ','))
                print()


            # take the union of all IDs, convert to list, and alpha sort it.
            sorted_index = sorted(list(union_index))

            if self.verbose > 0:
                print()
                print('First 5 index values:')
                print()
                for x in sorted_index[:5]:
                    print('  ', x)


            # apply index to all datasets
            for sw in self.dia_data_list:

                if level == 'Transition':
                    sw.df = sw.df.reindex(sorted_index, copy=False)
                elif level == 'Peptide':
                    sw.pep_df = sw.pep_df.reindex(sorted_index, copy=False)
                else:
                    sw.prot_df = sw.prot_df.reindex(sorted_index, copy=False)

        self.analytes_aligned = True

        self.tracker.capture('align rows')

        return




    def describe_sets(self, **kwargs):

        depth = kwargs.get('depth', 'basic')


        def print_level_table(level):

            if level not in self.analyte_alignment_dict:
                return

            headers = ['Data Set', 'Analytes', '% of Union']
            datatable = []
            l_dict = self.analyte_alignment_dict[level]
            for subdict in l_dict:
                row = [subdict['set_name'], format(subdict['total_in_set'], ','),
                       round(subdict['percent_intersection'], 2)]
                datatable.append(row)

            # move total row to end of list
            datatable = datatable[1:] + [datatable[0]]

            boiler_plates.section_header(title=(level + ' Alignment'), level=2)
            print(tabulate(datatable, headers=headers, tablefmt='simple'))

            return

        if self.verbose > 0:
            boiler_plates.section_header(title='Describe Comparison Sets', level=1)

            #for sw in self.dia_data_list:
            #    sw.describe(title=('DIA set: ' + sw.name_stem), depth=depth)

            if self.analytes_aligned is True:
                for level in ['Transition', 'Peptide', 'Protein']:
                    print_level_table(level)

        return




    def compare_aligned_samples(self, **kwargs):

        run_v_run_plots = kwargs.get('run_v_run_plots', False)
        aligned_analyte_correls = kwargs.get('aligned_analyte_correls', True)

        show_plots = kwargs.get('show_plots', self.show_plots)
        plots_to_file = kwargs.get('plots_to_file', self.plots_to_file)
        tables_to_file = kwargs.get('tables_to_file', self.tables_to_file)
        height_inches = kwargs.get('height_inches', self.plot_height_inches)
        width_inches = kwargs.get('width_inches', self.plot_width_inches)

        # cannot do these plots if samples are not aligned
        if self.samples_aligned is False:
            return



        # -------------------------------------------
        # run vs. run plotting

        if run_v_run_plots is True:

            if self.verbose > 0:
                boiler_plates.section_header(title='Run vs. Run Plotting', level=1)

            if plots_to_file is True:
                outpath = os.path.join(self.outpath, 'aligned_samples')
                if not os.path.exists(outpath):
                    os.makedirs(outpath)
            else:
                outpath = None

            if len(self.dia_data_list) > 2:
                print()
                print('WARNING!! Can only plot run vs. run intensity comparisons for two input data sets now and')
                print('         ', len(self.dia_data_list), 'dia data sets were specified. ')
                print('          Only the first two sets will be compared')

            for level in self.levels_present:

                for run in range(len(self.dia_data_list[0].intensity_columns)):

                        sw_x = self.dia_data_list[0]
                        sw_y = self.dia_data_list[1]

                        if level == 'Transition':
                            df_x = sw_x.df
                            df_y = sw_y.df
                        elif level == 'Peptide':
                            df_x = sw_x.pep_df
                            df_y = sw_y.pep_df
                        else:
                            df_x = sw_x.prot_df
                            df_y = sw_y.prot_df


                        x_col = sw_x.intensity_columns[run]
                        y_col = sw_y.intensity_columns[run]

                        x_name = (sw_x.name_stem + ':   ' + x_col)
                        y_name = (sw_y.name_stem + ':   ' + y_col)

                        df = pd.DataFrame({x_name: np.log10(df_x[x_col]), y_name: np.log10(df_y[y_col])})

                        """
                        moving_stats_df = plotting.scatter_sliding_stats(df, x=x_name, y=y_name,
                                                                         bin_divisor=50, num_steps=100, min_counts=10,
                                                                         analyses='scatter+percentiles', # 'scatter+log_log_y=mx+b'
                                                                         verbose=self.verbose)
                        """

                        plotting.pandas_plot(df, x=x_name, y=y_name, kind='scatter',
                                             # ylim=None, xlim=None,
                                             title=('Intensity vs. Intensity - ' + level + ' - Run: \n' + x_col),
                                             xlabel=(sw_x.name_stem + ' log10(intensity)'),
                                             ylabel=(sw_y.name_stem + ' log10(intensity)'),
                                             grid_maj_x=False,
                                             grid_maj_y=False,
                                             unmatched_as_offaxis_jitter=True,
                                             #legend=True,
                                             #analysis='scatter+median',  #'scatter+log_log_y=mx+b'
                                             #moving_stats_df=moving_stats_df,
                                             height_inches=height_inches,
                                             width_inches=width_inches,
                                             show_plots=show_plots,
                                             plots_to_file=plots_to_file,
                                             tables_to_file=tables_to_file,
                                             # metadata_files=metadata_files,
                                             outpath=outpath,
                                             file_stem=str('Intensity Comp - ' + level + ' - Run - ' + str(run)),
                                             verbose=self.verbose
                                             )

            self.tracker.capture('run vs. run plotting')



        # -------------------------------------------
        # correlation of analytes plotting

        if aligned_analyte_correls is True:

            if self.verbose > 0:
                boiler_plates.section_header(title='Correlation of Analytes', level=2)

            if plots_to_file is True:
                outpath = os.path.join(self.outpath, 'analyte_correl')
                if not os.path.exists(outpath):
                    os.makedirs(outpath)
            else:
                outpath = None

            # Get correlation of each row for each each pair of DiaData class input sets
            for i in range(len(self.dia_data_list)):
                for j in range(len(self.dia_data_list)):

                    # only do comparision when they are different and only once, obviously
                    if i > j:

                        sample_columns = self.dia_data_list[i].intensity_columns

                        i_name = self.dia_data_list[i].name_stem
                        j_name = self.dia_data_list[j].name_stem

                        df_tuples = [(self.dia_data_list[i].df, self.dia_data_list[j].df, 'Transition'),
                                     (self.dia_data_list[i].pep_df, self.dia_data_list[j].pep_df, 'Peptide'),
                                     (self.dia_data_list[i].prot_df, self.dia_data_list[j].prot_df, 'Protein'),
                                     ]


                        for (df_i, df_j, level) in df_tuples:

                            if self.verbose > 0:
                                print()
                                print()
                                print('Data level:                 ', level)
                                print('control:                    ', i_name)
                                print('test:                       ', j_name)
                                print('shape of each:              ', df_i[sample_columns].shape)

                            # start a dataframe with the average log intensities from i
                            corr_anal_df = df_i[['ave_log_int']].copy(deep=True)

                            # compute pair-wise correlation of each row between the two result sets
                            # this asks if they are saying the same thing about biology examined, if a real experiment

                            corr_anal_df['r'] = df_i[sample_columns].corrwith(df_j[sample_columns],
                                                                              method='pearson',
                                                                              axis=1)

                            if self.verbose > 0:
                                print()
                                print('count correlation measures: ', corr_anal_df['r'].count())

                            # Count number of non-nan values common for each row (intersection of data completeness)
                            # (determining have a value in both matrices by subtraction here - could be any operation)
                            corr_anal_df['intersection_counts'] = \
                                (df_i[sample_columns] - df_j[sample_columns]).notnull().sum(axis=1)

                            # Mask out correl values for any cases with low common n
                            corr_anal_df['r'].mask((corr_anal_df['intersection_counts'] < 10), inplace=True)

                            if self.verbose > 0:
                                print('count w at least 10 pairs:  ', corr_anal_df['r'].count())
                                print()

                            # for dev troubleshooting
                            #corr_outfile = os.path.join(outpath, (level + '_corr_table.txt'))
                            #corr_anal_df.to_csv(corr_outfile, sep='\t', index=False)


                            # -----------------------------------------------
                            # plot correl vs. log average intensity of the analyte

                            # get moving stats to overlay on plot
                            moving_stats_df = plotting.scatter_sliding_stats(corr_anal_df, x='ave_log_int', y='r',
                                                                             analyses='scatter+percentiles',
                                                                             verbose=self.verbose)


                            plotting.pandas_plot(corr_anal_df, x='ave_log_int', y='r', kind='scatter',
                                                 title=(level + ' Corr v Ave Int\n' + i_name + ' v ' + j_name),
                                                 xlabel='log10(ave intensity)',
                                                 ylabel='Pearson r',
                                                 xlim=None,
                                                 ylim=(-1, 1),
                                                 grid_maj_x=True,
                                                 grid_maj_y=True,
                                                 moving_stats_df=moving_stats_df,
                                                 analysis='scatter+percentiles',
                                                 height_inches=height_inches,
                                                 width_inches=width_inches,
                                                 show_plots=show_plots,
                                                 plots_to_file=plots_to_file,
                                                 tables_to_file=tables_to_file,
                                                 # metadata_files=metadata_files,
                                                 outpath=outpath,
                                                 file_stem=(level + ' Corr v Ave Int - ' + i_name + ' v ' + j_name),
                                                 verbose=self.verbose
                                                 )


                            # compute histogram of the correlation values

                            plotting.histogram_continuous(corr_anal_df, 'r',
                                                          hist=True,
                                                          norm_hist=True,
                                                          norm_hist_as_prob=False,
                                                          bins=plotting.generate_bin_series(-1, 1, 1 / 100),
                                                          kde=True,
                                                          rug=False,
                                                          title=(level + ' Corr Distrib\n' + i_name + ' v ' + j_name),
                                                          ylabel='frequency',
                                                          xlabel='Pearson r',
                                                          xlim=None,
                                                          grid_maj_x=True,
                                                          grid_maj_y=False,
                                                          height_inches=height_inches,
                                                          width_inches=width_inches,
                                                          show_plots=show_plots,
                                                          plots_to_file=plots_to_file,
                                                          tables_to_file=tables_to_file,
                                                          # metadata_files=metadata_files,
                                                          outpath=outpath,
                                                          file_stem=(level + ' Corr Distrib - ' + i_name + ' v ' + j_name),
                                                          verbose=self.verbose)

            self.tracker.capture('correlation of analytes')

        return




    def compare_aligned_analyte_aggregates(self, **kwargs):

        show_plots = kwargs.get('show_plots', self.show_plots)
        plots_to_file = kwargs.get('plots_to_file', self.plots_to_file)
        tables_to_file = kwargs.get('tables_to_file', self.tables_to_file)
        height_inches = kwargs.get('height_inches', self.plot_height_inches)
        width_inches = kwargs.get('width_inches', self.plot_width_inches)

        if self.verbose > 0:
            boiler_plates.section_header(title='Comparison of Aligned Analyte Aggregate Properties', level=2)


        if plots_to_file is True:
            outpath = os.path.join(self.outpath, 'aligned_analytes')
            if not os.path.exists(outpath):
                os.makedirs(outpath)
        else:
            outpath = None


        # Do analysis for each pair of DiaData class input sets
        for i in range(len(self.dia_data_list)):
            for j in range(len(self.dia_data_list)):

                # only do comparision when they are different and only once, obviously
                if i > j:

                    i_name = self.dia_data_list[i].name_stem
                    j_name = self.dia_data_list[j].name_stem

                    df_tuples = [(self.dia_data_list[i].df, self.dia_data_list[j].df, 'Transition'),
                                 (self.dia_data_list[i].pep_df, self.dia_data_list[j].pep_df, 'Peptide'),
                                 (self.dia_data_list[i].prot_df, self.dia_data_list[j].prot_df, 'Protein'),
                                 ]


                    # iterate over data levels
                    for (df_i, df_j, level) in df_tuples:

                        if self.verbose > 0:
                            print()
                            print()
                            print('Data level:                 ', level)
                            print('control (x-axis):           ', i_name)
                            print('test (y-axis):              ', j_name)

                        # iterate over properties to compare
                        for metric in ['cv_int', 'f_complete', 'ave_qval', 'ldr', 'ave_acc', 'ave_rms', 'ave_delta_rt']:

                            # skip the plots if the metric is suppressed for this level
                            # (ex. q-values are in data at trans level for filtering sake, but not worth plotting)
                            if level not in self.plot_param_dict[metric].get('suppress_levels', set()):

                                if self.verbose > 0:
                                    print('     metric:', metric)

                                # metrics where need to use max_int instead of ave_int (one never from replicates)
                                max_int_metrics = ['ldr', 'ave_acc', 'ave_rms']

                                if metric in max_int_metrics:
                                    log_int_metric = 'max_log_int'
                                else:
                                    log_int_metric = 'ave_log_int'

                                # only plot if metric present in both sets
                                if metric in df_i and metric in df_j:

                                    # start a dataframe with the metric from i and intensity
                                    metric_df = df_i[[metric, log_int_metric]].copy(deep=True)

                                    # rename the metric col as the set name
                                    metric_df.rename(columns={metric: i_name}, inplace=True)

                                    # add the metric from the second set
                                    metric_df[j_name] = df_j[metric]

                                    # compute the ratio of the metric between the two sets as (y / x)
                                    metric_df['y/x ratio'] = df_j[metric] / df_i[metric]
                                    metric_df['log y/x ratio'] = np.log10(metric_df['y/x ratio'])


                                    # -----------------------------------------------
                                    # plot metric in set a vs. set b

                                    title = level + ': Comparison of Metric: ' + metric + '\n' + i_name + ' v ' + j_name
                                    file_stem = (level + '_' + metric + '_comparison__' + i_name + '_v_' + j_name)

                                    plotting.pandas_plot(metric_df,
                                                         x=i_name, y=j_name, kind='scatter',
                                                         title=title,
                                                         xlabel=(metric + ' for ' + i_name),
                                                         ylabel=(metric + ' for ' + j_name),
                                                         xlim=self.plot_param_dict[metric].get('axis_range', None),
                                                         ylim=self.plot_param_dict[metric].get('axis_range', None),
                                                         grid_maj_x=True,
                                                         grid_maj_y=True,
                                                         equal_aspect=True,
                                                         unmatched_as_offaxis_jitter=True,
                                                         height_inches=height_inches,
                                                         width_inches=width_inches,
                                                         show_plots=show_plots,
                                                         plots_to_file=plots_to_file,
                                                         tables_to_file=tables_to_file,
                                                         outpath=outpath,
                                                         file_stem=file_stem,
                                                         verbose=self.verbose
                                                         )


                                    # -----------------------------------------------
                                    # plot metric ratio in set a vs. set b vs. intensity

                                    # get moving stats to overlay on plot
                                    moving_stats_df = plotting.scatter_sliding_stats(metric_df, x=log_int_metric,
                                                                                     y='log y/x ratio',
                                                                                     analyses='scatter+percentiles',
                                                                                     verbose=self.verbose)

                                    title = level + ': ' + metric + ' Ratio v Intensity\n' + j_name + ' : ' + i_name
                                    file_stem = (level + '_' + metric + '_ratio_v_int__' + j_name + '_v_' + i_name)

                                    plotting.pandas_plot(metric_df,
                                                         x=log_int_metric, y='log y/x ratio', kind='scatter',
                                                         title=title,
                                                         xlabel='average log10(intensity)',
                                                         ylabel=('log10 ' + j_name + ' / ' + i_name),
                                                         xlim=None,
                                                         ylim=self.plot_param_dict[metric].get('log_ratio_axis_range', None),
                                                         grid_maj_x=True,
                                                         grid_maj_y=True,
                                                         moving_stats_df=moving_stats_df,
                                                         analysis='scatter+percentiles',
                                                         height_inches=height_inches,
                                                         width_inches=width_inches,
                                                         show_plots=show_plots,
                                                         plots_to_file=plots_to_file,
                                                         tables_to_file=tables_to_file,
                                                         # metadata_files=metadata_files,
                                                         outpath=outpath,
                                                         file_stem=file_stem,
                                                         verbose=self.verbose
                                                         )


                                    # -----------------------------------------------
                                    # compute histogram of the ratios

                                    title = level + ': Ratio of Metric: ' + metric + '\n' + j_name + ' : ' + i_name
                                    file_stem = (level + '_' + metric + '_ratio__' + i_name + '_v_' + j_name)

                                    range_tuple = self.plot_param_dict[metric].get('log_ratio_axis_range', None)
                                    if range_tuple is not None:
                                        hist_min, hist_max = range_tuple
                                    else:
                                        hist_min = -3
                                        hist_max = 3


                                    plotting.histogram_continuous(metric_df, 'log y/x ratio',
                                                                  hist=True,
                                                                  norm_hist=True,
                                                                  norm_hist_as_prob=False,
                                                                  bins=plotting.generate_bin_series(hist_min, hist_max, 1 / 100),
                                                                  kde=True,
                                                                  rug=False,
                                                                  title=title,
                                                                  ylabel='frequency',
                                                                  xlabel=('log y/x ratio ' + metric),
                                                                  xlim=None,
                                                                  grid_maj_x=True,
                                                                  grid_maj_y=False,
                                                                  height_inches=height_inches,
                                                                  width_inches=width_inches,
                                                                  show_plots=show_plots,
                                                                  plots_to_file=plots_to_file,
                                                                  tables_to_file=tables_to_file,
                                                                  # metadata_files=metadata_files,
                                                                  outpath=outpath,
                                                                  file_stem=file_stem,
                                                                  verbose=self.verbose)


        self.tracker.capture('comparisons of analyte aggregate properties')

        return




    def compare_unaligned_analytes_aggregates(self, **kwargs):

        show_plots = kwargs.get('show_plots', self.show_plots)
        plots_to_file = kwargs.get('plots_to_file', self.plots_to_file)
        tables_to_file = kwargs.get('tables_to_file', self.tables_to_file)
        height_inches = kwargs.get('height_inches', self.plot_height_inches)
        width_inches = kwargs.get('width_inches', self.plot_width_inches)

        if self.verbose > 0:
            boiler_plates.section_header(title='Overlays of Unaligned Analyte Aggregate Properties', level=2)


        if plots_to_file is True:
            outpath = os.path.join(self.outpath, 'unaligned_analytes')
            if not os.path.exists(outpath):
                os.makedirs(outpath)
        else:
            outpath = None


        # iterate over properties to compare
        for metric in ['cv_int', 'f_complete', 'ave_qval', 'ldr', 'ave_acc', 'ave_rms', 'ave_delta_rt']:

            # check if cum hist should be descending
            descending = self.plot_param_dict[metric].get('descending', False)

            # skip the plots if the metric is suppressed for this level
            # (ex. q-values are in data at trans level for filtering sake, but not worth plotting)
            suppress_levels = self.plot_param_dict[metric].get('suppress_levels', set())

            # check if any max placed on cumulative range (to ignore outliers and focus the analysis)
            cum_max = self.plot_param_dict[metric].get('cum_max', None)

            # check if completeness weighted analysis should be skipped for this metric
            # ex. in spiked species analyte, the completeness isn't expected to be good so makes no sense to weight
            #     by completeness for metrics derived from these like LDR, acc, rms
            suppress_comp_weighting = self.plot_param_dict[metric].get('suppress_comp_weighting', False)


            # get each level

            # make a df with one col per input DiaData set
            trans_metric_df = pd.DataFrame()
            pep_metric_df = pd.DataFrame()
            prot_metric_df = pd.DataFrame()

            trans_weights_df = pd.DataFrame()
            pep_weights_df = pd.DataFrame()
            prot_weights_df = pd.DataFrame()

            # iterate over all DiaData input sets to build parallel metric dfs at all levels
            for i in range(len(self.dia_data_list)):

                name_stem = self.dia_data_list[i].name_stem

                if metric in self.dia_data_list[i].df:
                    trans_metric_df[name_stem] = self.dia_data_list[i].df[metric]
                    trans_weights_df[name_stem] = self.dia_data_list[i].df['counts']

                if metric in self.dia_data_list[i].pep_df:
                    pep_metric_df[name_stem] = self.dia_data_list[i].pep_df[metric]
                    pep_weights_df[name_stem] = self.dia_data_list[i].pep_df['counts']

                if metric in self.dia_data_list[i].prot_df:
                    prot_metric_df[name_stem] = self.dia_data_list[i].prot_df[metric]
                    prot_weights_df[name_stem] = self.dia_data_list[i].prot_df['counts']

            # compute cumulative distribution

            # initialize
            trans_cum_df = None
            trans_cum_df_w = None
            pep_cum_df = None
            pep_cum_df_w = None
            prot_cum_df = None
            prot_cum_df_w = None

            if len(trans_metric_df.columns) > 0 and 'Transition' not in suppress_levels:

                # determine operating x_max by considering max in data also
                if cum_max is not None:
                    max_eff = min(cum_max, trans_metric_df.max().max())
                else:
                    max_eff = None

                trans_cum_df = stat_methods.cum_freqs_df(trans_metric_df, num_thresh=200,
                                                         descending=descending, x_max=max_eff)

                if suppress_comp_weighting is False:
                    trans_cum_df_w = stat_methods.cum_freqs_df(trans_metric_df, weights_df=trans_weights_df,
                                                               num_thresh=200, descending=descending, x_max=max_eff)


            if len(pep_metric_df.columns) > 0 and 'Peptide' not in suppress_levels:

                # determine operating x_max by considering max in data also
                if cum_max is not None:
                    max_eff = min(cum_max, pep_metric_df.max().max())
                else:
                    max_eff = None

                pep_cum_df = stat_methods.cum_freqs_df(pep_metric_df, num_thresh=200,
                                                       descending=descending, x_max=max_eff)

                if suppress_comp_weighting is False:
                    pep_cum_df_w = stat_methods.cum_freqs_df(pep_metric_df, weights_df=pep_weights_df,
                                                             num_thresh=200, descending=descending, x_max=max_eff)


            if len(prot_metric_df.columns) > 0 and 'Protein' not in suppress_levels:

                # determine operating x_max by considering max in data also
                if cum_max is not None:
                    max_eff = min(cum_max, prot_metric_df.max().max())
                else:
                    max_eff = None

                prot_cum_df = stat_methods.cum_freqs_df(prot_metric_df, num_thresh=200,
                                                        descending=descending, x_max=max_eff)

                if suppress_comp_weighting is False:
                    prot_cum_df_w = stat_methods.cum_freqs_df(prot_metric_df, weights_df=prot_weights_df,
                                                              num_thresh=200, descending=descending, x_max=max_eff)


            # get just the upper x limit (lower includes allowance for jittered unmatched offaxis)
            x_lims = self.plot_param_dict[metric].get('axis_range', None)
            if x_lims is not None:
                xlim = (0, x_lims[1])
            else:
                xlim = (0, None)


            for cum_df, level, weighted in \
                    [(trans_cum_df, 'Transition', False), (pep_cum_df, 'Peptide', False),
                     (prot_cum_df, 'Protein', False),

                     (trans_cum_df_w, 'Transition', True), (pep_cum_df_w, 'Peptide', True),
                     (prot_cum_df_w, 'Protein', True)]:

                if cum_df is not None:

                    # -----------------------------------------------
                    # extract yields at critical points

                    data_set_list = cum_df.columns.tolist()
                    data_set_list.remove('threshold')

                    critical_values = self.plot_param_dict[metric]['critical_values']

                    if weighted is False:
                        weighting = ''
                    else:
                        weighting = 'completeness'

                    for set_name in data_set_list:
                        yields = stat_methods.find_x_yields_at_critical_ys(critical_values, cum_df[set_name],
                                                                           cum_df['threshold'])
                        for i in range(len(critical_values)):
                            self.capture_result_value('unaligned_cum_yields',
                                                      set_name=set_name, level=level, metric=metric, weighting=weighting,
                                                      thresh=critical_values[i], yield_value=yields[i])


                    # -----------------------------------------------
                    # plot of line from each input set overlaid

                    if plots_to_file is True or show_plots is True:

                        if weighted is False:
                            title = level + ': Cumulative Distribs: ' + metric
                            file_stem = (level + '_' + metric + '_cum_dists')
                            ylabel = 'cumulative freq'
                        else:
                            title = level + ': Completeness Weighted Cum. Distribs: ' + metric
                            file_stem = (level + '_' + metric + '_cum_dists_weighted')
                            ylabel = 'completeness weighted cumulative freq'

                        plotting.pandas_plot(cum_df, x='threshold',
                                             kind='line',
                                             title=title,
                                             xlabel=metric,
                                             ylabel=ylabel,
                                             xlim=xlim,
                                             ylim=(0, None),
                                             legend=True,
                                             grid_maj_x=True,
                                             grid_maj_y=True,
                                             equal_aspect=False,
                                             unmatched_as_offaxis_jitter=False,
                                             height_inches=height_inches,
                                             width_inches=width_inches,
                                             show_plots=show_plots,
                                             plots_to_file=plots_to_file,
                                             tables_to_file=tables_to_file,
                                             outpath=outpath,
                                             file_stem=file_stem,
                                             verbose=self.verbose
                                             )


        self.tracker.capture('overlays of unaligned analyte aggregate properties')

        return




    def capture_result_value(self, result_type, **kwargs):
        """
        Capture a results data point where the only required arg is the type of result, which determines which sub-dict
        the result goes to. After that, kwargs are used so different result structures can be used in each sub-dict.
        Could have made separate results objects but expect common mechanics here.

        """

        set_name = kwargs.get('set_name', None)

        level = kwargs.get('level', None)
        metric = kwargs.get('metric', None)
        weighting = kwargs.get('weighting', None)
        thresh = kwargs.get('thresh', None)

        yield_value = kwargs.get('yield_value', None)

        if metric is not None:
            category = self.plot_param_dict[metric].get('category', 'missing')

        if result_type == 'unaligned_cum_yields':

            id_value = level + metric + weighting + str(thresh)

            # initialize columns if first time
            if self.results['unaligned_cum_yields'].shape[0] == 0:

                self.results['unaligned_cum_yields'] = \
                    pd.DataFrame(columns=['category', 'metric', 'level', 'weighting', 'thresh', set_name])

            # if first row
            if self.results['unaligned_cum_yields'].shape[0] == 0 or \
                id_value not in self.results['unaligned_cum_yields'].index.values:

                row_series = pd.Series({'category': category, 'metric': metric, 'level': level, 'weighting': weighting,
                                        'thresh': thresh, set_name: yield_value},
                                       name=id_value)

                self.results['unaligned_cum_yields'] = \
                    self.results['unaligned_cum_yields'].append(row_series, ignore_index=False)

            elif set_name not in self.results['unaligned_cum_yields'].columns.values:
                self.results['unaligned_cum_yields'][set_name] = pd.Series()
                self.results['unaligned_cum_yields'].loc[id_value, set_name] = yield_value

            else:
                self.results['unaligned_cum_yields'].loc[id_value, set_name] = yield_value

        else:
            print('FATAL ERROR! Unknown results type:', result_type, 'given to capture_result_value() method.')
            quit()

        return




    def unaligned_cum_yields_table(self, **kwargs):

        include_completeness_weighting = kwargs.get('include_completeness_weighting', )

        df = self.results['unaligned_cum_yields']

        sort_columns = ['level', 'category', 'metric', 'weighting', 'thresh']
        set_name_list = [col for col in df.columns if col not in sort_columns]

        # sort the table
        df.sort_values(sort_columns, ascending=True, inplace=True)


        # -------------------------------
        # find and remove constant prefix, if present (assuming underscore separation of elements in string names)

        done = False
        prefix = None
        num_elements = 0
        while done is False:
            num_elements += 1
            if sum([1 for col in set_name_list if
                    col.split('_')[:num_elements] == set_name_list[0].split('_')[:num_elements]]) == len(set_name_list):
                element_list = set_name_list[0].split('_')[:num_elements]
                prefix = '_'
                prefix = prefix.join(element_list)
                prefix += '_'
            else:
                done = True

        # if found a prefix, remove it from all headers - both the list and actual df
        if len(prefix) > 0:
            renaming_dict = {x: x.replace(prefix, '') for x in set_name_list}
            set_name_list = [x.replace(prefix, '') for x in set_name_list]
            df.rename(columns=renaming_dict, inplace=True)


        # -------------------------------
        # compute relative performance columns

        # compute ave col
        df['ave_yield'] = df[set_name_list].mean(axis=1, skipna=True)

        diff_cols = [(set_name + '_%\u0394') for set_name in set_name_list]

        for set_name in set_name_list:
            df[(set_name + '_%\u0394')] = round(100 * (df[set_name] - df['ave_yield']) / df['ave_yield'], 1)


        # -------------------------------
        # clean up master table

        df.drop(['ave_yield'], axis=1, inplace=True)

        if include_completeness_weighting is False:
            df.drop(df[df['weighting'] == 'completeness'].index, inplace=True)
            df.drop(['weighting'], axis=1, inplace=True)

        df[set_name_list] = df[set_name_list].round(0)


        # -------------------------------
        # build subset dataframes by level

        trans_df = df[df['level'] == 'Transition'].copy(deep=True)
        pep_df = df[df['level'] == 'Peptide'].copy(deep=True)
        prot_df = df[df['level'] == 'Protein'].copy(deep=True)

        for df in [trans_df, pep_df, prot_df]:
            df.drop(['level'], axis=1, inplace=True)


        # -------------------------------
        # distill dfs by computing average relative yields within each category

        trans_sum_df = trans_df.groupby(['category']).mean()
        pep_sum_df = pep_df.groupby(['category']).mean()
        prot_sum_df = prot_df.groupby(['category']).mean()

        for df in [trans_sum_df, pep_sum_df, prot_sum_df]:
            df.reset_index(inplace=True)
            df.drop(set_name_list + ['thresh'], axis=1, inplace=True)


        # -------------------------------
        # output

        if self.verbose > 0:

            for df, sum_df, level in [(trans_df, trans_sum_df, 'Transition'),
                                      (pep_df, pep_sum_df, 'Peptide'),
                                      (prot_df, prot_sum_df, 'Protein')]:

                if prefix is None:
                    title = level + ' Level Yields from Unaligned Data'
                else:
                    title = level + ' Level Yields from Unaligned Data - ' + prefix.strip('_')

                boiler_plates.section_header(title = title, level=2)
                print(tabulate(df, headers='keys', tablefmt='simple', showindex=False))
                print()
                print()
                print(level + ' Level Relative Yield Averages')
                print()
                print(tabulate(sum_df, headers='keys', tablefmt='simple', showindex=False ))
                print()
                print()
                print()

        return




    def summarize(self):

        if self.verbose > 0:
            boiler_plates.section_header(title='Comparison Summary', level=1)

        # todo - create and call output for aligned analyte results

        self.unaligned_cum_yields_table(include_completeness_weighting=False)


        return




    def profile(self):

        self.tracker.report(title='Timing Profile for DiaComparison Class', verbose=self.verbose)

        return


    # ---------------- END CLASS ----------------------------------------------------




