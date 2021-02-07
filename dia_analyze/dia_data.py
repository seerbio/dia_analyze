__author__ = 'Seymour Data Science'


import os.path
from pathlib import Path

import numpy as np
import pandas as pd

import re

from tabulate import tabulate

from crozes_base.helpers import helper_functions, pandas_helpers
from crozes_base.io import boiler_plates
from crozes_base.stats import plotting, stat_methods

from seer_swath_analyze import modifications, quant


"""
GUIDANCE ON WHERE TO UPDATE TO SUPPORT ADDITIONAL FORMATS etc

- All relevant areas of code for update have been tagged with 'EXTENSION_ZONE'. Search on this to find areas mentioned
 in the following comments.

- The mapping dictionaries in the convert_format() method are critical. Note that there are separate dicts for results
vs libraries as the column headers can be different within the same software

- Presently this assumes a single dialect from each software. This may not be true if column renaming is used - ex in 
Spectronaut.

- Control of which removal of rows using 'exclude form assay' type indicator columns is controlled in the 
remove_excluded() method, which has a list that enumerates all such columns.

- This class must be able to determine library vs. results, data level (trans/pep/prot) and the dialect of the data
automatically. This is done in the method categorize_input(), which has lots of hard-coded logic to this end.

- Any new columns added to the std format must be included in the unstack() uses of pivot in order to keep these columns

"""




class SwathData(object):
    """
    A SWATH data object - can be a library (lacking quant data) or a results set (containing quant data on samples).

    """


    def __init__(self, swath_input, **kwargs):
        """
        Create a new SwathData class instance.

        :param input: A file path or pandas dataframe.
        :param: kwarg: in_delim: delimiter in input file, if filepath passed rather than dataframe.

        """

        self.name_stem = kwargs.get('name_stem', None)
        self.outpath = kwargs.get('outpath', None)
        self.verbose = kwargs.get('verbose', 0)
        in_delim = kwargs.get('in_delim', None)
        self.output_format = kwargs.get('output_format', None)

        self.include_species = kwargs.get('include_species', [])
        self.exclude_sample_strings = kwargs.get('exclude_sample_strings', [])

        self.ignore_provided_rollup = kwargs.get('ignore_provided_rollup', False)

        self.ignore_provided_norm = kwargs.get('ignore_provided_norm', False)

        self.expt_type = kwargs.get('expt_type', None)  # options: 'replicates', 'dilution_series_species_spike'
        self.species_column = kwargs.get('species_column', None)
        self.const_species = kwargs.get('const_species', None)
        self.var_species = kwargs.get('var_species', None)
        self.dil_concs = kwargs.get('dil_concs', {})
        self.max_percent_rms_ldr = kwargs.get('max_percent_rms_ldr', 20.0)
        self.acc_rms_fixed_dil = kwargs.get('acc_rms_fixed_dil', False)

        self.plot_height_inches = kwargs.get('plot_height_inches', 4)
        self.plot_width_inches = kwargs.get('plot_width_inches', 6)

        # instantiate the timing and RAM use tracker
        self.tracker = helper_functions.Time_and_RAM_Tracker()


        # instantiate analytetracker
        self.analyte_tracker = []       # list to collect lists of measurements

        # error trackers
        self.ldr_errors = 0
        self.ldr_successes = 0
        self.ldr_skipped = 0

        # Get mod catalog
        self.mod_dict = modifications.get_mod_catalog()

        # initializations
        self.prec_df = None
        self.pep_df = None
        self.prot_df = None
        self.mod_counts_dict = {}
        self.semantic_dict = {}     # given a std_name key, get the col name for where that info is currently
        self.sample_dilutions = None

        self.dialect = 'unknown'
        self.level = 'unknown'
        self.lib_or_results = 'unknown'
        self.stacked = False
        self.modified_seq_col = 'unknown'
        self.mod_syntax = 'unknown'
        self.provided_norm_present = 'unknown'
        self.normalized_status = False
        self.mz_indexing = False

        self.intensity_columns = []
        self.score_columns = []
        self.q_value_columns = []
        self.rt_columns = []
        self.pred_rt_columns = []
        self.delta_rt_columns = []
        self.sample_name_stems = []

        self.intensity_suffix = ' Intensity'
        self.score_suffix = ' Score'
        self.q_value_suffix = ' q-value'
        self.obs_rt_suffix = ' Obs_RT'
        self.pred_rt_suffix = ' Pred_RT'
        self.delta_rt_suffix = ' Delta_RT'

        # build a collection of these as tuples: (data_type, column_list, suffix) in
        self.data_arrays_list = []


        self.sample_norm_factors = {}   # dict just to be safe to avoid sorting out of sync issues
        self.char_dict = {}             # characterization stats dict


        # ---------------------------------
        # Get dataframe(s)

        # if input is a string, assume it's a file path
        if type(swath_input) == str:

            # if ends in .fcresult, is a Sciex OneOmics 'fold change' results zip file
            if swath_input.endswith('.fcresult') is True:
                self.read_one_omics_fcresult(swath_input)

            # otherwise assume is a single text file and read to a pandas dataframe
            else:
                self.read_input(swath_input, delimiter=in_delim)

            if self.name_stem is None:
                self.name_stem = Path(swath_input).stem

        elif type(swath_input) == pd.DataFrame:
            # todo - this usage path untested and unclear if needed.
            self.df = swath_input
            self.name_stem = 'swath_data_from_pkl'
        else:
            print('FATAL ERROR! Input to SwathData class must be a pandas dataframe or a file path')
            quit()

        if self.verbose > 0:
            print('Input df shape:                ', self.df.shape)


        # ---------------------------------
        # Basic properties

        # determine and initialize properties
        self.categorize_input()

        self.input_format = self.dialect
        self.converted = False

        # get the column that contains the sequence with modifications - ie complete peptide indication
        self.get_mod_seq_col()


        # ---------------------------------
        # Convert format

        self.describe(title='Describe: Initial SWATH data input before conversion')
        self.convert_format()
        self.describe(title='Describe: Immediately after conversion')


        # ---------------------------------
        # Non-optional post-conversion steps

        # remove rows that do not have required species if set
        self.filter_by_species()

        # determine experimental role of samples
        self.determine_experimental_structure()




    def describe(self, **kwargs):
        """
        This is meant to be a high level description, not detailed. Use self.characterize() for detailed view.

        """

        title = kwargs.get('title', 'Describe SWATH Data: ' + str(self.name_stem))

        if self.verbose > 0:
            boiler_plates.section_header(title=title, level=1)
            print('Data type (lib or result):    ', self.lib_or_results)
            print('Dialect:                      ', self.dialect)
            print('Mod syntax:                   ', self.mod_syntax)
            print('Stacked (long) format:        ', self.stacked)
            print('Sample normalization status:  ', self.normalized_status)
            print()
            print('Sample columns:               ', len(self.intensity_columns))
            print()
            print('Dataframes present:           ')
            print('                   transition:', self.df.shape)
            if self.prec_df is not None:
                print('                   precursor: ', self.prec_df.shape)
            if self.pep_df is not None:
                print('                   peptide:   ', self.pep_df.shape)
            if self.prot_df is not None:
                print('                   protein:   ', self.prot_df.shape)
            print()
            if self.verbose > 1 and len(self.intensity_columns) > 0:
                print('Sample columns')
                print()
                for x in self.intensity_columns:
                    print('  ', x)

            if self.verbose > 2:
                print('Sample of first 10 rows:')
                print()
                print(tabulate(self.df.head(10), headers='keys', tablefmt='simple'))
                print()
                print()

        # todo - add table with df type vs data array type and num columns, num non-empty value present - two tables?

        return




    def read_input(self, input_file, **kwargs):

        delimiter = kwargs.get('delimiter', None)

        if self.verbose > 0:
            boiler_plates.section_header(title='Read input from file', level=2)


        # -----------------------------------------------------------------------------------
        # Infer delimiter from filename (if None, file read is way slower when falling back to Python inferred delim)

        # Use explicit delimiter if provided. Otherwise infer from filename here
        if delimiter is None:
            delimiter, format_name = helper_functions.infer_delimiter_or_format(input_file)
        else:
            format_name = None


        # -----------------------------------------------------------------------------------
        # Read

        if delimiter is not None or ('.xls' not in input_file and '.xlsx' not in input_file):
            self.df = pd.read_csv(input_file, sep=delimiter)  # sep=None invokes automatic determination
        elif '.xls' in input_file or '.xlsx' in input_file:
            self.df = pd.read_excel(input_file)
        else:
            self.df = None
            print('FATAL ERROR!! File type could not be determined')
            quit()

        # this is just to grab a smaller sample of an input file too large to trim manually in excel
        #output_file = input_file.replace('.txt', '_SMALL.txt')
        #self.df.head(10000).to_csv(output_file, sep='\t')    # todo - remove!!!
        #quit()

        self.tracker.capture('read input')

        return




    def read_one_omics_fcresult(self, fcresult_file):
        """
        This method handles the output from the Sciex OneOmics cloud environment, captured as an .fcresult file.
        This file is a zip of a large number of text files. A small subset of these files are read, processeed,
        and structured into three data levels, allowing them to be processed the same way as data from other
        software.

        """

        def process_level(level):
            """
            Worker function to process each level of data

            """

            qval_data = None
            obs_rt_data = None
            pred_rt_data = None
            prefix = None

            if level == 'transition':
                area_data = 'TransitionQuantData-[Area].txt'
            elif level == 'peptide':
                area_data = 'PeptideQuantData-[Area].txt'
                qval_data = 'PeptideQuantData-[FDR].txt'
                obs_rt_data = 'PeptideQuantData-[ObservedRT].txt'
                pred_rt_data = 'PeptideQuantData-[PredictedRT].txt'
            else:
                area_data = 'ProteinQuantData-[Area].txt'


            combined_df = pd.DataFrame()

            for (data_type, file_name, suffix) in [('intensity', area_data, self.intensity_suffix),
                                                   ('qval', qval_data, self.q_value_suffix),
                                                   ('obs_rt', obs_rt_data, self.obs_rt_suffix),
                                                   ('pred_rt', pred_rt_data, self.pred_rt_suffix)]:

                if file_name is not None:

                    full_file_name = os.path.join(results_folder, file_name)

                    if os.path.exists(full_file_name):

                        # read data
                        df = pd.read_csv(full_file_name, sep='\t')

                        # get list of actual data cols
                        sample_cols = [x for x in df.columns if '.wiff' in x]


                        # if intensity array, check for norm'd data and select which to keep
                        if data_type == 'intensity':

                            raw_cols = [x for x in sample_cols if x.startswith('raw_')]
                            norm_cols = [x for x in sample_cols if x.startswith('normalized_')]

                            if len(norm_cols) > 0:
                                self.provided_norm_present = True
                            else:
                                self.provided_norm_present = False

                            # use the normalized columns if both present and instructed to do so
                            if self.provided_norm_present is True and self.ignore_provided_norm is False:
                                self.normalized_status = True
                                prefix = 'normalized_'
                                if len(raw_cols) > 0:
                                    df.drop(raw_cols, axis=1, inplace=True)
                                sample_cols = norm_cols
                            else:
                                self.normalized_status = False
                                prefix = 'raw_'
                                if len(norm_cols) > 0:
                                    df.drop(norm_cols, axis=1, inplace=True)
                                sample_cols = raw_cols


                        # convert all the quant cols from int64 to float64, setting zeros to nan except for qvals
                        if data_type != 'qval':
                            df[sample_cols] = \
                                df[sample_cols].astype('float64').replace(0, np.nan)
                        # will really have zeros in q-val and need to keep them!
                        # todo - check how blanks differ from zeros in qvals in the input files!
                        else:
                            df[sample_cols] = df[sample_cols].astype('float64')

                        # trim crap off col names and add suffix
                        renaming_dict = {x: x.replace('.qresult.export', suffix) for x in sample_cols}
                        sample_cols = [x.replace('.qresult.export', suffix) for x in sample_cols]

                        # remove raw or norm prefix if necessary
                        if data_type == 'intensity' and prefix is not None:
                            renaming_dict = {x: y.replace(prefix, '', 1) for (x, y) in renaming_dict.items()}
                            sample_cols = [x.replace(prefix, '', 1) for x in sample_cols]

                        df.rename(columns=renaming_dict, inplace=True)


                        # merge the different array types to one df for peptide data

                        # only have intensity data types for protein and transition levels w OneOmics data
                        if level != 'peptide':
                            combined_df = df

                        # for peptide level, need to potentially merge
                        else:
                            # if intensity, the just started so set combined equal to new df
                            if data_type == 'intensity':
                                combined_df = df

                            # if one of other levels, need to merge into combined df
                            else:
                                combined_df = pd.concat([combined_df, df[sample_cols]], axis=1, sort=False)

            return combined_df


        # ---------------------------------
        # Main data read


        # Verify have outpath
        if self.outpath is None:
            print('FATAL ERROR: An outpath argument is required when reading OneOmics .fcresult zips')
            print('             since a temp folder must be created to unzip the results.')
            print('             Pass an outpath argument when instantiating SwathData class.')
            results_folder = None
            quit()
        else:
            results_folder = os.path.join(self.outpath, 'UNZIPS', self.name_stem)


        # Unzip the .fcresult file
        helper_functions.unzip(fcresult_file, ziptopath=results_folder)

        # Extract the text file data to dataframes
        self.df = process_level('transition')
        self.pep_df = process_level('peptide')
        self.prot_df = process_level('protein')


        # ---------------------------------
        # Fill in transition level data not provided from peptide level (need to filter a trans level for stats)

        needed_cols = [x for x in self.pep_df.columns if x not in self.df.columns]


        # create a temporary peptide ID in both dfs
        self.df['prec_ID'] = \
            self.df['Protein Id'].astype(str) + '_' + self.df['Peptide HyperTag'].astype(str) + '_' \
            + self.df['Peptide Charge'].astype(str)

        self.pep_df['prec_ID'] = \
            self.pep_df['Protein Id'].astype(str) + '_' + self.pep_df['Peptide HyperTag'].astype(str) + '_' \
            + self.pep_df['Peptide Charge'].astype(str)


        # set these temp indices as the the df index in both
        self.df.set_index('prec_ID', drop=True, inplace=True)
        self.pep_df.set_index('prec_ID', drop=True, inplace=True)


        # fill missing data
        self.df[needed_cols] = self.pep_df[needed_cols]

        # delete the temp ID cols be resetting the indices on the dfs
        self.df.reset_index(drop=True, inplace=True)
        self.pep_df.reset_index(drop=True, inplace=True)


        # ---------------------------------
        # Final clean up

        # Set categorizations
        self.dialect = 'oneomics'
        self.level = 'transition'   # means the level of the lowest level, which is self.df
        self.lib_or_results = 'results'
        self.stacked = False


        self.tracker.capture('unzip fcresult file')

        return




    def determine_experimental_structure(self):


        if self.lib_or_results == 'library':
            return

        # if have dilution concentration metadata mapping info, map all samples
        if len(self.dil_concs) > 0:

            if self.verbose > 0:
                boiler_plates.section_header(title='Assign Dilution Concentrations to Samples', level=1)

            self.sample_dilutions = []  # this is None by default, so need to initialize dict

            for sample in self.intensity_columns:
                match_list = [k for k in self.dil_concs if k in sample]

                if len(match_list) == 1:

                    conc = self.dil_concs[match_list[0]]

                    # using a list structure because is in order of self.intensity_columns
                    self.sample_dilutions.append(conc)

                    if self.verbose > 0:
                        print('    Conc: ', conc, '   Sample:   ', sample)


                else:
                    print('FATAL ERROR! Dilution concentration mapping strings not specific enough!')
                    print('Sample:', sample, 'matches to multiple or none in match list:', match_list)
                    quit()


        return




    def categorize_input(self):
        """
        Note: Initially the design allowed input at any level to start as self.df, the main table. However, the design
        has evolved to assume self.df is always transition level, self.pep_df is always peptide, and self.prot_df is
        always the protein level

        """

        # todo - refactor to adapt for change in comment in doc string above

        cols = self.df.columns

        # EXTENSION_ZONE
        if self.dialect == 'unknown':

            if 'PG.ProteinAccessions' in cols or 'ExcludeFromAssay' in cols or 'R.Condition' in cols:
                self.dialect = 'spectronaut'

                # Get level of Spectronaut data
                if 'F.FrgMz' in cols or 'FragmentMz' in cols:
                    self.level = 'transition'
                    if 'F.PeakArea' in cols:
                        self.lib_or_results = 'results'
                        self.stacked = True
                    else:
                        self.lib_or_results = 'library'

                elif 'FG.PrecMz' in cols:
                    self.level = 'peptide'
                    self.lib_or_results = 'results'

                else:
                    self.level = 'protein'
                    self.lib_or_results = 'results'

            elif 'Protein Locator' in cols:   # this assumes using the custom reports I made to get feature tables
                self.dialect = 'skyline'

                # Get level of Skyline data
                if 'Product Mz' in cols:
                    self.level = 'transition'

                    sample_columns = [x for x in cols if 'Area' in x]

                    if len(sample_columns) > 0:
                        self.lib_or_results = 'results'
                    else:
                        self.lib_or_results = 'library'
                elif 'Precursor Mz' in cols:
                    self.level = 'peptide'
                    self.lib_or_results = 'results'
                else:
                    self.level = 'protein'
                    self.lib_or_results = 'results'

            elif 'Protein' in cols or 'uniprot_id' in cols:
                self.dialect = 'sciex'

                # Get level of PeakView SWATH export data
                if 'Fragment MZ' in cols or 'Q3' in cols:
                    self.level = 'transition'
                    # count columns w 'area' in them
                    col_w_area = sum([1 for x in self.df.columns if 'Area' in x])


                    if col_w_area > 0:
                        self.lib_or_results = 'results'
                    else:
                        self.lib_or_results = 'library'
                elif 'Precursor MZ' in cols:
                    self.level = 'peptide'
                    self.lib_or_results = 'results'
                else:
                    self.level = 'protein'
                    self.lib_or_results = 'results'


        if self.verbose > 2:
            print()
            print('Categorization of input:')
            print('   dialect:        ', self.dialect)
            print('   level:          ', self.level)
            print('   lib or results: ', self.lib_or_results)
            print('   stacked:        ', self.stacked)

        # catch currently non-functioning case
        if self.ignore_provided_norm is False and self.dialect == 'oneomics':
            print('FATAL ERROR! The data as normalized by OneOmics is presently not being handled correctly.')
            print('             Set -newnorm cmd flag to renormalize here or find the code problem.')
            quit()

        return




    def get_mod_seq_col(self):
        """
        This determines which column has the modified sequence data and also what mod syntax is used if it can.
        """

        # this depends on both dialect and library vs. results - ex spectronaut col headers are different
        if self.level == 'transition' or self.level == 'peptide':

            # EXTENSION_ZONE

            if self.dialect == 'spectronaut':
                if self.lib_or_results == 'library':
                    self.modified_seq_col = 'LabeledPeptide'
                else:
                    self.modified_seq_col = 'FG.LabeledSequence'

                # spectronaut supports different syntax, so determine by inspection
                self.mod_syntax = modifications.determine_mod_syntax(self.df[self.modified_seq_col])

            elif self.dialect == 'sciex':
                if self.lib_or_results == 'library':
                    self.modified_seq_col = 'modification_sequence'
                    self.mod_syntax = 'TLC'
                else:
                    self.modified_seq_col = 'modification_sequence'  # todo - is this the same?
            elif self.dialect == 'skyline':
                self.modified_seq_col = 'Peptide Modified Sequence Three Letter Codes'
                self.mod_syntax = 'TLC'
            elif self.dialect == 'std':
                self.modified_seq_col = 'seq_w_mods'
            elif self.dialect == 'oneomics':
                self.modified_seq_col = 'Peptide HyperTag'
                self.mod_syntax = 'TLC'
            else:
                self.modified_seq_col = 'unknown'
                print('FATAL ERROR: modified sequence column needs to be specified for dialect:', self.dialect)
                quit()

        return




    def convert_format(self, **kwargs):
        """
        Given a dialect dataframe - library or results (stacked or unstacked), produce a dataframe in a new format.
        Can also convert results dataframes from stacked to unstacked - aka feature table.

        :param output_format: Target format may be specified - currently only 'skyline' for library conversion or 'std'
         for results conversion. These are implied automatically now, but kwarg set up for future option.
        :return:
        """

        output_format = kwargs.get('output_format', self.output_format)

        if self.verbose > 0:
            boiler_plates.section_header(title='Convert format', level=1)

        if output_format is None:
            if self.lib_or_results == 'library':
                output_format = 'skyline'
                # todo - if support library conversion to non-skyline, CHANGE HERE!!
            else:
                output_format = 'std'

        # validations
        if self.lib_or_results == 'results' and output_format != 'std':
            print('ERROR! SwathData class currently only supports conversion of results to "std" format.')
            quit()
        elif self.lib_or_results == 'library' and output_format != 'skyline':
            print('ERROR! SwathData class currently only supports conversion of libraries to "skyline" format.')
            quit()

        if self.dialect == 'std':
            print('ERROR! SwathData class currently does not support results input already in std format.')
            print('       Uses cases starting from std format should probably be built on pkl SwathData class input.')
            quit()


        # check that conversion has not already been called (can only convert once to be safe at this point)
        if self.converted is True:
            print('ERROR! SwathData class currently only supports calling convert_format method once.')
            quit()

        """
        Mapping instructions:

        # EXTENSION_ZONE

        - always include important columns, even if they don't get renamed or transformed (ex. the mod sequence)
        - Transformations are processed in the order they are listed.
        - If any changes are made to any float or int type metadata, they must be updated in 
        get_sample_columns() or what is passed to it as metadata cols.

        """

        results_mapping_dict = {

            'spectronaut to std':
                {'R.Label': {'new_name': 'file_name', 'transformations': [strip_wiff]},
                 'R.Condition': {'new_name': 'run_name'},
                 'PG.ProteinGroups': {'new_name': 'protein'},   # was 'protein' before, but need species in Skyline
                 'PG.ProteinNames': {'new_name': 'protein_names'},  # was 'protein_names' before
                 'PEP.StrippedSequence': {'new_name': 'sequence'},
                 'PEP.NrOfMissedCleavages': {'new_name': 'missed_cleavages'},
                 # 'EG.ModifiedPeptide':
                 #    {'new_name': 'seq_w_mods_wo_labels', 'transformations': [strip_flanking]},
                 'FG.LabeledSequence':
                     {'new_name': 'seq_w_mods',
                      'transformations': [strip_flanking]},
                 'FG.Charge': {'new_name': 'prec_z'},
                 'FG.PrecMz': {'new_name': 'prec_mz'},
                 'F.Charge': {'new_name': 'frag_z'},
                 'F.FrgMz': {'new_name': 'frag_mz'},
                 'F.FrgType': {'new_name': 'frag_type'},
                 'F.FrgLossType': {'new_name': 'frag_loss_type', 'transformations': [spectronaut_nl_frag_map]},
                 'F.FrgNum': {'new_name': 'frag_index'},
                 'F.PeakArea': {'new_name': 'intensity'},
                 #'PEP.Quantity': {'new_name': 'pep_intensity'}, # todo - uncomment when hack fixed
                 'EG.Quantity': {'new_name': 'pep_intensity'},  # todo - change to = 'prec_intensity'
                 'PG.Quantity': {'new_name': 'prot_intensity'},
                 'EG.Qvalue': {'new_name': 'qval'},
                 'EG.ApexRT': {'new_name': 'rt_obs'},
                 'EG.RTPredicted': {'new_name': 'rt_pred'},
                 'F.ExcludedFromQuantification': {'new_name': 'excluded_quantification'},
                 },

            'skyline to std':
                {'Protein Locator': {'new_name': 'protein',
                                     'transformations': [strip_mol_group_tag]},
                 'Peptide Sequence': {'new_name': 'sequence'},
                 # todo - need to add column(s) that may contain species signatures here
                 # '<not offered>':
                 #    {'new_name': 'seq_w_mods_wo_labels', 'transformations': [strip_flanking]},
                 'Peptide Modified Sequence Three Letter Codes':
                     {'new_name': 'seq_w_mods', 'transformations': [strip_nl_mods]},
                 'Precursor Charge': {'new_name': 'prec_z'},
                 'Precursor Mz': {'new_name': 'prec_mz'},
                 'Product Charge': {'new_name': 'frag_z'},
                 'Product Mz': {'new_name': 'frag_mz'},
                 'Fragment Ion Ordinal': {'new_name': 'frag_index'},  # ex '6'
                 'Fragment Ion Type': {'new_name': 'frag_type'},  # ex 'y'
                 'Loss Formulas': {'new_name': 'frag_loss_type', 'transformations': [string_conv_blank_map, skyline_nl_frag_map]},
                 'Fragment Ion': {'new_name': 'frag_type_index'},  # ex 'y6'
                 'Is Decoy': {'new_name': 'decoy'}
                 },

            'sciex to std':     # refers to desktop SWATH PeakView plugin
                {'Protein': {'new_name': 'protein'},
                 #'<not offered>': {'new_name': 'sequence'},    # todo - need to derive this from seq_w_mods for IDs!
                 'Peptide':
                     {'new_name': 'seq_w_mods'},
                 #'RT': {'new_name'}     # not clear how to map this todo - sort this out
                 # todo - need to add column(s) that may contain species signatures here
                 'Precursor Charge': {'new_name': 'prec_z'},
                 'Precursor MZ': {'new_name': 'prec_mz'},
                 'Fragment Charge': {'new_name': 'frag_z'},
                 'Fragment MZ': {'new_name': 'frag_mz'},
                 'Residue': {'new_name': 'frag_index'},  # ex '6'
                 'Ion Type': {'new_name': 'frag_type', 'transformations': [sciex_frag_type_map]}
                # will be 'y - 18' (includes NL, int syntax) before conv
                # always exported as feature table, so no single intensity column to map
                 },

            'oneomics to std':  # refers to .fcresults bundle generated by Sciex OneOmics cloud processing
                {'Protein Id': {'new_name': 'protein'},
                 'Peptide HyperTag': {'new_name': 'seq_w_mods'},
                 'Peptide Sequence': {'new_name': 'sequence'},
                 'Peptide Charge': {'new_name': 'prec_z'},
                 'Peptide Mz': {'new_name': 'prec_mz'},
                 'Transition Charge': {'new_name': 'frag_z'},
                 'Transition Mz': {'new_name': 'frag_mz'},
                 },

        }

        library_mapping_dict = {

            # todo - any species needs impact to these? Don't think so

            'sciex to skyline':
                {'Q1': {'new_name': 'Precursor m/z', 'required': True, 'std_name': 'prec_mz'},
                 'Q3': {'new_name': 'Product m/z', 'required': True, 'std_name': 'frag_mz'},
                 'uniprot_id': {'new_name': 'Protein', 'required': True, 'std_name': 'protein'},
                 'modification_sequence':
                     {'new_name': 'modification_sequence', 'required': True, 'std_name': 'seq_w_mods'},

                 'iRT': {'new_name': 'iRT', 'required': True, 'std_name': 'iRT'},
                 # todo - is there a relative intensity column to add and require?

                 'stripped_sequence': {'remove'},
                 },

            'spectronaut to skyline':
                {'PrecursorMz': {'new_name': 'Precursor m/z', 'required': True, 'std_name': 'prec_mz'},
                 'PrecursorCharge': {'new_name': 'precursorcharge', 'required': True, 'std_name': 'prec_z'},

                 'iRTSourceSpecific': {'new_name': 'irt', 'required': True, 'std_name': 'iRT'},
                 'RelativeIntensity':
                     {'new_name': 'library_intensity', 'required': True, 'std_name': 'library_intensity'},
                 # require intensity for all scores via mProphet

                 'FragmentMz': {'new_name': 'Product m/z', 'required': True, 'std_name': 'frag_mz'},
                 'FragmentCharge': {'new_name': 'productcharge', 'required': False, 'std_name': 'frag_z'},

                 # these are ignored by Skyline - just included to help verification
                 'FragmentType': {'new_name': 'Fragment Ion Type', 'required': False, 'std_name': 'frag_type'},
                 'FragmentNumber': {'new_name': 'Fragment Ion Ordinal', 'required': False, 'std_name': 'frag_index'},
                 'FragmentLossType': {'new_name': 'FragmentLossType', 'required': False},

                 #'UniProtIds': {'new_name': 'Protein', 'required': True, 'std_name': 'protein'},   # was the protein col before
                 'Protein Name': {'new_name': 'proteinname', 'required': True, 'std_name': 'protein'}, # was 'protein_names'

                 'LabeledPeptide': {'new_name': 'Peptide Modified Sequence Three Letter Codes',
                                    'transformations': [strip_flanking],
                                    'required': True, 'std_name': 'seq_w_mods'},
                 'StrippedPeptide': {'remove'},
                 'ModifiedPeptide': {'remove'},
                 'IntModifiedPeptide': {'remove'},
                 'IntLabeledPeptide': {'remove'},

                 'ExcludeFromAssay': {'new_name': 'exclude_transition'},

                 },
        }

        conversion = self.input_format + ' to ' + output_format

        # stop if conversion not supported
        supported_conversions = {'spectronaut to skyline', 'sciex to skyline',              # library conversions

                                 'spectronaut to std', 'sciex to std', 'skyline to std',    # results conversions
                                 'oneomics to std'
                                 }

        if self.lib_or_results == 'results' and 'std' not in conversion or \
                conversion not in supported_conversions:
            print('FATAL ERROR! Unsupported conversion:', conversion)
            quit()



        # ---------------------------------
        # select appropriate mapping dictionary

        if self.lib_or_results == 'results':
            mapping_dict = results_mapping_dict[conversion]
        else:
            mapping_dict = library_mapping_dict[conversion]



        # ---------------------------------
        # Determine which columns to keep from the initial dataframe

        # determine which columns can be mapped to items in the mapping_dict
        available_mapping = {x: mapping_dict[x] for x in mapping_dict if x in self.df.columns}

        # get list of columns explicitly marked to remove in mapping_dict
        remove_cols = [x for x in mapping_dict if 'remove' in mapping_dict[x]]

        # remove these remove columns from available mapping so don't deal with them after this
        available_mapping = {x: available_mapping[x] for x in available_mapping if x not in remove_cols}


        # get quant columns if possible (the called method checks if possible)
        self.get_sample_columns()


        # exclusive (if not explicitly listed, removed at this point)
        final_cols = [x for x in available_mapping if x not in remove_cols]
        # todo - previously, exclusion of columns not explicitly included was conditional - kept for libraries, but skyline doesn't carry extra stuff
        # todo - cont - the only reason to keep stuff is verification, and should explicitly add anything we want for that.

        # add sample data columns if have them at this point (NOT derived ones like delta rts! Those are added later)
        for x in [self.intensity_columns, self.score_columns, self.q_value_columns,
                  self.rt_columns, self.pred_rt_columns]:
            if len(x) > 0:
                final_cols += x


        # ---------------------------------
        # KEY STEP - reduce to desired subset of columns

        self.df = self.df[final_cols]


        # ---------------------------------
        # Convert all the quant cols from int64 to float64 if needed

        if self.dialect == 'skyline' and len(self.intensity_columns) > 0:  # is results data if have sample cols
            # Convert all the quant cols from int64 to float64, setting zeros to nan
            self.df[self.intensity_columns] = self.df[self.intensity_columns].astype('float64').replace(0, np.nan)

        # (this is already done upon read for OneOmics data)


        # ---------------------------------
        # Ensure types of potential problematic columns

        if self.verbose > 2:
            print()
            print('Types after:')
            self.get_dtypes()
            print()


        # ---------------------------------
        # Drop rows w missing metadata (incomplete description of analyte cannot be used)

        if 'FG.LabeledSequence' in self.df.columns:     # observed some NAN floats in seq col w Spectronaut
            shape_before = self.df.shape
            self.df.dropna(subset=['FG.LabeledSequence'], inplace=True)
            shape_after = self.df.shape
            if self.verbose > 0 and shape_before != shape_after:
                boiler_plates.section_header(title='Removal of rows with missing values in metadata columns', level=2)
                print('   shape before:', shape_before)
                print('   shape after: ', shape_after)
                print()



        # ---------------------------------
        # Transformations

        for col in available_mapping:

            if 'transformations' in available_mapping[col]:

                for trans in available_mapping[col]['transformations']:
                    self.df[col] = self.df[col].apply(trans)


        # ---------------------------------
        # Convert column names

        renaming_dict = {x: available_mapping[x]['new_name'] for x in available_mapping}

        self.df.rename(columns=renaming_dict, inplace=True)

        # update modified seq col
        self.modified_seq_col = available_mapping[self.modified_seq_col]['new_name']

        # also rename metadata cols if high levels already present
        if self.pep_df is not None:
            self.pep_df.rename(columns=renaming_dict, inplace=True)

        if self.prot_df is not None:
            self.prot_df.rename(columns=renaming_dict, inplace=True)


        """ ============================================================================================

            KEY BOUNDARY POINT IN THE CODE!!! - Reference by post-conversion names hereafter!!

            ============================================================================================
        """

        # ---------------------------------
        # Mod syntax standardization

        if self.mod_syntax == 'PSI-Mod':
            self.df[self.modified_seq_col] = self.df[self.modified_seq_col].apply(self.long_mod_to_tlc)

            if self.pep_df is not None:
                self.pep_df[self.modified_seq_col] = self.pep_df[self.modified_seq_col].apply(self.long_mod_to_tlc)

            # todo - add prec_df support here or pep_df

            self.mod_syntax = 'TLC'

        elif self.mod_syntax == 'integer':
            print('FATAL ERROR! Conversion of integer mod syntax to TLC is not supported yet!')
            quit()


        # ---------------------------------
        # Frag type and loss type standardization

        if 'frag_loss_type' in self.df.columns:

            # if Spectronaut results, merge the frag type and loss type to one frag type field
            #if self.dialect == 'spectronaut':

            # merge frag_type and frag_loss_type to single field frag_type (already together in sciex)
            # this assumes that the mapping has already been applied to the loss field in spectronaut_nl_frag_map()
            self.df['frag_type'] = self.df['frag_type'].astype(str) + self.df['frag_loss_type'].astype(str)
            self.df.drop(['frag_loss_type'], axis=1, inplace=True)


        # ---------------------------------
        # Get species info for mixed species spike type experiment

        if self.lib_or_results == 'results' and self.expt_type == 'dilution_series_species_spike':

            # Validate have necessary info
            if self.const_species is None or self.var_species is None or self.species_column is None:
                print('FATAL ERROR!! When experiment type = dilution_series_species_spike, you must provide')
                print('              the species_column, var_species, and const_species parameters.')
                quit()

            new_species_col = mapping_dict[self.species_column]['new_name']

            self.df['spike_status'] = \
                self.df.apply(lambda x: categorize_species_spike_status(x[new_species_col],
                                                                        self.const_species,
                                                                        self.var_species), axis=1)

            if self.pep_df is not None:
                self.pep_df['spike_status'] = \
                    self.pep_df.apply(lambda x: categorize_species_spike_status(x[new_species_col],
                                                                                self.const_species,
                                                                                self.var_species), axis=1)

            if self.prot_df is not None:
                self.prot_df['spike_status'] = \
                    self.prot_df.apply(lambda x: categorize_species_spike_status(x[new_species_col],
                                                                                 self.const_species,
                                                                                 self.var_species), axis=1)

        # ---------------------------------
        # Set column types

        # todo - seq w mods was getting interpreted as non string in a concat - any others to force type on?
        if self.lib_or_results == 'results':

            # set all the integer metadata to this type (for some reason read as floats sometimes)
            if col in self.df.columns:
                for col in ['prec_z', 'frag_z', 'frag_index']:
                    self.df[col] = self.df[col].astype('int')


        # ---------------------------------
        # Unstack if needed or get quant columns if already in analyte x samples format

        # todo - add other ways to recognize when need to unstack data - ex some skyline exports look like this
        # EXTENSION_ZONE
        if self.lib_or_results == 'results':

            if self.dialect == 'spectronaut':

                self.unstack()


        # ---------------------------------
        # Create semantic lookup dict - where to find information

        self.get_semantic_dict(mapping_dict)


        # ---------------------------------
        # Validate library

        if self.lib_or_results == 'library':

            self.validate_library()


        # ---------------------------------
        # Update attributes

        # change the present dialect to the converted format
        self.dialect = output_format

        # change this flag to True so can prevent calling convert() more than once
        self.converted = True

        # update col w species info
        if self.species_column is not None:
            self.species_column = mapping_dict[self.species_column]['new_name']


        # ---------------------------------
        # Create an ID for each analyte in all data levels present, setting the appropriate one as the df index

        self.get_indices()

        # ---------------------------------
        # Rollup and complete the data at all levels

        if self.lib_or_results == 'results':

            # if we don't have higher levels, compute them now that we have indices
            if self.pep_df is None or self.ignore_provided_rollup is True:

                # todo - make sure this logic works for skyline case trying to use higher level data provided (not supported yet)
                self.rollup(method='sum')

            # Add derived data, regardless of how we got the complete data above
            self.add_derived_data_arrays()


        # ---------------------------------
        # Characterize analytes (this is the first point we can do this since uses indices

        self.count_analytes('After conversion')


        self.tracker.capture('convert format')

        return




    def get_semantic_dict(self, mapping_dict):

        # build a semantic mapping dict
        # (inverting the mapping dict so can find data knowing the std name for semantic content)

        for key in mapping_dict:
            if 'std_name' in mapping_dict[key] or self.lib_or_results == 'results':

                if self.lib_or_results == 'results':
                    std_name = mapping_dict[key]['new_name']
                    column = std_name
                else:
                    std_name = mapping_dict[key]['std_name']
                    column = mapping_dict[key]['new_name']

                required = mapping_dict[key].get('required', False)

                self.semantic_dict.update({std_name: {'column': column,
                                                      'required': required}})

        return




    def validate_library(self):

        initial_rows = self.df.shape[0]

        # get list of required columns using std_name for semantic content per library mapping dict
        required_columns = [key for key in self.semantic_dict if self.semantic_dict[key].get('required', False) is True]

        # trackers
        missing_counts = []
        current_cols_w_gaps = []
        col_completely_missing = []

        # pass through just counting
        for std_name in required_columns:

            current_name = self.semantic_dict[std_name]['column']

            if current_name in self.df.columns:

                # cont rows where missing values for this key column
                values_present = self.df[current_name].count()
                missing_values = initial_rows - values_present
                missing_counts.append(missing_values)

                if missing_values > 0:
                    current_cols_w_gaps.append(current_name)

            else:
                col_completely_missing.append(current_name)
                missing_counts.append(initial_rows)


        # remove rows where nan in required col
        self.df.dropna(subset=current_cols_w_gaps, inplace=True)


        if self.verbose > 0:
            boiler_plates.section_header(title='Check for missing required content in library', level=2)
            print('Missing values in required columns:')
            print()
            for i in range(len(required_columns)):
                print('    ', required_columns[i], ':   ', missing_counts[i])
            if len(current_cols_w_gaps) > 0:
                before_after_row_output(initial_rows, self.df.shape[0])
            print()
            if len(col_completely_missing) > 0:
                print()
                print('FATAL ERROR! The following column(s) are required and not present:')
                for x in col_completely_missing:
                    print('    ', x)
                print()
                quit()

        return




    def get_indices(self):


        # Get column headers in current dialect for semantic content necessary to build indices

        protein = self.semantic_dict['protein'].get('column', 'protein')
        #sequence = self.semantic_dict['sequence'].get('column', 'sequence')
        seq_w_mods = self.semantic_dict['seq_w_mods'].get('column', 'seq_w_mods')
        prec_mz = self.semantic_dict['prec_mz'].get('column', 'prec_mz')
        frag_mz = self.semantic_dict['frag_mz'].get('column', 'frag_mz')
        prec_z = self.semantic_dict['prec_z'].get('column', 'prec_z')
        frag_z = self.semantic_dict['frag_z'].get('column', 'frag_z')
        if 'frag_type' in self.semantic_dict:
            frag_type = self.semantic_dict['frag_type'].get('column', 'frag_type')
        if 'frag_index' in self.semantic_dict:
            frag_index = self.semantic_dict['frag_index'].get('column', 'frag_index')


        # EXTENSION_ZONE

        # transition dataframe
        if self.df is not None:

            # transition indexing system 1 - have fragment indices so can uniquely indicate trans without masses
            if 'frag_index' in self.df.columns:

                # note: frag_type and frag_loss_type are combined prior to here
                #       they are combined like y-H20 so that no NL gives just 'y_x', not y__xx
                #       and more importantly, so that the ID is the same vs. data that does not use NL fragments

                self.df['trans_ID'] = \
                    self.df[protein].astype(str) + '_' + self.df[seq_w_mods].astype(str) + '_' + \
                    self.df[prec_z].astype(str) + '_' + \
                    self.df[frag_type].astype(str) + '_' + self.df[frag_index].astype(str) + '_' + \
                    self.df[frag_z].astype(str)

                # if change ID format later and need to add: # + self.df['frag_loss_type'].astype(str) + \

                # also capture the m/z-based transition ID to enable comparison across systems
                self.df['trans_ID_mz'] = \
                    self.df[protein].astype(str) + '_' + '_' + self.df[seq_w_mods].astype(str) + '_' + \
                    + self.df[prec_mz].round(4).astype(str) + '_' + self.df[frag_mz].round(4).astype(str)

            # transition indexing system 2 where have no choice but use masses to indicate transition
            else:
                self.df['trans_ID'] = \
                    self.df[protein].astype(str) + '_' + '_' + self.df[seq_w_mods].astype(str) + '_' + \
                    + self.df[prec_mz].round(4).astype(str) + '_' + self.df[frag_mz].round(4).astype(str)

                self.mz_indexing = True

            self.df.set_index('trans_ID', drop=False, inplace=True)

            # add the precursor index
            self.df['prec_ID'] = \
                self.df[protein].astype(str) + '_' + self.df[seq_w_mods].astype(str) + '_' + self.df[prec_z].astype(str)

            # add the peptide index
            self.df['pep_ID'] = self.df[protein].astype(str) + '_' + self.df[seq_w_mods].astype(str)

            # add the protein index
            self.df['prot_ID'] = \
                self.df[protein].astype(str)

        # peptide dataframe
        if self.pep_df is not None:

            # todo - resolve this temporary hack (this will screw up counts of pep per protein)

            # TEMP HACK - SET THE PEPTIDE ID AS THE PREC ID SINCE THAT'S WHAT THE TABLE REALLY IS
            self.pep_df['pep_ID'] = \
                self.pep_df[protein].astype(str) + '_' + self.pep_df[seq_w_mods].astype(str) + '_' + \
                self.pep_df[prec_z].astype(str)

            """  GO BACK TO THIS ONCE FIXED
            self.pep_df['pep_ID'] = \
                self.pep_df[protein].astype(str) + '_' + self.pep_df[seq_w_mods].astype(str)
            """

            self.pep_df.set_index('pep_ID', drop=False, inplace=True)

            # always compute the prot_ID - may not be synonymous with single col in future
            self.pep_df['prot_ID'] = self.pep_df[protein].astype(str)

        # protein dataframe
        if self.prot_df is not None:

            self.prot_df['prot_ID'] = self.prot_df[protein].astype(str)

            self.prot_df.set_index('prot_ID', drop=False, inplace=True)


        # --------------------------------------
        # Check all index values are unique at all levels present

        if self.verbose > 0:
            boiler_plates.section_header(title='Generate Unique Analyte Identifiers', level=2)

        for df, level in [(self.df, 'transition'), (self.pep_df, 'peptide'), (self.prot_df, 'protein')]:

            if df is not None:

                index_unique = df.index.nunique()

                if self.verbose > 0:
                    print()
                    print('Verify index identifiers are unique at level:', level)
                    print()
                    print('   Unique index values:   ', index_unique)
                    print('   Rows in table:         ', df.shape[0])
                    print()
                    if df.shape[0] == index_unique:
                        print('CONFIRMED. All index values are unique.')


                # If dups detected, first try dropping identical duplicates (in all columns)
                if df.shape[0] != index_unique:

                    rows_before_identical_dup_drop = df.shape[0]

                    df.drop_duplicates(inplace=True, keep='first') # uses all columns when subset= not specified

                    if self.verbose > 0:
                        print()
                        print('Dropping identical rows (identical in all columns)')
                        print()
                        print('   Unique index values:   ', index_unique)
                        print('   Rows in table after:   ', df.shape[0])
                        if df.shape[0] < rows_before_identical_dup_drop:
                            print()
                            print('WARNING!!!!!!! Identical rows were detected and dropped!!')
                            print('rows dropped (keeping first dup):  ',
                                  rows_before_identical_dup_drop - df.shape[0])


                # If dups still left, flag and crash
                if index_unique < df.shape[0]:

                    # get a new df with duplicated index rows, keeping all instances of duplicates to find the problem
                    dups_df = df[df.index.duplicated(keep=False)].copy(deep=True)

                    # sort by index to group the dups together in output
                    dups_df.sort_index(inplace=True)

                    print()
                    print('ERROR - DUPLICATE INDEX ROWS!!!! This should never happen - sort out why in output below!')
                    print()
                    # write out the dups with just a subset of columns to inspect to find the cause
                    print(tabulate(dups_df[['prec_mz', 'prec_z',
                                            'frag_mz', 'frag_z',
                                            'frag_type', 'frag_index']],  # 'frag_loss_type'
                                   headers='keys', tablefmt='simple'))
                    quit()

        return




    def get_sample_columns(self, **kwargs):
        """
        The standardized results format contains two kinds of columns:
        - metadata describing the analyte monitored in a given row
        - the names of each sample or run that the quantitative measurements pertain to in that column, of which
        there can be multiple data types - intensity, score, rt, etc

        This function gets a lists of the latter by type so the quantitative matrices within the df can be analyzed.

        This function can also be used on non-standardized df's to find which columns are quant columns. In this case,
        it must be passed a metadata_cols list telling which non-standard cols of the same type are metadata, not quant
        values.

        :return:
        """

        type_for_quant_data = kwargs.get('type_for_quant_data', np.float64)


        if self.lib_or_results != 'results' or self.stacked is True or self.dialect == 'unknown':
            return

        if self.verbose > 0:
            boiler_plates.section_header(title='Determine Sample Data Columns', level=2)

        # remove columns containing exclude signature string # todo - could be too brute force doing this here?
        if len(self.exclude_sample_strings) > 0:

            cols_to_exclude = []

            types_dict = self.get_dtypes()

            for col in self.df.columns:
                match_count = \
                    sum([1 for x in self.exclude_sample_strings if x in col]) # and types_dict[col] == np.float64])
                if match_count > 0:
                    cols_to_exclude.append(col)

            if len(cols_to_exclude) > 0:

                if self.verbose > 0:
                    print()
                    print('Excluding', len(cols_to_exclude), 'columns with exclude strings:',
                          self.exclude_sample_strings)

                    for x in cols_to_exclude:
                        print('   ', x)

                self.df.drop(cols_to_exclude, axis=1, inplace=True)

                if self.pep_df is not None:
                    self.pep_df.drop(cols_to_exclude, axis=1, inplace=True, errors='ignore')

                if self.prot_df is not None:
                    self.prot_df.drop(cols_to_exclude, axis=1, inplace=True, errors='ignore')

        # EXTENSION_ZONE
        data_array_mapping_dict = {
            'std':
                {self.intensity_suffix: self.intensity_suffix,
                 self.score_suffix: self.score_suffix,
                 self.obs_rt_suffix: self.obs_rt_suffix,
                 self.pred_rt_suffix: self.pred_rt_suffix,
                 self.q_value_suffix: self.q_value_suffix
                 },
            'skyline':
                {' Area': self.intensity_suffix,
                 ' Library Dot Product': self.score_suffix,
                 ' Best Retention Time': self.obs_rt_suffix,
                 ' Predicted Result Retention Time': self.pred_rt_suffix,
                 ' Detection Q Value': self.q_value_suffix
                 },
            'spectronaut':
                {' Intensity': self.intensity_suffix,
                 ' Score': self.score_suffix,
                 ' Obs_RT': self.obs_rt_suffix,
                 ' Pred_RT': self.pred_rt_suffix,
                 ' q-value': self.q_value_suffix
                 },
            'oneomics':  # will already be converted by this point on import
                {self.intensity_suffix: self.intensity_suffix,
                 self.score_suffix: self.score_suffix,
                 self.obs_rt_suffix: self.obs_rt_suffix,
                 self.pred_rt_suffix: self.pred_rt_suffix,
                 self.q_value_suffix: self.q_value_suffix
                 },
        }


        # throw error if don't have a mapping dict for the results type
        if self.dialect not in data_array_mapping_dict:
            map_dict = {}
            print('FATAL ERROR! No column content mapping dict for dialect:', self.dialect)
            print('             in SwathData method get_sample_columns().')
            quit()

        else:

            map_dict = data_array_mapping_dict[self.dialect]

            dialect_suffix_list = list(map_dict.keys())

            # convert all sample-related data cols to std names
            renaming_dict = {}

            for col in self.df.columns:

                known_suffix = [x for x in dialect_suffix_list if col.endswith(x)]

                if len(known_suffix) > 0:
                    renaming_dict.update({col: col.replace(known_suffix[0], map_dict[known_suffix[0]])})


            if len(renaming_dict) > 0:

                self.df.rename(columns=renaming_dict, inplace=True)

                if self.pep_df is not None:
                    self.pep_df.rename(columns=renaming_dict, inplace=True)

                if self.pep_df is not None:
                    self.prot_df.rename(columns=renaming_dict, inplace=True)


            self.intensity_columns = [x for x in self.df.columns if x.endswith(self.intensity_suffix)]
            self.score_columns = [x for x in self.df.columns if x.endswith(self.score_suffix)]
            self.q_value_columns = [x for x in self.df.columns if x.endswith(self.q_value_suffix)]
            self.rt_columns = [x for x in self.df.columns if x.endswith(self.obs_rt_suffix)]
            self.pred_rt_columns = [x for x in self.df.columns if x.endswith(self.pred_rt_suffix)]
            self.sample_name_stems = [x.replace(self.intensity_suffix, '') for x in self.intensity_columns]

            # add the delta rt columns if can compute them
            # note that the actual data columns are not added to the dataframes until later, but populating the column
            # list must be done here and we can tell if it will be possible to compute later them here.
            if len(self.pred_rt_columns) == len(self.rt_columns) and len(self.rt_columns) > 0:

                self.delta_rt_columns = \
                    [x.replace(self.pred_rt_suffix, self.delta_rt_suffix) for x in self.pred_rt_columns]

            # define a collection of these now that some are not empty
            self.data_arrays_list = [('intensity', self.intensity_columns, self.intensity_suffix),
                                     ('q_val', self.q_value_columns, self.q_value_suffix),
                                     ('score', self.score_columns, self.score_suffix),
                                     ('obs_rt', self.rt_columns, self.obs_rt_suffix),
                                     ('pred_rt', self.pred_rt_columns, self.pred_rt_suffix),
                                     ('delta_rt', self.delta_rt_columns, self.delta_rt_suffix)]


        if self.verbose > 0:
            print('Sample data arrays detected and converted from original dialect:', self.dialect, 'to standard.')

            for (data_type, column_list, suffix) in self.data_arrays_list:

                if len(column_list) > 0:
                    print()
                    print('    ', data_type, 'columns:', len(column_list), '     suffix: ', suffix)

                    if self.verbose > 0:
                        for x in column_list:
                            print('        ', x)


        # check the existence and ordering of parallel array data
        self.check_sample_data_col_alignment()


        # initialize sample norm factors to 1.0 for all samples (now that we know the sample names)
        self.sample_norm_factors = {sample: 1.0 for sample in self.intensity_columns}


        return




    def add_derived_data_arrays(self):

        # add the delta rt columns if can compute them
        if len(self.delta_rt_columns) > 0:

            for df in [self.df, self.pep_df]:

                if df is not None:

                    for col in self.delta_rt_columns:
                        df[col] = abs(df[col.replace(self.delta_rt_suffix, self.obs_rt_suffix)] -
                                      df[col.replace(self.delta_rt_suffix, self.pred_rt_suffix)])


        return




    def check_sample_data_col_alignment(self):
        """
        Make sure have corresponding data for each result data type present for all samples
        """

        # nothing to check if no intensity cols yet, so return
        if len(self.intensity_columns) == 0:
            return


        # first make sure the intensity cols are the issue by checking it has the max of any len
        if len(self.intensity_columns) < len(self.rt_columns) or \
            len(self.intensity_columns) < len(self.pred_rt_columns) or \
            len(self.intensity_columns) < len(self.q_value_columns):

            print('FATAL ERROR!!: The there are fewer intensity columns than other data types:')
            print('len intensities:', len(self.intensity_columns))
            print('len obs rts:    ', len(self.rt_columns))
            print('len pred rts:   ', len(self.pred_rt_columns))
            print('len q-values:   ', len(self.q_value_columns))
            # todo - solve this below and remove quit here. Change to warning
            quit()

        # if pass this first check, then check that there is a corresponding col in other types for each intensity col
        for (data_type, column_list, suffix) in self.data_arrays_list:

            # if have addition data array columns, make sure they are all in the same order as intensities
            # (name stem is derived from intensity list, so is necessarily in the same order already)
            # (this doesn't matter for pandas since label-based generally, but being safe to allow matrix ops in numpy)
            if len(column_list) > 0:
                count_unmatched = 0

                for s in range(len(self.intensity_columns)):

                    # check the analog is there anywhere (regardless of order)
                    missing = False
                    analog = self.intensity_columns[s].replace(self.intensity_suffix, suffix)
                    if analog not in column_list:
                        print('ERROR!! Missing analogous data:')
                        print('          intensity col:         ', self.intensity_columns[s])
                        print('          analogous', data_type, 'col:    ', analog)
                        print()
                        print('If you are seeing this error, a likely source is inclusion of blank runs.')
                        print('You can easily ignore blank runs using -xsamp arg from cmd line or passing a list to')
                        print('the kwarg "exclude_sample_strings" when instantiating the SwathData class.')
                        print("For example, if the word 'blank' is the file name for all blanks, pass ['blank']")
                        print('or cmd arg -xsamp "blank".')
                        print()
                        missing = True
                        quit()

                        # todo - track stems that are incomplete and remove them at end for all col types
                        # todo - OR instead insert the missing col w all nan's

                    if missing is False:
                        if column_list[s].replace(suffix, self.intensity_suffix) != self.intensity_columns[s]:
                            count_unmatched += 1

                if count_unmatched > 0:
                    print('ERROR!!', data_type, 'are not sorted same as intensity cols - NEED TO WRITE SORT CODE!')
                    quit()

        return




    def get_dtypes(self):

        types_dict = self.df.dtypes.to_dict()

        if self.verbose > 1:
            print()
            print('types for all columns in dataframe:')
            for x in types_dict:
                print('col:', x, '  type:', types_dict[x])

        return types_dict




    def filter_by_species(self):

        if len(self.include_species) == 0:
            return

        if self.verbose > 0:
            boiler_plates.section_header(title='Filter rows without required species', level=1)
            print('Species required strings:   ', self.include_species)

        self.df = \
            pandas_helpers.exclude_rows_without_substring_in_col(self.df,
                                                                 [(self.species_column, self.include_species)],
                                                                 title='Filter to Included Species: Transition Level',
                                                                 verbose=self.verbose)

        if self.prec_df is not None:
            self.prec_df = \
                pandas_helpers.exclude_rows_without_substring_in_col(self.prec_df,
                                                                     [(self.species_column, self.include_species)],
                                                                     title='Filter to Included Species: Peptide Level',
                                                                     verbose=self.verbose)

        if self.pep_df is not None:
            self.pep_df = \
                pandas_helpers.exclude_rows_without_substring_in_col(self.pep_df,
                                                                     [(self.species_column, self.include_species)],
                                                                     title='Filter to Included Species: Peptide Level',
                                                                     verbose=self.verbose)

        if self.prot_df is not None:
            self.prot_df = \
                pandas_helpers.exclude_rows_without_substring_in_col(self.prot_df,
                                                                     [(self.species_column, self.include_species)],
                                                                     title='Filter to Included Species: Protein Level',
                                                                     verbose=self.verbose)

        # get counts after filtering
        self.count_analytes('filter by species')

        return




    def filter_mods(self, **kwargs):
        """

        :param kwargs:
        filter_type: a string indicating a set of mods - ex. 'high' (currently only supported one)
        allowed_mods: pass a set in TLC format like 'M[Oxi]' and '[1Ac]-' to filter rows w mods not in this set.
        :return:
        """

        filter_type = kwargs.get('filter_type', 'high')
        allowed_mods = kwargs.get('allowed_mods', set())     # allow alternate route passing in specific list

        if self.verbose > 0:
            boiler_plates.section_header(title='Filter Rows with Unwanted Modifications', level=1)


        # check mods have been converted to right syntax
        if self.mod_syntax != 'TLC':
            print('Cannot filter mods yet as not converted to TLC syntax yet.')
            return

        # get mod_counts_dict if don't have already
        if self.mod_counts_dict is None:
            self.characterize_mods()

        # check have a mods counts dict now
        if self.mod_counts_dict is None:
            print('Mod frequencies cannot be computed yet.')
            # todo - it would be unnecessary to have this if I change how the removal step works
            return


        # if explicit set of mods passed in, use that
        if filter_type is None and len(allowed_mods) > 0:
            pass    # already have allowed_mods

        # otherwise an allowed set is indicated by a controlled list of names (currently just 'high')
        elif filter_type == 'high':
            allowed_mods = {'M[Oxi]', 'N[Dea]', 'Q[Dea]', 'C[CAM]', '[1Ac]-'}  # hardcoded here for now

        else:
            allowed_mods = set()
            print('FATAL ERROR!! filter_type argument not recognized:', filter_type)


        # hard-coded software specific logic here here!
        if self.dialect == 'skyline':   # self.dialect will be it's post-conversion value when this method called
            must_drop_mods = {'E[Dhy]'} # these aren't in skyline dictionary so case failure, not jus error
        else:
            must_drop_mods = set()

        # find mods present in but not allowed
        want_to_drop_mods = {m for m in self.mod_counts_dict if m not in allowed_mods}

        remove_mods_list = list(must_drop_mods | want_to_drop_mods)

        if self.verbose > 0:
            print('Allowed mods:')
            print(allowed_mods)
            print()
            print('Mods to be removed:')
            print(remove_mods_list)
            print()


        # remove the rows with these mods
        self.remove_if_contains(self.modified_seq_col, remove_mods_list, str)

        # get counts after filtering
        self.count_analytes('filter mods')

        self.tracker.capture('filter mods')

        return




    def filter_mz(self, **kwargs):

        min_ms1 = kwargs.get('min_ms1', None)
        max_ms1 = kwargs.get('max_ms1', None)
        min_ms2 = kwargs.get('min_ms2', None)
        max_ms2 = kwargs.get('max_ms2', None)

        if self.verbose > 0:
            boiler_plates.section_header(title='Filter Rows Outside m/z Ranges', level=1)

        prec_mz = self.semantic_dict['prec_mz'].get('column', 'prec_mz')
        frag_mz = self.semantic_dict['frag_mz'].get('column', 'frag_mz')

        if min_ms1 is not None:
            self.df.drop(self.df[self.df[prec_mz] < min_ms1].index, inplace=True)

            if self.prec_df is not None:
                self.prec_df.drop(self.prec_df[self.prec_df[prec_mz] < min_ms1].index, inplace=True)

        if max_ms1 is not None:
            self.df.drop(self.df[self.df[prec_mz] > max_ms1].index, inplace=True)

            if self.prec_df is not None:
                self.prec_df.drop(self.prec_df[self.prec_df[prec_mz] < max_ms1].index, inplace=True)

        if min_ms2 is not None:
            self.df.drop(self.df[self.df[frag_mz] < min_ms2].index, inplace=True)

        if max_ms2 is not None:
            self.df.drop(self.df[self.df[frag_mz] > max_ms2].index, inplace=True)



        # get counts after filtering
        self.count_analytes('filter by m/z ranges')

        self.tracker.capture('filter by m/z ranges')

        return




    def characterize_mods(self):
        """
        This is functional, not just descriptive.

        :return:
        """

        if self.mod_syntax == 'TLC':

            if self.pep_df is not None:
                self.mod_counts_dict = modifications.get_mod_freqs(self.pep_df,
                                                                   self.modified_seq_col, 'TLC',
                                                                   title='Peptide Level Modification Counts',
                                                                   verbose=self.verbose)
            elif self.df is not None:
                self.mod_counts_dict = modifications.get_mod_freqs(self.df,
                                                                   self.modified_seq_col, 'TLC',
                                                                   title='Transition Level Modification Counts',
                                                                   verbose=self.verbose)

        return




    def filter_excluded(self):
        """
        This method removes rows that have 'True' in any column indicating to exclude from the assay. Rather than
        trying to sort out which dialect is in use, a list simply includes the names of all such columns in any dialect
        and the method looks for their presence.

        :return:
        """

        # EXTENSION_ZONE
        headers_w_exclude_flags = ['ExcludeFromAssay', 'exclude_transition', 'F.ExcludedFromQuantification', 'excluded_quantification']

        before_rows = self.df.shape[0]

        if self.verbose > 0:
            boiler_plates.section_header(title='Filter Rows with Excluded Indicators from Analysis Software', level=1)

        # remove the rows with these flags
        for col in headers_w_exclude_flags:
            if col in self.df.columns:
                self.remove_if_contains(col, [True], bool, verbose=0)   # suppress output here and give net below

        if self.verbose > 0:
            before_after_row_output(before_rows, self.df.shape[0])


        # get counts after filtering
        self.count_analytes('filter excluded items')


        self.tracker.capture('filter excluded rows')

        return




    def filter_swath_decoys(self):
        """
        This method removes rows that are SWATH decoys (not decoy proteins!) from standard results.

        :return:
        """

        if self.dialect != 'std' or self.lib_or_results != 'results':
            return

        # EXTENSION_ZONE
        headers_w_exclude_flags = ['decoy']

        before_rows = self.df.shape[0]

        if self.verbose > 0:
            boiler_plates.section_header(title='Filter SWATH decoy rows', level=1)

        # remove the rows with these flags
        for col in headers_w_exclude_flags:
            if col in self.df.columns:
                self.remove_if_contains(col, [True], bool, verbose=0)   # suppress output here and give net below

        if self.verbose > 0:
            before_after_row_output(before_rows, self.df.shape[0])


        # get counts after filtering
        self.count_analytes('filter swath decoy rows')


        self.tracker.capture('filter swath decoy rows')

        return




    def output(self, outpath):

        if self.verbose > 0:
            boiler_plates.section_header(title='Output text file table(s)', level=1)


        if self.name_stem is not None:
            if self.lib_or_results == 'results':
                outfile = os.path.join(outpath, self.name_stem + '_TRANS.txt')
            else:
                outfile = os.path.join(outpath, self.name_stem + '_LIB.txt')

        # Don't include the index if this is a library, as never used on the front end, but include if results
        if self.lib_or_results == 'results':
            include_index = True
            drop_cols = []
        else:
            include_index = False
            # columns generated but not wanted in library output
            drop_cols = ['exclude_transition',
                         'trans_ID', 'prec_ID', 'pep_ID', 'prot_ID',
                         'trans_per_prec', 'pep_per_prot']

        keep_cols = [x for x in self.df.columns if x not in drop_cols]
        self.df[keep_cols].to_csv(outfile, sep='\t', index=include_index)

        if self.pep_df is not None:

            keep_cols = [x for x in self.pep_df.columns if x not in drop_cols]
            self.pep_df[keep_cols].to_csv(outfile.replace('_TRANS.txt', '_PEP.txt'), sep='\t', index=include_index)

            keep_cols = [x for x in self.prot_df.columns if x not in drop_cols]
            self.prot_df[keep_cols].to_csv(outfile.replace('_TRANS.txt', '_PROT.txt'), sep='\t', index=include_index)

        self.tracker.capture('output')

        return




    def save(self, outfile):

        # make sure extension is pkl
        if '.pkl' not in outfile:
            outfile = outfile.replace('.txt', '').replace('.csv', '').replace('.xls', '').replace('.xlsx', '').replace(
                '.xls', '') + '.pkl'

        # dump the object to pkl
        helper_functions.save_object(self, outfile)

        self.tracker.capture('save')

        return




    def profile(self):

        self.tracker.report(title='Timing Profile for SwathData Class', verbose=self.verbose)

        return




    def unstack(self):
        """
        Takes stacked quant data and unstacks it to a feature table. The input MUST be a standardized format results df.

        :return:
        """

        if self.verbose > 0:
            boiler_plates.section_header(title='Unstacking Sample Data Array(s)', level=2)
            print('Columns before unstacking:  ', self.df.shape[1])
            print()

        # Need to capture columns before to determine which are new sample columns added by unstacking
        before_cols = self.df.columns.tolist()

        # EXTENSION_ZONE (intentionally a list at trans level since filter based on presence next)
        trans_metadata = ['protein', 'protein_names',
                          'sequence', 'seq_w_mods', 'missed_cleavages',
                          'prec_mz', 'prec_z',
                          'frag_mz', 'frag_z', 'frag_index', 'frag_type', 'frag_loss_type',
                          'spike_status', 'excluded_quantification']

        pep_metadata = ['protein', 'protein_names',
                        'sequence', 'seq_w_mods', 'missed_cleavages',
                        'spike_status']

        prot_metadata = ['protein', 'protein_names', 'spike_status']

        df_tuples = [(self.prot_df, 'prot_intensity', prot_metadata),
                     (self.pep_df, 'pep_intensity', pep_metadata),
                     (self.df, 'intensity', trans_metadata)]    # MUST do this last since will overwrite effectively


        # try to get all three data levels from stacked df
        for (level_df, intensity_col, metadata_cols) in df_tuples:

            if intensity_col in self.df.columns:

                # make sure in df
                metadata_cols = [x for x in metadata_cols if x in self.df.columns]

                if intensity_col == 'prot_intensity':
                    sample_data_types = ['prot_intensity']      # no rt or q-val properties at protein level

                elif intensity_col == 'pep_intensity':
                    sample_data_types = \
                        [x for x in ['pep_intensity', 'qval', 'rt_obs', 'rt_pred'] if x in self.df.columns]

                else:
                    sample_data_types = \
                        [x for x in ['intensity', 'qval', 'rt_obs', 'rt_pred'] if x in self.df.columns]

                # this should work regardless of what level the data are
                level_df = pd.pivot_table(self.df, values=sample_data_types, columns=['file_name'],
                                          index=metadata_cols)

                # reset the index to flatten multiple index columns
                level_df.reset_index(inplace=True)

                # python brain fart on why pivot_table assignment above does not achieve code below, but it does not
                if intensity_col == 'prot_intensity':
                    self.prot_df = level_df
                elif intensity_col == 'pep_intensity':
                    self.pep_df = level_df
                else:
                    self.df = level_df


                # ----------------------------------
                # if multiple sample data types are unstacked - ex int and obs_rt, will create a column multi-index

                # first, flatten column index
                level_df.columns = level_df.columns.to_flat_index()

                # then rename tuple header types created by unstacking

                renaming_dict = {}

                name_map = {'qval': self.q_value_suffix, 'intensity': self.intensity_suffix,
                            'rt_pred': self.pred_rt_suffix, 'rt_obs': self.obs_rt_suffix,
                            'pep_intensity': self.intensity_suffix, 'prot_intensity': self.intensity_suffix}

                for col in level_df.columns:
                    if isinstance(col, tuple):
                        if str(col[1]) == '':
                            renaming_dict.update({col: str(col[0])})
                        else:
                            renaming_dict.update({col: (str(col[1]) + name_map[col[0]])})

                if len(renaming_dict) > 0:
                    level_df.rename(columns=renaming_dict, inplace=True)


        # change the stacking status
        self.stacked = False

        self.tracker.capture('unstack')     # capture at this point before next call so no overlap

        # get the sample data columns
        self.get_sample_columns()

        # compare cols after w original to find new sample data columns
        new_cols = [x for x in self.df.columns if x not in before_cols]

        if self.verbose > 0:
            boiler_plates.section_header(title='Stacking summary', level=2)
            print('Columns before unstacking:  ', len(before_cols))
            print('Newly created columns:      ', len(new_cols))
            print('Columns after unstacking:   ', self.df.shape[1])
            print()
            print('Columns removed by unstacking:')
            print()
            removed_cols = [x for x in before_cols if x not in self.df.columns]
            for x in removed_cols:
                print('   ', x)
            print()

        return




    def remove_if_contains(self, column, remove_list, remove_type, **kwargs):
        """
        This method removes all rows if a specified column has any one of a list of values that trigger removal. This
        is used to filter on mods, things flagged for exclusion by Spectronaut, and SWATH decoys.

        This method attempts to do this filtering on all data levels that have existing dataframes (ie all of
        self.df, self.prec_df, self.pep_df, self.prot_df. If the df is None or the column is not in the df, it skips.

        """

        verbose = kwargs.get('verbose', self.verbose)

        before_rows = self.df.shape[0]

        if self.prec_df is not None:
            before_rows_prec = self.prec_df.shape[0]
        else:
            before_rows_prec = np.nan

        if self.pep_df is not None:
            before_rows_pep = self.pep_df.shape[0]
        else:
            before_rows_pep = np.nan

        if self.prot_df is not None:
            before_rows_prot = self.prot_df.shape[0]
        else:
            before_rows_prot = np.nan

        if self.verbose > 0:
            print()
            print('Removing rows with column', column, 'containing:')
            print(remove_list)

        if remove_type is str:

            # brackets in mods have meaning in regular expression, so escape to match to treat literally
            remove_list = [re.escape(m) for m in remove_list]

            self.df = self.df[~self.df[column].str.contains('|'.join(remove_list))]

            if self.prec_df is not None:
                if column in self.prec_df.columns:
                    self.prec_df = self.prec_df[~self.prec_df[column].str.contains('|'.join(remove_list))]

            if self.pep_df is not None:
                if column in self.pep_df.columns:
                    self.pep_df = self.pep_df[~self.pep_df[column].str.contains('|'.join(remove_list))]

            if self.prot_df is not None:
                if column in self.prot_df.columns:
                    self.prot_df = self.prot_df[~self.prot_df[column].str.contains('|'.join(remove_list))]


        elif remove_type is bool:
            if len(remove_list) != 1 or (remove_list[0] is not True and remove_list[0] is not False):
                print('ERROR!   A boolean filter is applied to col:', column, 'expecting only True or False')
                print('         However list passed in is:', remove_list[0])
                quit()

            self.df = self.df[self.df[column] != remove_list[0]]

            if self.prec_df is not None:
                if column in self.prec_df.columns:
                    self.prec_df = self.prec_df[self.prec_df[column] != remove_list[0]]

            if self.pep_df is not None:
                if column in self.pep_df.columns:
                    self.pep_df = self.pep_df[self.pep_df[column] != remove_list[0]]

            if self.prot_df is not None:
                if column in self.prot_df.columns:
                    self.prot_df = self.prot_df[self.prot_df[column] != remove_list[0]]



        if verbose > 0:
            before_after_row_output(before_rows, self.df.shape[0], title='Transition level')

            if self.prec_df is not None:
                before_after_row_output(before_rows_prec, self.prec_df.shape[0], title='Precursor level')
            if self.pep_df is not None:
                before_after_row_output(before_rows_pep, self.pep_df.shape[0], title='Peptide level')
            if self.prot_df is not None:
                before_after_row_output(before_rows_prot, self.prot_df.shape[0], title='Protein level')


        return




    def long_mod_to_tlc(self, cell_string):
        """
        This is a worker function drive by pandas.apply. It converts mods in long form like XXN[Deamidated (N)]XXXXX to
        three letter code (TLC) form like XXN[Dea]XXXXX

        It also converts the terminal mod syntax.
        - Spectronaut denotes mods consistently to the right of the residue: Y[Acetyl (Protein N-term)]XXXXXX
        - The TLC syntax denotes terminal status like:                       [1Ac]-YXXXXXX where "[1Ac]-" is the mod
        - A mod that does not actually require being on the term would be:      Y[TLC]XXXXXX
        - This distinction is more subtle at the C-terminal end:
            ex. real terminal mod:                  AEISFEDRK-[-1K]         mod = "-[-1K]"
            ex. non-terminal mod on the terminus:   GPDVLTATVSGK[2Ox]       mod = "K[2Ox]
        - Some mods have both terminal and residue specificity!
            ex. pyroglutamate  [PGQ]-QLGLPGP[2Ox]PDVPDHAAYHPFR              mod = "[PGQ]-Q"
        - The only way to convert is to look up the mod in the mod_dict to get its specificity.
        - Mods can have multiple specificities. Use terminal if any are terminal or protein terminal

        (placed this method as a method, not external function, so mod_dict visible easily)

        :param cell_string:
        :return:
        """

        if "[" not in cell_string:
            pass
        else:

            if '-term' in cell_string:
                term_present = True
            else:
                term_present = False


            # keep this development output in code at high verbose level since mods are always problematic!
            if self.verbose > 2 and term_present is True:
                print()
                print('INITIAL:', cell_string)


            parts_list = cell_string.replace(']', '&[').split('[')

            # N-term mods - check if the second part is a mod and determine what kind
            if "&" in parts_list[1] and 'N-term' in parts_list[1]:

                if self.verbose > 2:
                    print('BEFORE:', parts_list)

                parts_list.insert(0, parts_list.pop(1))
                parts_list.insert(1, '-')

                if self.verbose > 2:
                    print('AFTER :', parts_list)

            # C-term mods - check if the last part is a mod and determine what kind
            if "&" in parts_list[-1] and 'C-term' in parts_list[1]:

                if self.verbose > 2:
                    print('BEFORE:', parts_list)

                # insert a hyphen part before the last part to indicate C-terminal in TLC system
                if parts_list[-2] != '-':   # this should never be the case in Spectronaut, but being paranoid here
                    parts_list.insert(-1, '-')

                if self.verbose > 2:
                    print('AFTER :', parts_list)


            new_list = [('[' + self.mod_dict[x.replace('&', '')] + ']') if '&' in x else x for x in parts_list]

            cell_string = ''.join(new_list)

            if self.verbose > 2:
                if term_present is True:
                    print('FINAL:', cell_string)

        return cell_string




    def analyze_dilution_series(self, row):

        ldr = np.nan
        ave_acc = np.nan
        ave_rms = np.nan

        # skip if not the variable component or fraction complete very low
        if row['spike_status'] == 'var' and row['f_complete'] > 0.1:

            ints = row[self.intensity_columns].tolist()

            try:
                ldr, ave_acc, ave_rms = \
                    quant.analyze_dilution_series(self.sample_dilutions, ints,
                                                  method='linear from top',
                                                  max_percent_rms_ldr=self.max_percent_rms_ldr,
                                                  acc_rms_fixed_dil=self.acc_rms_fixed_dil,
                                                  show_plots=False, verbose=self.verbose)

                self.ldr_successes += 1

            except:
                self.ldr_errors += 1

        else:
            self.ldr_skipped += 1

        return ldr, ave_acc, ave_rms




    def characterize(self, **kwargs):

        verbose = kwargs.get('verbose', self.verbose)
        full_analysis = kwargs.get('full_analysis', False)
        title = kwargs.get('title', '')
        show_plots = kwargs.get('show_plots', False)
        plots_to_file = kwargs.get('plots_to_file', False)
        tables_to_file = kwargs.get('tables_to_file', False)
        outpath = kwargs.get('outpath', os.path.join(self.outpath, 'CHARACTERIZE'))


        if verbose > 0:
            title = 'Characterize - Stage: ' + str(title)
            boiler_plates.section_header(title=title, level=1)


        # make sure outpath exists if needed - default is a new sub-folder in class-level output folder
        if outpath is not None and (plots_to_file is True or tables_to_file is True) \
                and not os.path.exists(outpath):
            os.makedirs(outpath)


        # ---------------------------------
        # Determine dimensions of transition table

        num_proteins = self.df['prot_ID'].nunique()  # always have this
        num_peptides = self.df['pep_ID'].nunique()
        num_prec = self.df['prec_ID'].nunique()
        num_trans = self.df['trans_ID'].nunique()
        matrix_size = num_trans * len(self.intensity_columns)


        self.char_dict.update({'level': self.level,

                               'intensity_columns': self.intensity_columns,
                               'num_samples': len(self.intensity_columns),

                               'num_proteins': num_proteins,
                               'num_peptides': num_peptides,
                               'num_prec': num_prec,
                               'num_trans': num_trans,
                               })


        # ---------------------------------
        # Stats on three levels if possible

        for (df, level) in [(self.df, 'trans'), (self.pep_df, 'pep'), (self.prot_df, 'prot')]:

            if df is not None and self.lib_or_results == 'results':

                # ---------------------------------
                # Basics

                self.char_dict.update({level + '_max_intensity': df[self.intensity_columns].max().max(),
                                       level + '_shape': df.shape
                                       })


                # ---------------------------------
                # Row-wise completeness

                df['counts'] = df[self.intensity_columns].notnull().sum(axis=1)
                df['f_complete'] = df['counts'] / len(self.intensity_columns)


                # ---------------------------------
                # Aggregate completeness

                matrix_size = df.shape[0] * len(self.intensity_columns)    # max possible values
                total_measurements = df['counts'].sum()                     # non-missing values (not nan)

                self.char_dict.update({level + '_matrix_size': matrix_size,
                                       level + '_total_measurements': total_measurements,
                                       level + '_missing_values': matrix_size - total_measurements,
                                       level + '_percent_complete': 100 * total_measurements / matrix_size})


                # ---------------------------------
                # Ave int, stdev, and CVs

                df['ave_int'] = df[self.intensity_columns].mean(axis=1, skipna=True)
                df['ave_log_int'] = df[self.intensity_columns].apply(np.log10).mean(axis=1, skipna=True)
                df['stdev_int'] = df[self.intensity_columns].std(axis=1, skipna=True)
                df['cv_int'] = df['stdev_int'] / df['ave_int']

                df['max_int'] = df[self.intensity_columns].max(axis=1, skipna=True)
                df['max_log_int'] = np.log10(df['max_int'])

                # if dilution expt, mask out data for variable component (not labeled const - also mixed and unknown)
                if self.sample_dilutions is not None and self.expt_type == 'dilution_series_species_spike':

                    df['ave_int'].mask((df['spike_status'] != 'const'), inplace=True)
                    df['ave_log_int'].mask((df['spike_status'] != 'const'), inplace=True)
                    df['stdev_int'].mask((df['spike_status'] != 'const'), inplace=True)
                    df['cv_int'].mask((df['spike_status'] != 'const'), inplace=True)


                # ---------------------------------
                # Intensity percentiles

                df['int_percentile'] = 100 * df['ave_int'].rank(pct=True)


                # ---------------------------------
                # Metrics from dilution experiments - ie w know relative conc for some analytes - ex LDR, acc, rms

                if full_analysis is True and self.sample_dilutions is not None:

                    if self.verbose > 0:
                        print()
                        print('LDR and accuracy computation for dilution experiments')
                        print('Level table shape:', df.shape)

                    # reset counter for this level
                    self.ldr_errors = 0
                    self.ldr_successes = 0
                    self.ldr_skipped = 0

                    df[['ldr', 'ave_acc', 'ave_rms']] = \
                        df[['spike_status', 'f_complete'] + self.intensity_columns].apply(
                            lambda row: self.analyze_dilution_series(row),
                            axis=1, result_type='expand')

                    # force types and replace non-numericals (polishes annoying bug I can't find)
                    df[['ldr', 'ave_acc', 'ave_rms']] = \
                        df[['ldr', 'ave_acc', 'ave_rms']].astype('float64') #.replace(0, np.nan)

                    df['log_ldr'] = np.log10(df['ldr'])

                    if self.verbose > 0:
                        print('Level:', level, '  Errors:', self.ldr_errors, ' Successes:', self.ldr_successes,
                              'Skipped:', self.ldr_skipped)


                # ---------------------------------
                # Compute scores metrics - Need these scores at both trans and precursor levels for filtering purposes

                for (metric_cols, suffix) in [(self.score_columns, '_score'), (self.q_value_columns, '_qval'),
                                              (self.rt_columns, '_obs_rt'), (self.delta_rt_columns, '_delta_rt')]:

                    # make sure cols are in df - handles absence of data at some levels
                    cols_in_df = sum([1 for x in metric_cols if x in df.columns])

                    if cols_in_df == len(metric_cols):
                        df[('ave' + suffix)] = df[metric_cols].mean(axis=1, skipna=True)
                        df[('stdev' + suffix)] = df[metric_cols].std(axis=1, skipna=True)
                        df[('max' + suffix)] = df[metric_cols].max(axis=1, skipna=True)
                        df[('min' + suffix)] = df[metric_cols].min(axis=1, skipna=True)



        # ---------------------------------
        # Spike status - for spike in experiments - is the analyte constant, varied, mixed or unknown

        if 'spike_status' in self.df.columns:

            spike_status_counts_trans = self.df['spike_status'].value_counts(dropna=False)


        # ---------------------------------
        # Level counts histograms

        hist_trans_per_prec, hist_pep_per_prot = self.get_analytes_per_analyte()


        # ---------------------------------
        # console output

        if verbose > 0:

            boiler_plates.section_header(title='Dimensions of Transitions Table', level=2)
            print('Samples in table:', len(self.intensity_columns))
            print()
            print('Proteins         ', format(num_proteins, ','))
            print('Peptides         ', format(num_peptides, ','))
            print('Precursors:      ', format(num_prec, ','))
            print('Transitions      ', format(num_trans, ','))

            self.show_analyte_count_progression()
            print()


            for (stem, name) in [('trans', 'Transition'), ('pep', 'Peptide'), ('prot', 'Protein')]:

                if (stem + '_total_measurements') in self.char_dict:    # check data level present by any level key

                    # Dimensions
                    boiler_plates.section_header(title=('Table Dimensions:' + name), level=2)
                    print('Shape:                ', self.char_dict[(stem + '_shape')])
                    print('Max intensity         ', format(round(self.char_dict[(stem + '_max_intensity')], 0), ','))
                    print()

                    # Data completeness
                    boiler_plates.section_header(title=('Data Completeness:' + name), level=2)
                    print('Total measurements:   ', format(self.char_dict[(stem + '_total_measurements')], ','))
                    print('Missing values:       ', format(self.char_dict[(stem + '_missing_values')], ','))
                    print('Data completeness:    ', round(self.char_dict[(stem + '_percent_complete')], 2), '%')
                    print()


            # histogram tables - trans per prec and pep per prot
            if self.level == 'transition':
                boiler_plates.section_header(title='Frequency of Transitions per Precursor', level=2)
                print(tabulate(hist_trans_per_prec, headers='keys', tablefmt='simple'))
                print()
                print('    ', hist_trans_per_prec['prec'].sum(), 'total precursors')
                print()

            if self.level != 'protein':
                boiler_plates.section_header(title='Frequency of Peptides per Protein (truncated to 30)', level=2)
                print(tabulate(hist_pep_per_prot.iloc[:31], headers='keys', tablefmt='simple'))
                print()
                print('max peptides/protein: ', hist_pep_per_prot.index.values.max())
                print()

            # spike status breakdown
            if 'spike_status' in self.df.columns:
                boiler_plates.section_header(title='Frequency of Transition Spike Statuses', level=2)
                print(spike_status_counts_trans)
                print()


            # -----------------------------------------------------
            # Characterize mods (must be TLC mod system) - done in console output as function provides output

            self.characterize_mods()

        self.tracker.capture('characterize')

        # -----------------------------------------------------
        # Generate histogram analysis and plots

        if full_analysis is True and self.lib_or_results == 'results':
            self.hist_and_plots(show_plots=show_plots,
                                plots_to_file=plots_to_file, tables_to_file=tables_to_file,
                                outpath=outpath)

        return




    def hist_and_plots(self, **kwargs):

        verbose = kwargs.get('verbose', self.verbose)
        show_plots = kwargs.get('show_plots', False)
        plots_to_file = kwargs.get('plots_to_file', False)
        tables_to_file = kwargs.get('tables_to_file', False)
        height_inches = kwargs.get('height_inches', self.plot_height_inches)
        width_inches = kwargs.get('width_inches', self.plot_width_inches)
        outpath = kwargs.get('outpath', os.path.join(self.outpath, 'CHARACTERIZE'))

        if verbose > 0:
            boiler_plates.section_header(title='Get histograms and plots', level=2)

        # ---------------------------------
        # Stats on three levels if possible

        for (df, level, name) in [(self.df, 'trans', 'Transitions'),
                                  (self.pep_df, 'pep', 'Peptides'),
                                  (self.prot_df, 'prot', 'Proteins')]:

            if df is not None:

                plotting.histogram_continuous(df, 'counts',
                                              hist=True,
                                              norm_hist=False,
                                              norm_hist_as_prob=False,
                                              bins=plotting.generate_bin_series(-0.5,
                                                                                len(self.intensity_columns) + 0.5,
                                                                                1),
                                              kde=False,
                                              rug=False,
                                              title=('Analyte Data Completeness Distribution - ' + name),
                                              ylabel=name,
                                              xlabel='measurements',
                                              xlim=None,
                                              grid_maj_x=False,
                                              grid_maj_y=False,
                                              show_plots=show_plots,
                                              height_inches=height_inches,
                                              width_inches=width_inches,
                                              plots_to_file=plots_to_file,
                                              tables_to_file=tables_to_file,
                                              # metadata_files=metadata_files,
                                              outpath=outpath, file_stem=('Completeness Distribution - ' + name),
                                              verbose=verbose)

                # plot data completeness vs. average intensity of the analyte

                # get moving stats to overlay on plot
                moving_stats_df = plotting.scatter_sliding_stats(df, x='ave_log_int', y='f_complete',
                                                                 analyses='scatter+percentiles',
                                                                 verbose=verbose)

                plotting.pandas_plot(df, x='ave_log_int', y='f_complete', kind='scatter',
                                     ylim=(0, 1.0),
                                     xlim=None,
                                     title=('Data completeness vs. Average Intensity - ' + name),
                                     xlabel='average log10(intensity)',
                                     ylabel='completeness fraction',
                                     grid_maj_x=False,
                                     grid_maj_y=False,
                                     moving_stats_df=moving_stats_df,
                                     analysis='scatter+percentiles',
                                     height_inches=height_inches,
                                     width_inches=width_inches,
                                     show_plots=show_plots,
                                     plots_to_file=plots_to_file,
                                     tables_to_file=tables_to_file,
                                     # metadata_files=metadata_files,
                                     outpath=outpath, file_stem=('Completeness v Ave Intensity - ' + name),
                                     verbose=verbose
                                     )

                # compute histogram of the quant values for the average log intensity

                max_intensity = self.char_dict[(level + '_max_intensity')]

                plotting.histogram_continuous(df, 'ave_log_int',
                                              hist=True,
                                              norm_hist=True,
                                              norm_hist_as_prob=False,
                                              bins=plotting.generate_bin_series(0, np.log10(max_intensity),
                                                                                np.log10(max_intensity) / 50),
                                              kde=True,
                                              rug=False,
                                              title=('Intensity Distribution - ' + name),
                                              ylabel='frequency',
                                              xlabel='average log10(intensity)',
                                              xlim=None,
                                              grid_maj_x=False,
                                              grid_maj_y=False,
                                              height_inches=height_inches,
                                              width_inches=width_inches,
                                              show_plots=show_plots,
                                              plots_to_file=plots_to_file, tables_to_file=tables_to_file,
                                              # metadata_files=metadata_files,
                                              outpath=outpath, file_stem=('Intensity Distribution - ' + name),
                                              verbose=verbose)


                # -----------------------------------------------
                # plot cv vs. log average intensity of the analyte

                # get moving stats to overlay on plot
                moving_stats_df = plotting.scatter_sliding_stats(df, x='ave_log_int', y='cv_int',
                                                                 analyses='scatter+percentiles',
                                                                 verbose=verbose)


                plotting.pandas_plot(df, x='ave_log_int', y='cv_int', kind='scatter',
                                     title=('cv vs. Average Intensity - ' + name),
                                     xlabel='log10(ave intensity)',
                                     ylabel='cv',
                                     xlim=None,
                                     ylim=(0, 2),
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
                                     outpath=outpath, file_stem=('cv vs. Average Intensity - ' + name),
                                     verbose=verbose
                                     )



                # -------------------------------------------------
                # compute histogram of the linear dynamic range

                if 'log_ldr' in df.columns:
                    
                    plotting.histogram_continuous(df, 'log_ldr',
                                                  hist=True,
                                                  norm_hist=True,
                                                  norm_hist_as_prob=False,
                                                  bins=plotting.generate_bin_series(0, 2, 0.02),
                                                  kde=True,
                                                  rug=False,
                                                  title=('LDR Distribution - ' + name),
                                                  ylabel='frequency',
                                                  xlabel='log10(ldr)',
                                                  xlim=None,
                                                  grid_maj_x=False,
                                                  grid_maj_y=False,
                                                  height_inches=height_inches,
                                                  width_inches=width_inches,
                                                  show_plots=show_plots,
                                                  plots_to_file=plots_to_file, tables_to_file=tables_to_file,
                                                  # metadata_files=metadata_files,
                                                  outpath=outpath, file_stem=('LDR Distribution - ' + name),
                                                  verbose=verbose)


                    # -----------------------------------------------
                    # plot ldr vs. log max intensity of the analyte

                    # get moving stats to overlay on plot
                    moving_stats_df = plotting.scatter_sliding_stats(df, x='max_log_int', y='log_ldr',
                                                                     analyses='scatter+percentiles',
                                                                     verbose=verbose)


                    plotting.pandas_plot(df, x='max_log_int', y='log_ldr', kind='scatter',
                                         title=('LDR vs. Max Intensity - ' + name),
                                         xlabel='log10(max intensity)',
                                         ylabel='log10(linear dynamic range)',
                                         xlim=None,
                                         ylim=(0, 2),
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
                                         outpath=outpath, file_stem=('LDR vs. Max Intensity - ' + name),
                                         verbose=verbose
                                         )



                # -------------------------------------------------
                # plots only at peptide level

                if level == 'pep':

                    # -----------------------------------------------
                    # delta rt histogram analysis

                    plotting.histogram_continuous(df, 'ave_delta_rt',
                                                  hist=True,
                                                  norm_hist=True,
                                                  norm_hist_as_prob=False,
                                                  bins=plotting.generate_bin_series(0, 2, (2 / 100)),
                                                  kde=True,
                                                  rug=False,
                                                  title=('Average Delta RT Distribution - ' + name),
                                                  ylabel='frequency',
                                                  xlabel='abs(obs rt - pred rt)',
                                                  xlim=None,
                                                  grid_maj_x=False,
                                                  grid_maj_y=False,
                                                  height_inches=height_inches,
                                                  width_inches=width_inches,
                                                  show_plots=show_plots,
                                                  plots_to_file=plots_to_file, tables_to_file=tables_to_file,
                                                  # metadata_files=metadata_files,
                                                  outpath=outpath, file_stem=('Delta RT Distribution - ' + name),
                                                  verbose=verbose)


                    # -----------------------------------------------
                    # delta rt vs. log average peptide intensity

                    # get moving stats to overlay on plot
                    moving_stats_df = plotting.scatter_sliding_stats(df, x='ave_log_int', y='ave_delta_rt',
                                                                     analyses='scatter+percentiles',
                                                                     verbose=verbose)


                    plotting.pandas_plot(df, x='ave_log_int', y='ave_delta_rt', kind='scatter',
                                         title=('Ave Delta RT vs. Average Intensity - ' + name),
                                         xlabel='log10(ave intensity)',
                                         ylabel='abs(obs_rt - pred_rt)',
                                         xlim=None,
                                         ylim=(0, 2),
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
                                         outpath=outpath, file_stem=('Delta RT vs. Average Intensity - ' + name),
                                         verbose=verbose
                                         )


        # -------------------------------------------------
        # plots only at transition level

        if 'ave_score' in self.intensity_columns:

            moving_stats_df = plotting.scatter_sliding_stats(self.df, x='ave_log_int', y='ave_score',
                                                             analyses='scatter+percentiles',
                                                             verbose=verbose)


            plotting.pandas_plot(self.df, x='ave_log_int', y='ave_score', kind='scatter',
                                 title='cv vs. Average Intensity',
                                 xlabel='log10(ave intensity)',
                                 ylabel='ave_score',
                                 xlim=None,
                                 ylim=None,
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
                                 outpath=outpath, file_stem='ave_score vs. Average Intensity',
                                 verbose=verbose
                                 )

        if 'max_score' in self.intensity_columns:

            plotting.pandas_plot(self.df, x='max_score', y='cv_int', kind='scatter',
                                 title='cv vs. Max Score',
                                 xlabel='max score',
                                 ylabel='cv',
                                 xlim=None,
                                 ylim=None,
                                 grid_maj_x=True,
                                 grid_maj_y=True,
                                 #moving_stats_df=moving_stats_df,
                                 #analysis='scatter+percentiles',
                                 height_inches=height_inches,
                                 width_inches=width_inches,
                                 show_plots=show_plots,
                                 plots_to_file=plots_to_file,
                                 tables_to_file=tables_to_file,
                                 # metadata_files=metadata_files,
                                 outpath=outpath, file_stem='cv vs. max_score',
                                 verbose=verbose
                                 )

        self.tracker.capture('hist and plots')

        return




    def get_analytes_per_analyte(self):
        """
        Get columns and histograms for transitions per precursor and peptides per protein.

        """

        # distribution of transitions per precursor

        if 'prec_ID' in self.df.columns:
            # get the counts of each unique precursor
            trans_counts_for_each_prec = self.df['prec_ID'].value_counts()

            # write these back as a new column that we can filter on
            self.df['trans_per_prec'] = trans_counts_for_each_prec.loc[self.df['prec_ID']].values

            hist_trans_per_prec = pd.DataFrame(
                {'prec': trans_counts_for_each_prec.value_counts(dropna=True, normalize=False)})

            # reindex to span full range with no gaps in order (round avoids float to int conversion issues)
            complete_count_series = range(0, (int(round(hist_trans_per_prec.index.values.max(), 0)) + 1), 1)
            hist_trans_per_prec = hist_trans_per_prec.reindex(complete_count_series, fill_value=0)

            # add normalized column
            hist_trans_per_prec['fract_prec'] = hist_trans_per_prec['prec'] / hist_trans_per_prec['prec'].sum()


        # distribution of peptides per protein (don't really care about precursors per peptide distrib)
        # this is slightly more complicated than the above hist since can have multiple prec per pep

        if 'pep_ID' in self.df.columns:
            # make a temp df and reduce it to unique peptides
            temp_df = self.df[['prot_ID', 'pep_ID']].copy(deep=True)
            temp_df.drop_duplicates(subset=['pep_ID'], inplace=True)

            # get the counts of each protein
            pep_counts_for_each_prot = temp_df['prot_ID'].value_counts()

            # write these back as a new column that we can filter on
            self.df['pep_per_prot'] = pep_counts_for_each_prot.loc[self.df['prot_ID']].values

            hist_pep_per_prot = pd.DataFrame(
                {'prot': pep_counts_for_each_prot.value_counts(dropna=True, normalize=False)})

            # reindex to span full range with no gaps in order (round avoids float to int conversion issues)
            complete_count_series = range(0, (int(round(hist_pep_per_prot.index.values.max(), 0)) + 1), 1)
            hist_pep_per_prot = hist_pep_per_prot.reindex(complete_count_series, fill_value=0)

            # add normalized column
            hist_pep_per_prot['fract_prot'] = hist_pep_per_prot['prot'] / hist_pep_per_prot['prot'].sum()


        return hist_trans_per_prec, hist_pep_per_prot




    def filter_by_completeness(self, **kwargs):

        min_count = kwargs.get('min_count', 1)
        min_fract_complete = kwargs.get('min_fract_complete', None)


        if self.lib_or_results != 'results':
            return

        if self.verbose > 0:
            boiler_plates.section_header(title='Filter by Completeness', level=1)

        for (df, level) in [(self.df, 'trans'), (self.prec_df, 'pep'), (self.pep_df, 'pep'), (self.prot_df, 'prot')]:

            if df is not None:

                initial_rows = df.shape[0]

                # recompute counts to be sure current - faster than calling characterize()
                df['counts'] = df[self.intensity_columns].notnull().sum(axis=1)
                df['f_complete'] = df['counts'] / len(self.intensity_columns)

                df.drop(df[df['counts'] < min_count].index, inplace=True)

                if min_fract_complete is not None:
                    df.drop(df[df['f_complete'] < min_fract_complete].index, inplace=True)

                if self.verbose > 0:
                    print(level, ':  Initial rows:', initial_rows, ':  Final rows:', df.shape[0], 'min:', min_count)

        self.count_analytes('filter by completeness')

        self.tracker.capture('filter by completeness')

        return




    def filter_by_scores(self, **kwargs):
        """
        This method removes...

        :param kwargs:
        :return:
        """

        max_qvalue = kwargs.get('max_qvalue', 0.01)
        min_score = kwargs.get('min_score', None)

        if self.lib_or_results == 'library':
            return

        if self.verbose > 0:
            boiler_plates.section_header(title='Filter by Scores', level=1)


        def process_level(df, level):

            # --------------------------------------
            # Mask peak groups below the min score

            # check if score columns in this data level
            score_col_count = sum([1 for x in self.score_columns if x in df.columns])

            if score_col_count > 0 and min_score is not None:

                non_empty_values_before = df[self.intensity_columns].count().sum()

                for name_stem in self.sample_name_stems:

                    area_col = name_stem + self.intensity_suffix
                    score_col = name_stem + self.score_suffix

                    df[area_col].mask((df[score_col] < min_score), inplace=True)

                non_empty_values_after = df[self.intensity_columns].count().sum()

                if self.verbose > 0:
                    title = 'Filter by score at level:' + level
                    boiler_plates.section_header(title=title, level=2)
                    print('Min score allowed:        ', min_score)
                    print('Values before filter:     ', format(non_empty_values_before, ','))
                    print('Values after filter:      ', format(non_empty_values_after, ','))
                    print('Removed percentage:       ',
                          round(100 * (non_empty_values_before - non_empty_values_after) / non_empty_values_before, 2))


            # --------------------------------------
            # Mask peak groups above the max q-value

            # check if q-val columns in this data level
            qval_col_count = sum([1 for x in self.q_value_columns if x in df.columns])

            if qval_col_count > 0 and max_qvalue is not None:

                non_empty_values_before = df[self.intensity_columns].count().sum()

                for name_stem in self.sample_name_stems:

                    area_col = name_stem + self.intensity_suffix
                    qval_col = name_stem + self.q_value_suffix

                    df[area_col].mask(df[qval_col].isna(), inplace=True)        # nans not covered by line below
                    df[area_col].mask(df[qval_col] > max_qvalue, inplace=True)

                non_empty_values_after = df[self.intensity_columns].count().sum()

                if self.verbose > 0:
                    title = 'Filter by q-value at level:' + level
                    boiler_plates.section_header(title=title, level=2)
                    print()
                    print('Max q-value allowed:      ', max_qvalue)
                    print('Values before filter:     ', format(non_empty_values_before, ','))
                    print('Values after filter:      ', format(non_empty_values_after, ','))
                    print('Removed percentage:       ',
                          round(100 * (non_empty_values_before - non_empty_values_after) / non_empty_values_before, 2))


            # --------------------------------------
            # Mask other data arrays after all intensity data masking done

            # todo - could make this optional in case you want to look at data outside bounds, but prob better done via threshold change

            for (data_type, column_list, suffix) in self.data_arrays_list[1:]:

                if len(column_list) > 0:

                    count = 0

                    for int_col in self.intensity_columns:

                        corresponding_col = int_col.replace(self.intensity_suffix, suffix)

                        if corresponding_col in df.columns:
                            df[corresponding_col].mask(df[int_col].isna(), inplace=True)
                            count += 1

                    # count because you won't have some cols at higher levels - ex no RT or q-val for protein level
                    if count > 0:
                        print('Values after filter:      ', format(df[column_list].count().sum(), ','), 'in array:', data_type)

            return df


        # --------------------------------------
        # process each level if present

        if self.df is not None:

            self.df = process_level(self.df, 'Transition')

        if self.pep_df is not None:
            self.pep_df = process_level(self.pep_df, 'Peptide')

        if self.prot_df is not None:
            self.prot_df = process_level(self.prot_df, 'Protein')


        # get counts (shouldn't change based on current logic removing only single points, not whole rows)
        self.count_analytes('filter by q-values and scores')


        self.tracker.capture('filter by scores')

        return




    def sample_normalization(self, **kwargs):
        """
        Determines sample normalization via various approaches - all at the transition level.

        IMPORTANT: Relation with provided peptide and protein rollup: if ignore_provided_rollup is False, as is the
        default, it will use the provided rollup but apply the normalization determined at the transition level for each
        sample in this method, if called. If ignore_provided_rollup is True, then we're going to recompute the rollup
        later (after the normalization), in which case, the higher levels will automatically acquire the normalization
        results without any necessary action in this method.

        Logic on normalization already present in data controlled by ignore_provided_norm:
        If the sample norm method is called, it's assumed you want to normalized the data. However, if you also pass
        ignore_provided_norm = False, then if the input data was normalized or provided a normalized form (OneOmics
        provides both norm'd and un-norm'd selected upon instantiation of SwathData class), normalization will NOT
        be triggered here if we can both tell the input data was normalized AND you prefer to use that norm if
        present by passing ignore_provided_norm = False.

        If you want it to normalize no matter what calling this method, then no need to pass anything as
        ignore_provided_norm is assumed True by default if you call this method. However, upon import to SwathData
        class, it is NOT assumed True.

        :param kwargs:
        min_percentile: exclude data from below this percentile in determining the normalization factor
        max_percentile: exclude data from above this percentile in determining the normalization factor
        :return:
        """

        ignore_provided_norm = kwargs.get('ignore_provided_norm', True)
        method = kwargs.get('method', 'ratio_to_analyte_ave')   # choices: 'ratio_to_analyte_ave', None,
        apex_method = kwargs.get('apex_method', 'mode_by_ave_of_top')   # 'mode'
        min_percentile = float(kwargs.get('min_percentile', 50))
        max_percentile = float(kwargs.get('max_percentile', 100))
        min_f_complete = kwargs.get('min_f_complete', 0.80)
        min_score = kwargs.get('min_score', 0.8)
        max_qval = kwargs.get('max_qval', 0.01)
        ignore_provided_rollup = kwargs.get('ignore_provided_rollup', False)

        show_plots = kwargs.get('show_plots', False)
        plots_to_file = kwargs.get('plots_to_file', False)
        tables_to_file = kwargs.get('tables_to_file', False)
        height_inches = kwargs.get('height_inches', self.plot_height_inches)
        width_inches = kwargs.get('width_inches', self.plot_width_inches)
        outpath = kwargs.get('outpath', os.path.join(self.outpath, 'CHARACTERIZE'))


        if self.lib_or_results == 'library' or\
                method is None or (ignore_provided_norm is False and self.provided_norm_present is True):
            return


        # make sure outpath exists if needed - default is a new sub-folder in class-level output folder
        if outpath is not None and (plots_to_file is not None or tables_to_file is not None) \
                and not os.path.exists(outpath):
            os.makedirs(outpath)


        # recompute this always (in light mode) to be sure current
        # todo - refactor necessary bits from this function so don't need to call whole thing!
        self.characterize(verbose=0, full_analysis=False, outpath=None)    # not passing any output generating kwargs


        if self.verbose > 0:
            boiler_plates.section_header(title='Sample Normalization', level=1)
            print('Method:               ', method)
            print()
            print('Parameters for selection of analytes to use for normalization:')
            print('   min percentile:    ', min_percentile)
            print('   max percentile:    ', max_percentile)
            print('   min f complete:    ', min_f_complete)
            print('   max q-value:       ', max_qval)
            print('   min score:         ', min_score)
            print()


        if method == 'ratio_to_analyte_ave':

            col_list = ['ave_int', 'cv_int', 'stdev_int',
                        'int_percentile',
                        'f_complete',
                        ]

            if 'spike_status' in self.df.columns:
                col_list.append('spike_status')

            if len(self.score_columns) > 0:
                col_list += ['min_score', 'max_score', 'ave_score', 'stdev_score']

            if len(self.q_value_columns) > 0:
                col_list += ['min_qval', 'max_qval', 'ave_qval', 'stdev_qval']

            norm_df = self.df[col_list].copy(deep=True)


            # ---------------------------------------------
            # Determine which analytes contribute to normalization by inclusion in ref_series

            # create initial ref series
            values_initial = norm_df['ave_int'].count()
            norm_df['ref_series'] = norm_df['ave_int']


            # filter to constant species if mixed species spike type expt
            if 'spike_status' in norm_df.columns and self.expt_type == 'dilution_series_species_spike':

                norm_df['ref_series'].mask((norm_df['spike_status'] != 'const'), inplace=True)

            values_after_species_requirements = norm_df['ref_series'].count()


            # remove values below min or above max percentile, setting to nan
            norm_df['ref_series'].mask(
                (norm_df['int_percentile'] < min_percentile) | (norm_df['int_percentile'] > max_percentile),
                inplace=True)

            values_after_intensity_requirements = norm_df['ref_series'].count()


            # further require a high detection rate (high data completeness), which biases toward lower variability
            norm_df['ref_series'].mask((norm_df['f_complete'] < min_f_complete), inplace=True)

            values_after_high_completeness_requirement = norm_df['ref_series'].count()


            # further require all high quality detections
            if 'max_qval' in norm_df.columns:
                norm_df['ref_series'].mask((norm_df['max_qval'] > max_qval), inplace=True)

            values_after_low_qval_requirement = norm_df['ref_series'].count()

            """
            # further require all high quality detections
            if 'min_score' in norm_df.columns:
                norm_df['ref_series'].mask((norm_df['min_score'] < min_score), inplace=True)

            values_after_high_score_requirement = norm_df['ref_series'].count()
            """

            # further explicitly select the lower variability half of analytes still remaining
            median_cv_current_set = np.nanmedian(norm_df['stdev_int'] / norm_df['ref_series'])

            norm_df['ref_series'].mask((norm_df['cv_int'] > median_cv_current_set), inplace=True)

            values_after_low_variability_requirement = norm_df['ref_series'].count()

            if self.verbose > 0:
                print()
                print('Selection for normalization calculation:')
                print()
                print('   Initial rows:                             ', values_initial)
                print('   After species requirements:               ', values_after_species_requirements)
                print('   After intensity quantile requirements:    ', values_after_intensity_requirements)
                print('   After high data completeness requirement: ', values_after_high_completeness_requirement)
                print('   After low max q-value req for trans row:  ', values_after_low_qval_requirement)
                #print('   After high min score req for trans row:   ', values_after_high_score_requirement)
                print('   After low variation requirement:          ', values_after_low_variability_requirement)
                print('         (cv less than median cv:', median_cv_current_set, ')')
                print()


            hists_df = pd.DataFrame()

            # go through all columns
            for sample in self.intensity_columns:

                norm_df[sample] = np.log10(self.df[sample] / norm_df['ref_series'])

                # compute histogram of the quant values for the first replicate

                table_df = \
                    plotting.histogram_continuous(norm_df, sample,
                                                  hist=True,
                                                  norm_hist=True,
                                                  norm_hist_as_prob=False,
                                                  bins=plotting.generate_bin_series(-2.01, 2.01, 0.02),
                                                  kde=True,
                                                  rug=False,
                                                  title=('Pre-Norm Selected Set Distrib: \n' + str(sample)),
                                                  ylabel='frequency',
                                                  xlabel='log10(sample / sample ave)',
                                                  xlim=None,
                                                  ylim=None,
                                                  grid_maj_x=True,
                                                  grid_maj_y=False,
                                                  height_inches=height_inches,
                                                  width_inches=width_inches,
                                                  show_plots=show_plots,
                                                  plots_to_file=plots_to_file, tables_to_file=tables_to_file,
                                                  # metadata_files=metadata_files,
                                                  outpath=outpath, file_stem=('Norm ratios - ' + str(sample)),
                                                  verbose=self.verbose)

                # get norm factor
                if apex_method == 'mode_by_ave_of_top':
                    hist_max = table_df['density'].max()
                    log_factor = table_df.loc[table_df['density'] > hist_max * 0.95, 'bin_center'].mean()

                elif apex_method == 'mode':
                    hist_max = table_df['density'].max()
                    log_factor = table_df.loc[table_df['density'] == hist_max, 'bin_center'].mean()

                else:
                    log_factor = None
                    print('FATAL ERROR! Apex method unrecognized')
                    quit()

                # convert out of log space to actual factor
                factor = 10**log_factor

                # update the factor - assume existing value present - 1.0 if not norm'd before or prior correction
                self.sample_norm_factors.update({sample: factor * self.sample_norm_factors[sample]})

                table_df.set_index('bin_center', inplace=True, drop=True)
                hists_df[sample] = table_df['density']


            # todo - need to create df aligned on same x series - not ensured w kde currently
            #print('Histogram data')
            #print(hists_df)


            # get the kde histogram traces for all samples to char dict
            self.char_dict.update({'pre_norm_distributions': hists_df})


        # apply normalization to transition data
        self.df[self.intensity_columns] = self.df[self.intensity_columns] / pd.Series(self.sample_norm_factors)

        # apply same normalization factors to peptide table if appropriate
        if ignore_provided_rollup is False and self.pep_df is not None:
            self.pep_df[self.intensity_columns] = self.pep_df[self.intensity_columns] / pd.Series(self.sample_norm_factors)

        # apply same normalization factors to protein table if appropriate
        if ignore_provided_rollup is False and self.prot_df is not None:
            self.prot_df[self.intensity_columns] = self.prot_df[self.intensity_columns] / pd.Series(self.sample_norm_factors)


        # descriptive stats on factors
        factor_list = [self.sample_norm_factors[x] for x in self.sample_norm_factors]
        ave_norm_factor = float(np.mean(factor_list))
        stdev_norm_factor = float(np.std(factor_list))
        cv_norm_factor = stdev_norm_factor / ave_norm_factor

        if self.verbose > 0:
            print()
            print('Determined normalization factors:')
            print()
            for sample in self.intensity_columns:
                print('factor:', round(self.sample_norm_factors[sample], 3), '   sample:', sample)
            print()
            print('average factor:    ', round(ave_norm_factor, 3))
            print('stdev factor:      ', round(stdev_norm_factor, 3))
            print('cv factor:         ', round(cv_norm_factor, 3))
            print()

        self.normalized_status = True

        self.tracker.capture('sample normalization')

        return




    def count_analytes(self, state):

        trans = self.df['trans_ID'].nunique()
        prec = self.df['prec_ID'].nunique()
        pep = self.df['pep_ID'].nunique()
        prot = self.df['prot_ID'].nunique()

        self.analyte_tracker.append([state, trans, prec, pep, prot])

        return




    def show_analyte_count_progression(self):

        if self.verbose == 0:
            return

        headers = ['State',
                   'Transitions', 'Precursors', 'Peptides', 'Proteins',
                   'Delta trans', 'Delta prec', 'Delta pep', 'Delta prot']

        first_row = self.analyte_tracker[0] + ['', '', '', '']
        datatable = [first_row]

        for r in range(1, len(self.analyte_tracker)):
            delta_counts = np.array(self.analyte_tracker[r][1:]) - np.array(self.analyte_tracker[r - 1][1:])
            row = self.analyte_tracker[r] + delta_counts.tolist()
            datatable.append(row)

        boiler_plates.section_header(title='Analyte Triage Progression', level=2)
        print(tabulate(datatable, headers=headers, tablefmt='simple'))

        return




    def filter_analytes(self, **kwargs):
        """
        Remove transitions, peptides, and proteins that do not meet requirements - ex min transitions per peptide or
        min peptides per protein.

        Currently works only for libraries to test impact of removing transitions or peptides.

        Note that this was created because Skyline fails to generate one of the important scores - the dot product of
        transition intensities vs. library intensities - when there is only one transition in a peak group.

        :return:
        """

        min_trans_per_pg = kwargs.get('min_trans_per_pg', 1)
        min_pep_per_prot = kwargs.get('min_pep_per_prot', 1)

        max_trans_per_pg = kwargs.get('max_trans_per_pg', None)
        max_pep_per_prot = kwargs.get('max_pep_per_prot', None)

        max_trans_method = kwargs.get('max_trans_method', 'intensity') # alternatives: 'low_variance'
        max_pep_method = kwargs.get('max_pep_method', 'intensity')


        # This method relies on having columns present that give the number of trans per peak group, pep per protein,
        # Make sure this is the most current info, call
        self.get_analytes_per_analyte()


        # todo - MAJOR: need to make this function work through higher levels including propagation of filtering
        # todo - cont. if already have pep_df and prot_df, need to filter those too. This gets messy w pepdf vs.precdf

        # ------------------------------------------
        # Transition/peak group filtering

        if min_trans_per_pg > 1:
            self.df.drop(self.df[self.df['trans_per_prec'] < min_trans_per_pg].index, inplace=True)
            self.count_analytes('min trans per peak group: ' + str(min_trans_per_pg))

            # todo - LOGIC GAP: need to track pep and protein present before and after as a result of this step and
            # todo - (cont) remove any pep and protein from those dfs if removed by this step. Right now, this is a
            # todo - (con) cause of the tables becoming inconsistent! There will be pep and prot in those tables that
            # todo - (cont) have no supporting transitions!

        # todo - later add filter on max trans maybe (faster to process if little quality value added)
        if max_trans_per_pg is not None:
            pass


        # ------------------------------------------
        # Peptide/protein filtering

        if min_pep_per_prot > 1:
            self.df.drop(self.df[self.df['pep_per_prot'] < min_pep_per_prot].index, inplace=True)
            self.count_analytes('min pep per prot: ' + str(min_pep_per_prot))


        # todo - later add filter on max pep maybe (value in faster to process even if little quality value added?)
        # alternative idea is to make pseudo proteins grouping peptides that act the same, but this could be dynamic
        if max_pep_per_prot is not None:
            pass

        # ------------------------------------------
        # output

        if self.verbose > 0 and (min_trans_per_pg > 1 or max_trans_per_pg > 1):
            boiler_plates.section_header(title='Filter Analytes', level=1)

            self.show_analyte_count_progression()


        self.tracker.capture('filter analytes')

        return




    def rollup(self, **kwargs):
        """
        Compute peptide and protein levels from current transition data

        :return:
        """

        method = kwargs.get('method', 'sum')
        ignore_provided_rollup = kwargs.get('ignore_provided_rollup', self.ignore_provided_rollup)


        if self.lib_or_results != 'results':
            return


        if self.verbose > 0:
            boiler_plates.section_header(title='Rollup to Peptide and Protein Level', level=1)

        # if not told to ignore provided rollup to pep and protein level in input and have it aleady, stop
        if ignore_provided_rollup is False and self.prot_df is not None and self.pep_df is not None:
            if self.verbose > 0:
                print('Peptide and protein data provided in input and not told to ignore it, so rollup NOT re-computed')
            return

        pep_meta_cols = ['seq_w_mods', 'sequence', 'missed_cleavages']
        prot_meta_cols = ['protein', 'protein_names', 'spike_status']

        # remove any missing cols
        pep_meta_cols = [x for x in pep_meta_cols if x in self.df.columns]
        prot_meta_cols = [x for x in prot_meta_cols if x in self.df.columns]

        # build final lists with relevant ID at front
        pep_meta_cols = ['pep_ID'] + prot_meta_cols + pep_meta_cols
        prot_meta_cols = ['prot_ID'] + prot_meta_cols

        if method == 'sum':
            # note: the min_count requirement of 1 returns nan if all nan vs. 0 if allowed to default to min_count=0
            self.pep_df = \
                self.df[pep_meta_cols + self.intensity_columns].groupby(pep_meta_cols).sum(min_count=1)

            self.prot_df = \
                self.df[prot_meta_cols + self.intensity_columns].groupby(prot_meta_cols).sum(min_count=1)

            # flatten the multi-indices created by the groupby operation
            self.pep_df.reset_index(inplace=True)
            self.prot_df.reset_index(inplace=True)

            # set the id columns as the indices and remove from the columns
            self.pep_df.set_index('pep_ID', drop=True, inplace=True)
            self.prot_df.set_index('prot_ID', drop=True, inplace=True)

        else:
            print('summary method unknown:', method)
            quit()


        # -------------------------------------------
        # rollup other data arrays (require different aggregation functions)

        for metric_cols in [self.q_value_columns, self.score_columns,
                            self.rt_columns, self.pred_rt_columns, self.delta_rt_columns]:

            # make sure cols are in df - handles absence of data
            cols_in_df = sum([1 for x in metric_cols if x in self.df.columns])

            if cols_in_df == len(metric_cols) and cols_in_df > 0:

                metric_df = self.df[pep_meta_cols + metric_cols].groupby(pep_meta_cols).mean()
                metric_df.reset_index(inplace=True)
                metric_df.set_index('pep_ID', drop=True, inplace=True)

                # add to main df
                self.pep_df[metric_cols] = metric_df[metric_cols]

        self.tracker.capture('rollup data levels')

        return





    # future method placeholders/ideas:

    def append(self, sw2):
        """
        This method takes another instance of SwathData and merges it into this instance to allow aggregation.
        """

        # todo - write append method

        return



    def workflow(self, wf_name):
        """
        This method takes a single workflow name argument and invokes a workflow progression of operations all done
        post-conversion (ex a series of filters, rollup, characterize). Envision using this to make it easy for anyone
        to invoke an exact sequence of processing without errors.
        """

        # todo - write workflow method

        return




    def save_library(self, filename):
        """
        This method dumps a library file based on the current transition table. Could also require a format arg. Would
        likely use mapping data in convert_format() method, so that should go somewhere more general.
        """

        # todo - write save library method

        return

    # ---------------- END CLASS ----------------------------------------------------




def strip_flanking(cell_string):
    """
    This is a worker function drive by pandas.apply. It removes flanking underscores from sequences
    :param cell_string:
    :return:
    """

    cell_string = cell_string.strip('_')

    return cell_string




def string_conv_blank_map(cell_string):
    """
    This is a worker function drive by pandas.apply. It forces the cell to str type and converts any nan to blank.
    :param cell_string:
    :return:
    """

    cell_string = str(cell_string).replace('nan', '')

    return cell_string




def strip_mol_group_tag(cell_string):

    cell_string = cell_string.replace('MoleculeGroup:/', '')

    return cell_string



def strip_wiff(cell_string):
    cell_string = str(cell_string)
    cell_string = cell_string.replace('.wiff', '')

    return cell_string



def strip_nl_mods(cell_string):
    """
    This is a worker method called by pandas apply. Its role is to remove mods that are not part of the TLC standard
    that are a result of Skyline's unusual way of dealing with neutral loss fragments where they are defined as
    mods. Mechanically, they must be set as fixed mods to allow the possibility of a NL fragment on any residue
    that could potentially lose the loss moiety. This results in absurd looking mod sequences. For example:

    D[Water Loss (D, E, S, T)]S[Water Loss (D, E, S, T)]N[Ammonia Loss (K, N, Q, R)]R[Ammonia Loss (K, N, Q, R)]PS[Water Loss (D, E, S, T)]GIPE[Water Loss (D, E, S, T)]R[Ammonia Loss (K, N, Q, R)]

    This method strips this down to the TLC standard version of this, which is just DSNRPSGIPER

    """


    water_nl = '[Water Loss (D, E, S, T)]'
    ammonia_nl = '[Ammonia Loss (K, N, Q, R)]'

    cell_string = cell_string.replace(water_nl, '').replace(ammonia_nl, '')

    return cell_string




def spectronaut_nl_frag_map(cell_string):
    """
    This is a worker function drive by pandas.apply. It converts Spectronaut's controlled vocabulary for fragment
    neutral loss types to std.
    :param cell_string:
    :return:
    """

    # EXTENSION_ZONE - may be missing some here!
    mapping_dict = {'noloss': '', 'H2O': '-H2O', 'NH3': '-NH3'}

    # if new value in dict, get it, otherwise use the same value
    cell_string = mapping_dict.get(cell_string, cell_string)

    return cell_string




def skyline_nl_frag_map(cell_string):
    """
    This is a worker function drive by pandas.apply. It converts Skyline's controlled vocabulary for fragment
    neutral loss types to std.
    :param cell_string:
    :return:
    """

    # EXTENSION_ZONE - may be missing some here!
    mapping_dict = {'H4COS': '-H4COS', 'H2O': '-H2O', 'NH3': '-NH3'}

    # if new value in dict, get it, otherwise use the same value
    cell_string = mapping_dict.get(cell_string, cell_string)

    return cell_string




def sciex_frag_type_map(cell_string):
    """
    This is a worker function drive by pandas.apply. It converts Sciex controlled vocabulary for fragment types which
    includes neutral loss types to std.
    :param cell_string:
    :return:
    """

    # EXTENSION_ZONE - may be missing some here!

    # if new value in dict, get it, otherwise use the same value
    cell_string = cell_string.replace(' - 18', '-H2O').replace(' - 17', '-NH3')

    return cell_string




def categorize_species_spike_status(cell_string, const_species, var_species):

    # const, var, unknown, mixed
    if const_species in cell_string:
        if var_species in cell_string:
            spike_status = 'mixed'
        else:
            spike_status = 'const'
    elif var_species in cell_string:
        spike_status = 'var'
    else:
        spike_status = 'unknown'

    return spike_status




def before_after_row_output(rows_before, rows_after, **kwargs):
    title = kwargs.get('title', None)
    removed_rows = (rows_before - rows_after)
    percent_removal = 100 * removed_rows / rows_before
    print()
    if title is not None:
        print(title)
    print('Rows before filter:   ', rows_before)
    print('Rows after filter:    ', rows_after)
    print('Rows removed:         ', removed_rows, '(', round(percent_removal, 2), '%)')
    print()

    return

