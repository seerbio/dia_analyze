import os.path
import numpy as np
import pandas as pd

from crozes_base.helpers import helper_functions
from crozes_base.io import boiler_plates




def get_mod_catalog(**kwargs):
    """
    This function gets a modification catalog necessary for mapping answers across dialects. It first tries to get a
    previously saved mod catalog stored as a pkl file in an expected project reference files folder that is hardcoded
    relative to the location of this module. If this file has not be generated or is not present for whatever reason,
    it looks for the Sciex Excel version of their mod catalog in the same folder and reads it to create this mod
    dictionary and then also saves a pkl version of this file for faster function next time.

    :param kwargs:
    :return:
    """

    verbose = kwargs.get('verbose', 0)

    parent_dir = os.path.dirname(__file__)

    references_path = os.path.join(parent_dir, 'reference_files')

    mod_pkl = os.path.join(references_path, 'mod_catalog.pkl')
    mod_xlsx = os.path.join(references_path, 'Sciex_Unified_Modification_Catalog.xlsx')


    # read in pickled mod catalog if present
    if os.path.exists(mod_pkl):

        mod_dict = helper_functions.read_pickled_object(mod_pkl)

        if verbose > 1:
            print()
            print('mod catalog read from prior pkl file')
            print()

    else:

        if verbose > 1:
            print()
            print('creating mod catalog from Sciex excel mod catalog...')
            print()

        # read excel catalog to dictionary
        if os.path.exists(mod_xlsx):

            # read excel catalog
            mod_df = pd.read_excel(mod_xlsx, sheet_name='Modification Catalog Table', skiprows=[0])

            # drop empty rows
            mod_df.dropna(how='all', inplace=True)

            useful_columns = {'PSI MS Name    (Standard name, except where in red)': 'psi_name',
                              'Site': 'site',
                              'Position': 'term_specificity',
                              'UniMod: Classification': 'unimod_class',
                              'TLC in AB SCIEX Data Dictionary': 'tlc'}

            # reduce to relevant columns
            mod_df = mod_df[useful_columns.keys()]

            # rename columns to more manageable names
            mod_df = mod_df.rename(columns=useful_columns)

            # remove substitutions
            mod_df = mod_df[mod_df['unimod_class'] != 'AA substitution']

            # remove rows with no three letter code
            mod_df = mod_df[pd.notnull(mod_df['tlc'])]

            # get the name spectronaut would call each mod
            mod_df['spectronaut_name'] = mod_df.apply(lambda row: get_spectronaut_name(row), axis=1)

            # create two dictionaries from excel table
            mod_dict = mod_df.set_index('spectronaut_name').to_dict()['tlc']

            # write out as pkl file for future faster use
            helper_functions.save_object(mod_dict, mod_pkl)

        else:
            mod_dict = None
            print('FATAL ERROR!!: No modification catalog file present in reference_files folder')
            quit()


    return mod_dict




def get_spectronaut_name(row):
    """
    This is a worker function drive by pandas.apply.

    Returns the modification name as it would appear in Spectronaut results.

    :param row: a row in a dataframe with expected column headers
    :return:
    """

    if row['term_specificity'] == 'Anywhere':
        name = str(row['psi_name']) + ' (' + str(row['site']) + ')'

    # ----------------
    # Protein N-term cases

    elif row['term_specificity'] == 'Protein N-term' and row['site'] == 'N-term':
        name = str(row['psi_name']) + ' (Protein N-term)'

    # todo - make sure this case is not treated differently - ie does Spectronaut syntax include the residue specificity?
    elif row['term_specificity'] == 'Protein N-term' and row['site'] != 'N-term':
        name = str(row['psi_name']) + ' (Protein N-term)'

    # ----------------
    # Protein C-term cases

    elif row['term_specificity'] == 'Protein C-term' and row['site'] == 'C-term':
        name = str(row['psi_name']) + ' (Protein C-term)'

    # todo - make sure this case is not treated differently - ie does Spectronaut syntax include the residue specificity?
    elif row['term_specificity'] == 'Protein C-term' and row['site'] != 'C-term':
        name = str(row['psi_name']) + ' (Protein C-term)'

    # ----------------
    # Peptide N-term cases

    elif row['term_specificity'] == 'Any N-term' and row['site'] == 'N-term':
        name = str(row['psi_name']) + ' (N-term)'

    # todo - make sure this case is not treated differently - ie does Spectronaut syntax include the residue specificity?
    elif row['term_specificity'] == 'Any N-term' and row['site'] != 'N-term':
        name = str(row['psi_name']) + ' (N-term)'

    # ----------------
    # Peptide C-term cases

    elif row['term_specificity'] == 'Any C-term' and row['site'] == 'C-term':
        name = str(row['psi_name']) + ' (C-term)'

    # todo - make sure this case is not treated differently - ie does Spectronaut syntax include the residue specificity?
    elif row['term_specificity'] == 'Any C-term' and row['site'] != 'C-term':
        name = str(row['psi_name']) + ' (C-term)'

    else:
        name = None
        print('FATAL ERROR!! Could not compute spectronaut name for the following input row:')
        print(row)
        quit()


    return name




def determine_mod_syntax(mod_series):
    """

    :param mod_series: a pandas series - the column with modified peptide sequences from a dataframe
    :return: mod_syntax: a string indicating the mod syntax present in the input data
    """

    # initialize
    mod_syntax = 'unknown'

    # controlled vocab for syntax names: 'TLC', 'PSI-Mod', 'integer', 'unknown'

    # determination approach assumes Cysteine are going to be CAM alkylated
    have_counts = False
    test_len = min(1000, mod_series.shape[0])

    while have_counts is False:

        # make a list of first rows in series
        sample_list = mod_series.iloc[:test_len].tolist()

        TLC_count = sum([1 for x in sample_list if '[CAM]' in x])
        long_count = sum([1 for x in sample_list if '[Carbamidomethyl (C)]' in x])
        int_count = sum([1 for x in sample_list if '58' in x])

        if TLC_count > 0 or long_count > 0 or int_count > 0:
            have_counts = True
        elif test_len == mod_series.shape[0]:
            have_counts = True  # don't have counts, but checked whole list so stop while loop and crash next
        else:
            test_len = min(test_len * 10, mod_series.shape[0])

    if TLC_count > 0 and long_count == 0 and int_count == 0:
        mod_syntax = 'TLC'
    elif TLC_count == 0 and long_count > 0 and int_count == 0:
        mod_syntax = 'PSI-Mod'
    elif TLC_count == 0 and long_count == 0 and int_count > 0:
        mod_syntax = 'integer'
    else:
        print('FATAL ERROR!! Mod syntax cannot be determined')
        print('The following counts were measured in a sample of', test_len, 'rows.')
        print('TLC:     ', TLC_count)
        print('PSI-Mod: ', long_count)
        print('Integer: ', int_count)
        quit()

    return mod_syntax




def get_mod_freqs(df, col_name, mod_system, **kwargs):
    """
    This function assumes the input is fully converted to the specified mod_system. For example, if set to 'TLC' system,
    the peptides must both use the TLC codes and also the positioning of terminal mods in the Sciex TLC syntax. See
    the long_mod_to_tlc() method of DiaData class for details.

    :param df:
    :param col_name:
    :param mod_system: Only 'TLC' supported for now.
    :param kwargs:
    :return:
    """

    title = kwargs.get('title', 'Modification Counts')
    verbose = kwargs.get('verbose', 0)

    if verbose > 0:
        boiler_plates.section_header(title=title, level=2)

    mod_counts_dict = {}


    if mod_system == 'TLC':

        for seq in df[col_name]:

            _, mod_list = get_seq_and_mods_from_combined(seq)

            # add to counting dict
            for mod in mod_list:
                if mod in mod_counts_dict:
                    mod_counts_dict[mod] += 1
                else:
                    mod_counts_dict.update({mod: 1})

    else:

        print('mod syntax argument unknown or unsupported:', mod_system)
        quit()

    # todo - change this to output a nice table sorted by descending counts
    if verbose > 0:
        for mod in mod_counts_dict:
            print('count: ', mod_counts_dict[mod], 'mod: ', mod, )


    return mod_counts_dict




def get_seq_and_mods_from_combined(seq):

    mod_list = []

    if '[' in seq:

        aa_sequence = ''

        # convert the sequence (1) so that can be split easily w '[' both before and after mods
        #                      (2) also and mark mods w '%' so can tell mod string apart from sequence string
        seq_conv = seq.replace('[', '[%').replace(']', '[')
        parts = seq_conv.split('[')
        parts = [x for x in parts if len(x) > 0]

        for p in range(len(parts)):
            part = parts[p]

            # find the mods (tagged with %)
            if '%' in part:

                # unambiguously indicated N-terminal mod case per Sciex TLC syntax - ex. [1Ac]-XXXXX
                if p == 0:
                    if parts[p + 1][0] == '-':
                        mod = '[' + part[1:] + ']-'
                    else:
                        mod = None
                        print('Mod case w unexpected n-term mod syntax!')
                        print(seq)
                        quit()

                # C-terminal cases
                elif p == len(parts) - 1:
                    # these two are the same I think, but need to check syntax
                    # - ie if side chain mod, should not use - even when on C-term res
                    # no! they aren't the same but the else clause would cover it all

                    # unambiguous case - has a hyphen indicating terminal, not side chain
                    # (check if the prior part is a hyphen)
                    if parts[p - 1][-1] == '-':
                        mod = '-[' + part[1:] + ']'

                    # ambiguous case - need to look up in mod catalog to know if side chain or terminal
                    # - side chain mods coincidentally on the C-term denoted XXXK[1Ac] where "K[1Ac]" is the mod
                    # - true terminal mods (MUST be on the term) denoted XXXK-[Amm] where "-[Amm]" is the mod
                    else:
                        mod = parts[p - 1][-1] + '[' + part[1:] + ']'
                else:
                    mod = parts[p - 1][-1] + '[' + part[1:] + ']'

                mod_list.append(mod)

            else:
                aa_sequence += part

    else:
        aa_sequence = seq

    return aa_sequence, mod_list
