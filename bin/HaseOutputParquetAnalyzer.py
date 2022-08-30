import numpy as np
import timeit
from scipy import stats
import pandas as pd
import os
from collections import OrderedDict
import argparse
import glob
import math
import re
import pyarrow

parser = argparse.ArgumentParser(description = "Analyze HASE output parquet files.")

parser.add_argument('-i', '--InputFile', required = True,
                    help = "One or multiple input parquet arrays. Usage of wildcards is supported but then the argument has to be quoted.")

parser.add_argument('-o', '--OutputFile', type = str,
                    required = True,
                    help = "Path to tab-separated output file.")

parser.add_argument('-nc', '--NumberOfCovariates', type = int,
                    required = False,
                    help = "Number of covariates, used for calculating P-values. This is number of independent variables in the model (# of covariates + tested) - 1.")

parser.add_argument('-p', '--PThresh', type = float,
                    required = False,
                    help = "P-value threshold for filtering the results. Defaults None.", 
                    default = None)

parser.add_argument('-phen', '--PhenFilter', required = False, default = None, nargs = '+',
                    help = "One of three: first, file with the list of phenotypes to include. In that case, no header expected. Second, individual phenotype IDs specified and separated by space. Third, specified as None. If not specified, defaults to None.")

parser.add_argument('-snps', '--SnpFilter', required = False, default = None, nargs = '+',
                    help = "One of three: first, file with the list of SNPs/variants to include. In that case, no header expected. Second, individual SNP IDs specified and separated by space. Third, specified as None. If not specified, defaults to None.")

parser.add_argument('-sref', '--SnpRef', required = True,
                    help = "Reference for variant mapper, as specified in HASE documentation. Has to be gzipped and space-delimited.")

# parser.add_argument('-prec', '--RoundingPrecision', type = int,
#                     required = False,
#                     help = "Rounding precision for association statistics. If specified, output file is smaller but statistic precision is reduced. Defaults to None.", 
#                     default = None)

args = parser.parse_args()

# List of numpy arrays:
GeneFilter = args.PhenFilter
SnpFilter = args.SnpFilter
PThresh = args.PThresh
NumberOfCovariates = args.NumberOfCovariates
SnpRef = args.SnpRef

if ((GeneFilter is not None) and (str(GeneFilter).strip('[]"\'') != 'None')):
    if len(GeneFilter) == 1:
        print GeneFilter
        if re.match(r"[a-zA-Z0-9_.-/]*.txt$", str(GeneFilter).strip('[]"\'')):
            GeneFilter = str(GeneFilter).strip('[]"\'')
            GeneFilter = pd.read_csv(GeneFilter, header = None, delimiter = '\t')
            GeneFilter = set(GeneFilter.iloc[:,0].tolist())
    else:
        GeneFilter = GeneFilter

if ((SnpFilter is not None) and (str(SnpFilter).strip('[]"\'') != 'None')):
    if len(SnpFilter) == 1:
        print SnpFilter
        if re.match(r"[a-zA-Z0-9_.-/]*.txt$", str(SnpFilter).strip('[]"\'')):
            SnpFilter = str(SnpFilter).strip('[]"\'')
            SnpFilter = pd.read_csv(SnpFilter, header = None, delimiter = '\t')
            SnpFilter = set(SnpFilter.iloc[:,0].tolist())
    else:
        SnpFilter = SnpFilter


# Report some info
print "Input is: " + args.InputFile
print "Output is: " + args.OutputFile
if ((GeneFilter is not None) and (str(GeneFilter).strip('[]"\'') != 'None')):
    print "Output will be filtered to " + str(len(list(GeneFilter))) + " phenotypes."
if ((SnpFilter is not None) and (str(SnpFilter).strip('[]"\'') != 'None')):
    print "Output will be filtered to " + str(len(list(SnpFilter))) + " SNPs."
if ((PThresh is not None) and (str(PThresh).strip('[]"\'') != 'None')):
    print "Output will be filtered to association P<" + str(PThresh)

input_parquet = glob.glob(args.InputFile)

if len(input_parquet) == 0:
    print("Input path is empty!")
    exit()

# Function for P-value calculation
def get_p_value(t_stat, df = None):
    if df is None:
        return stats.norm.sf(np.abs(t_stat))*2
    else:
        return stats.t.sf(np.abs(t_stat), df)*2


# Read in reference
ref = pd.read_csv(SnpRef, compression = 'gzip', sep = ' ')
ref = ref.loc[:, ['ID', 'str_allele1', 'str_allele2']]
if ((SnpFilter is not None) and (str(SnpFilter).strip('[]"\'') != 'None')):
    ref = ref[ref['ID'].isin(list(SnpFilter))]

print "SNP reference loaded."

with open(os.path.join(args.OutputFile), 'w') as f:
    f.write('phenotype\tSNP\tref_all\talt_all\tbeta\tse\tP\ti_squared\tN\n')
    for filename in input_parquet:
        print "Processing file: " + filename
        df = pd.read_parquet(filename)

        df['df'] = df['sample_size'] - NumberOfCovariates - 2

        # Filtering
        if ((GeneFilter is not None) and (str(GeneFilter).strip('[]"\'') != 'None')):
            print "Gene filter active." 
            df = df[df['phenotype'].isin(list(GeneFilter))]
        if ((SnpFilter is not None) and (str(SnpFilter).strip('[]"\'') != 'None')):
            print "SNP filter active."
            df = df[df['SNP'].isin(list(SnpFilter))]

        # Calculate P-values and betas
        T = df['beta']/df['standard_error']
        Pvals = get_p_value(T, df['df'])


        # Put into pandas:
        inp = OrderedDict()
        inp["phenotype"] = df["phenotype"]
        inp["SNP"] = df["variant"]
        inp["beta"] = df['beta']
        inp["se"] = df["standard_error"]
        inp["P"] = Pvals
        inp["i_squared"] = df["i_squared"]
        inp["N"] = df["sample_size"]

        df_out = pd.DataFrame(inp)

        if ((PThresh is not None) and (str(PThresh).strip('[]"\'') != 'None')):
            print "P-value filter: " + str(PThresh)
            df_out = df_out[df_out['P'] < PThresh]


        temp_ref = ref[ref['ID'].isin(list(df_out['SNP']))]
        df_out = pd.merge(df_out, temp_ref, left_on = 'SNP', right_on = 'ID')
        df_out = df_out.rename(columns = {'str_allele1' : 'ref_all', 'str_allele2' : 'alt_all'})

        df_out = df_out[['phenotype', 'SNP', 'ref_all', 'alt_all', 'beta', 'se', 'P', 'i_squared', 'N']]

        if len(df_out.index) > 0:
            df_out.to_csv(f, sep = "\t", header = False, index = None)
        else:
            print "No rows to include from this array."

print "Output written."
