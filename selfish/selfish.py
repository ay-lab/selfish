#!/usr/bin/env python3
import argparse
import os
import sys
import re
import math
from collections import defaultdict

import numpy as np
from scipy.stats import norm
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import straw
import cooler
import pandas as pd
import statsmodels.stats.multitest as smm

# Read the version file and get the value.
dir = os.path.dirname(__file__)
version_py = os.path.join(dir, "_version.py")
exec(open(version_py).read())


# Set the arguments required for the program.
def parse_args(args):
    parser = argparse.ArgumentParser(description="Check the help flag")

    parser.add_argument("-f1",
                        "--file1",
                        dest="f1_path",
                        help="REQUIRED: Contact map 1",
                        required=False)

    parser.add_argument("-f2",
                        "--file2",
                        dest="f2_path",
                        help="REQUIRED: Contact map 2",
                        required=False)

    parser.add_argument("-o",
                        "--outfile",
                        dest="outdir",
                        help="REQUIRED: Name of the output file.\
                       Output is a numpy binary.",
                        required=True)

    parser.add_argument("-c",
                        "--changefile",
                        dest="changedir",
                        help="Output file with the contact map fold changes",
                        required=False,
                        default="")
    parser.add_argument(
                        "-st",
                        "--sparsityThreshold",
                        dest="st",
                        type=float,                        
                        help="OPTIONAL: Selfish filters out interactions in sparse areas. The default value is 0.7.",
                        required=False)
    parser.add_argument("-ns",
                        '--no-sparsityCheck',
                         dest='st',
                         action='store_false',
                         required=False,
                         help="OPTIONAL: Do not apply any sparsity checking.")
    parser.set_defaults(st=0.7)
    parser.add_argument(
        "-norm",
        "--normalization",
        default=False,
        dest="norm_method",
        help="RECOMMENDED: Normalization method (KR, VC,...).",
        required=False)
    parser.add_argument(
        "-ch",
        "--chromosome",
        dest="chromosome",
        help="REQUIRED: Specify which chromosome to run the program for.",
        default='n',
        required=True)
    parser.add_argument(
        "-ch2",
        "--chromosome2",
        dest="chromosome2",
        help="Optional: Specify the second chromosome for interchromosomal analysis.",
        default='n',
        required=False)

    parser.add_argument(
        "-t",
        "--tsvout",
        dest="tsvout",
        help="Should the program output a tab-seperated output rather than the numpy array. Value must be the p-value cutoff (i.e. 0.05)",
        default=0,
        required=False)

    parser.add_argument("-r",
                        "--resolution",
                        dest="resolution",
                        help="REQUIRED: Resolution used for the contact maps",
                        required=True)

    parser.add_argument("-bed1", "--bed1", dest="bed1",
                        help="BED file 1 for HiC-Pro type input",
                        default="",
                        required=False)
    parser.add_argument("-bed2", "--bed2", dest="bed2",
                        help="BED file 2 for HiC-Pro type input",
                        default="",
                        required=False)
    parser.add_argument("-m1", "--matrix1", dest="mat1",
                        help="MATRIX file 1 for HiC-Pro type input",
                        default="",
                        required=False)
    parser.add_argument("-m2", "--matrix2", dest="mat2",
                        help="MATRIX file 2 for HiC-Pro type input",
                        default="",
                        required=False)

    parser.add_argument("-b1", "--biases1", dest="biasfile1",
                        help="RECOMMENDED: biases calculated by\
                        ICE or KR norm for each locus for contact map 1 are read from BIASFILE",
                        required=False)

    parser.add_argument("-b2", "--biases2", dest="biasfile2",
                        help="RECOMMENDED: biases calculated by\
                        ICE or KR norm for each locus for contact map 2 are read from BIASFILE",
                        required=False)

    parser.add_argument(
        "-sz",
        "--sigmaZero",
        dest="s_z",
        type=float,
        default=1.6,
        help="OPTIONAL: sigma0 value for the method. DEFAULT is 1.6. \
        Experimentally chosen for 5Kb resolution",
        required=False)

    parser.add_argument("-i", "--iterations", dest="s", default=10,
                        type=int,
                        help="OPTIONAL: iteration count for the method. \
                        DEFAULT is 10. Experimentally chosen for \
                        5Kb resolution",
                        required=False)

    parser.add_argument(
        "-d",
        "--distanceFilter",
        dest="distFilter",
        help="OPTIONAL: If the data is too sparse for distant \
        locations(ie. distance > 5Mb), you can filter them. \
        DEFAULT is None",
        required=False)

    parser.add_argument("-v",
                        "--verbose",
                        dest="verbose",
                        type=bool,
                        default=True,
                        help="OPTIONAL: Verbosity of the program",
                        required=False)
#    parser.add_argument("-mn",
#                        "--mutual",
#                        dest="mutual_nz",
#                        type=bool,
#                        default=True,
#                        help="OPTIONAL: Only report interactions that are non-zero in both contact maps.",
#                        required=False)
    parser.add_argument("-nm",
                        '--no-mutual',
                         dest='mutual_nz',
                         action='store_false',
                         required=False,
                         help="OPTIONAL: Consider also interactions that are non-zero in one of the contact maps.")
    parser.set_defaults(mutual_nz=True)
    #parser.add_argument("-nb",
    #                    '--no-balance',
    #                     dest='cooler_do_balance',
    #                     action='store_false',
    #                     required=False,
    #                     help="OPTIONAL: Set if the cooler data was normalized prior to creating the .cool file.")
    parser.add_argument("-p",
                        "--plot",
                        dest="plot",
                        help="OPTIONAL: use this flag for generating plots. \
                      DEFAULT is False.",
                        default=False,
                        required=False)
    parser.add_argument(
        "-lm",
        "--lowmem",
        dest="lowmemory",
        help="OPTIONAL: Uses Float32 instead of Float64, halves the memory usage. Default is False",
        default=False,
        required=False)

    parser.add_argument("-V", "--version", action="version",
                        version="Selfish {}".format(__version__), help="Print version and exit")

    return parser.parse_args()


def read_map_pd(path, res, bias, dt, ch):
    """
    :param path: File path for the contact map
    :param res: Resolution of the contact map
    :param bias: Bias dictionary to multiply the values
    :param dt: Data type to use as the matrix, np.float64 OR np.float32
    :return: A square matrix with contact counts as values.
    """
    sep = get_sep(path)
    o = pd.read_csv(path, sep=sep, header=None)
    o.dropna(inplace=True)
    o = o[np.vectorize(isChr)(o[0], ch)]
    o = o[np.vectorize(isChr)(o[2], ch)]
    o[1] //= res
    o[3] //= res
    n = max(max(o[1]), max(o[3])) + 1
    out = np.zeros((n, n), dtype=dt)
    out[o[1], o[3]] = o[4]
    out[o[3], o[1]] = o[4]
    if bias:
        factors = np.vectorize(bias.get)(o[1], 1)
        o[4] = np.divide(o[4], factors)
        factors = np.vectorize(bias.get)(o[3], 1)
        o[4] = np.divide(o[4], factors)

    np.nan_to_num(o, copy=False, nan=0, posinf=0, neginf=0)
    return np.triu(out)


def get_diags(map):
    """
    :param map: Contact map, numpy matrix
    :return: 2 Dictionaries where keys are the diagonal number and values are the mean of that diagonal in one dictionary and the std. in the other dictionary.
    """
    means = {}
    stds = {}
    for i in range(len(map)):
        diag = map.diagonal(i)
        diag = diag[diag != 0]
        if len(diag) > 0:
            mean = np.mean(diag)
            std = np.std(diag) if np.std(diag) != 0 else 1
            if math.isnan(mean):
                means[i] = 0
            else:
                means[i] = mean
            if math.isnan(std):
                stds[i] = 1
            else:
                stds[i] = std
        else:
            means[i] = 0
            stds[i] = 1
    return means, stds


def normalize_map(cmap):
    """
    :param cmap: contact map
    :param non_zero: same size matrix with map that has True as value for indices to be normalized.
    :return: Changes matrix in place to have each diagonal sum = 0
    """
    means, stds = get_diags(cmap)
    x, y = np.nonzero(cmap)
    distance = np.abs(x - y)
    m = np.vectorize(means.get)(distance)
    s = np.vectorize(stds.get)(distance)
    cmap[x, y] -= m
    cmap[x, y] /= s

    x, y = np.where(cmap == 0)
    distance = np.abs(x-y)
    m = np.vectorize(means.get)(distance)
    s = np.vectorize(stds.get)(distance)
    cmap[x, y] -= m
    cmap[x, y] /= s
    np.nan_to_num(cmap, copy=False, posinf=0, neginf=0, nan=0)


def inter_normalize_map(cmap):
    vals = cmap[cmap > 0]
    m = np.mean(vals)
    s = np.std(vals)
    cmap -= m
    cmap /= s
    np.nan_to_num(cmap, copy=False, nan=0, posinf=0, neginf=0)


def get_sep(f):
    """
    :param f: file path
    :return: Guesses the value separator in the file.
    """
    with open(f) as file:
        for line in file:
            if "\t" in line:
                return '\t'
            if " " in line.strip():
                return ' '
            if "," in line:
                return ','
            break
    raise FileNotFoundError


def read_bias(f, chr, res):
    """
    :param f: Path to the bias file
    :return: Dictionary where keys are the bin coordinates and values are the bias value to multiply them with.
    """
    d = defaultdict(lambda: 1.0)
    if f:
        sep = get_sep(f)
        with open(f) as file:
            for line in file:
                line = line.strip().split(sep)
                if isChr(line[0], chr):
                    val = float(line[2])
                    if not np.isnan(val):
                        d[float(line[1]) // res] = val
                    else:
                        d[float(line[1]) // res] = np.Inf
        return d
    return False


def DCI(f1,
        f2,
        bed1='',
        bed2='',
        res=100000,
        sigma0=1.6,
        s=10,
        plot_results=False,
        changes="",
        verbose=True,
        mutual_nz=True,
        #cooler_do_balance=True,
        distance_filter=5000000,
        bias1=False,
        bias2=False,
        chromosome='n',
        chromosome2=None,
        tsvout=None,
        st=False,
        low_memory=False,
        norm_method= None):
    """
    :param f1: Path to contact map 1
    :param f2: Path to contact map 2
    :param res: Resolution of the contact maps
    :param sigma0: Sigma0 parameter of the SELFISH method
    :param s: Iteration count parameter of the SELFISH method
    :param plot_results: Boolean parameter to draw plots of the files
    :param verbose: Boolean parameter to make program print information during execution
    :param distance_filter: Distance where interactions that are further would be discarded
    :param bias1: Path to a bias file 1
    :param bias2: Path to a bias file 2
    :param chromosome: Chromosome for which the program will run
    :param low_memory: Whether to halve the memory usage by using 32 bit precision instead of 64
    :return: Matrix of p-values
    """
    #if not chromosome2 or chromosome2 == 'n':
    #    chromosome2 = chromosome

    if (chromosome != chromosome2) and not ((('.hic' in f1) or ('.cool' in f1) or ('.mcool' in f1)) and (('.hic' in f2) or ('.cool' in f2) or ('.mcool' in f2))):
        print(
            "Interchromosomal analysis is only supported for .hic and .cool input formats.")
        raise FileNotFoundError

    if plot_results:
        import seaborn as sns
    scales = [(sigma0 * (2**(i / s))) for i in range(1, s + 3)]

    dt = np.float32 if low_memory else np.float64

    biasDict1 = read_bias(bias1, chromosome, res)
    biasDict2 = read_bias(bias2, chromosome, res)

    if type(f1) == str and type(f2) == str:
        if verbose:
            print("Reading Contact Map 1")
        if bed1 != '' and f1 != '':
            if chromosome == 'n':
                print("You need to specify a chromosome for matrix files.")
                raise FileNotFoundError
            a = readBEDMAT(bed1, f1, res, chromosome, biasDict1)
        elif f1.endswith('.hic'):
            if chromosome == 'n':
                print("You need to specify a chromosome for hic files.")
                raise FileNotFoundError
            a = readHiCFile(f1, chromosome, chromosome2,
                            res, distance_filter // res, norm_method)
        elif f1.endswith('.cool'):  # cooler.fileops.is_cooler(f1):
            a = readCoolFile(f1, chromosome, chromosome2, norm_method)
        elif f1.endswith('.mcool'):  # cooler.fileops.is_multires_file(f1):
            a = readMultiCoolFile(f1, chromosome, chromosome2, res, norm_method)
        else:
            a = read_map_pd(f1, res, biasDict1, dt, chromosome)

        if verbose:
            print("Reading Contact Map 2")
        if bed2 != '' and f2 != '':
            if chromosome == 'n':
                print("You need to specify a chromosome for matrix files.")
                raise FileNotFoundError
            b = readBEDMAT(bed2, f2, res, chromosome, biasDict2)
        elif f2.endswith('.hic'):
            if chromosome == 'n':
                print("You need to specify a chromosome for hic files.")
                raise FileNotFoundError
            b = readHiCFile(f2, chromosome, chromosome2,
                            res, distance_filter // res, norm_method)
        elif f2.endswith('.cool'):  # cooler.fileops.is_cooler(f1):
            b = readCoolFile(f2, chromosome, chromosome2, norm_method)
        elif f2.endswith('.mcool'):  # cooler.fileops.is_multires_file(f1):
            b = readMultiCoolFile(f2, chromosome, chromosome2, res, norm_method)
        else:
            b = read_map_pd(f2, res, biasDict2, dt, chromosome)

        f1 = f1.split('.')[0] if '.' in f1 else f1
        f2 = f2.split('.')[0] if '.' in f2 else f2
    elif type(f1) == np.ndarray and type(f2) == np.ndarray \
            and f1.shape[0] == f1.shape[1] \
            and f2.shape[0] == f2.shape[1]:
        a = f1.copy()
        f1 = "Map 1"
        b = f2.copy()
        f2 = "Map 2"
    else:
        print(
            "Error: inputs should either be file names or square numpy matrices"
        )
        return

    if chromosome == chromosome2:
        if verbose:
            print("Applying distance filter")
        if distance_filter > 0:
            a = np.tril(a, distance_filter // res)
            b = np.tril(b, distance_filter // res)
    
    temp_x = max(a.shape[0],b.shape[0])
    temp_y = max(a.shape[1],b.shape[1])
	
    if chromosome == chromosome2:
	    temp_x = max(temp_x,temp_y)
	    temp_y = temp_x
    
    #tmp_n = max(max(a.shape), max(b.shape))
    tmp = np.zeros((temp_x,temp_y))
    tmp[:a.shape[0], :a.shape[1]] = a
    a = tmp.copy()
    tmp = np.zeros((temp_x, temp_y))
    tmp[:b.shape[0], :b.shape[1]] = b
    b = tmp.copy()
    tmp = None
    
    if mutual_nz:
        non_zero_indices = np.logical_and(a != 0, b != 0)
    else:
        non_zero_indices = np.logical_or(a != 0, b != 0)
    #non_zero_indices1 = (a != 0)
    #non_zero_indices2 = (b != 0)
       
    if plot_results:
        plt.clf()
        p_f = a.copy()
        p_f[p_f == 0] = 1
        sns.heatmap(np.abs(np.log10(p_f)))
        plt.title("Contact counts (KR normalised, log scale) " + f1)
        plt.savefig('f1_1_kr.png')
        plt.clf()
        p_f = b.copy()
        p_f[p_f == 0] = 1
        sns.heatmap(np.log10(p_f))
        plt.title("Contact counts (KR normalised, log scale) " + f2)
        plt.savefig('f2_2_kr.png')

    # np.save("kr_a.dat",a)
    # np.save("kr_b.dat",b)    

    changes_array = np.zeros_like(a)
    changes_array[non_zero_indices] = np.log2(
        np.divide(a + 1, b + 1))[non_zero_indices]
    if changes != "":
        np.save(changes, changes_array)

    if chromosome == chromosome2:
        if verbose:
            print("Diagonal Normalizing Map 1")
        normalize_map(a)
        if verbose:
            print("Diagonal Normalizing Map 2")
        normalize_map(b)
    else:
        if verbose:
            print("Normalizing Map 1")
        inter_normalize_map(a)
        if verbose:
            print("Normalizing Map 2")
        inter_normalize_map(b)
    # normalize the zero interactions using calculated means and sigmas for each contact map

    if plot_results:
        plt.clf()
        sns.heatmap(np.abs(a))
        plt.title("Contact counts (Diagonal normalised) " + f1)
        plt.savefig('f1_1_diag.png')
        plt.clf()
        sns.heatmap(np.abs(b))
        plt.title("Contact counts (Diagonal normalised) " + f2)
        plt.savefig('f2_2_diag.png')

        plt.clf()
        sns.heatmap(np.abs(np.abs(a) - np.abs(b)))
        plt.title("Contact difference (Diagonal normalised)")
        plt.savefig('f12_diff_diag.png')

    # assign zero to the first four diagonals
    il = np.tril_indices(a.shape[0], 2)
    a[il] = 0
    il = np.tril_indices(b.shape[0], 2)
    b[il] = 0
    # np.save("diag_a.dat",a)
    # np.save("diag_b.dat",b)
    if verbose:
        print("Applying gaussians")

    final_p = np.ones(len(a[non_zero_indices]))
    final_scale = np.ones(len(a[non_zero_indices]))

    size = a.shape[0]
    b -= a
    #a = None
    diff = b.copy()
    #b = None
    np.nan_to_num(diff, copy=False, posinf=0, neginf=0, nan=0)
    d_pre = gaussian_filter(diff, scales[0])
    count = 1
    for scale in scales[1:]:
        if verbose:
            print("Gaussian:", count, "of", len(scales) - 1)
        d_post = gaussian_filter(diff, scale)
        d_diff = d_post - d_pre
        params = norm.fit(d_diff[non_zero_indices])
        p_vals = norm.cdf(d_diff[non_zero_indices],
                          loc=params[0],
                          scale=params[1])
        np.nan_to_num(p_vals, copy=False, posinf=1, neginf=1, nan=1)
        p_vals[p_vals > 0.5] = 1 - p_vals[p_vals > 0.5]
        p_vals *= 2
        final_scale[p_vals < final_p] = scale
        final_p[p_vals < final_p] = p_vals[p_vals < final_p]        
        d_pre = d_post.copy()
        count += 1

    d_diff = None
    _, out_p = smm.multipletests(final_p, method='fdr_bh')[:2]

    o = np.ones(diff.shape, dtype=dt)
    o[non_zero_indices] = out_p
    if plot_results:                               
        thr = 10
        o[o == 0] = 10**-1*thr
        o = o[:temp_x,:temp_y]

        plt.clf()
        o = np.abs(np.log(o))
        #params = norm.fit(d_diff[non_zero_indices])

        o[o > thr] = thr
        o[np.isinf(o)] = thr  # log(0) = -inf
        sns.heatmap(np.abs(o))
        # sns.heatmap(np.abs(np.log10(o)))
        plt.title("Differential analysis")
        plt.savefig("f1f2_selfish.png",dpi=200)                    
    
    if st:         
        index = o < tsvout
        sig_2d_indices = np.argwhere(index)
        q_val = o[index]
        osc = np.zeros(o.shape, dtype=dt)
        osc[non_zero_indices] = final_scale 
        scale = osc[index]
        changes_array = changes_array[index]
        x = sig_2d_indices[:,0] #[x[0] for x in sig_2d_indices]
        y = sig_2d_indices[:,1] #[x[1] for x in sig_2d_indices] 
        print("number of differential interactions before sparsity check: ",len(x)) 
        os = o < 1  #np.zeros(diff.shape, dtype=dt)   
        #os[non_zero_indices] = 1                  
        nonsparse = x != 0
        #print(f"ns={nonsparse}")
        for i in range(len(x)):
            s = math.ceil(scale[i])
            c1 = np.sum(os[x[i]-s:x[i]+s+1, y[i]-s:y[i]+s+1]) / \
            ((2*s+1)**2)
            s = 2*s
            c2 = np.sum(os[x[i]-s:x[i]+s+1, y[i]-s:y[i]+s+1]) / \
            ((2*s+1)**2)
            if c1 < st or c2 < 0.6:
                nonsparse[i] = False

        x = x[nonsparse]
        y = y[nonsparse]
        q_val= q_val[nonsparse]
        diff = None
        print("number of differential interactions after sparsity check: ",len(x))
        return _, x, y, q_val, changes_array
 
    diff = None
    return  o, _, _, _, changes_array


def sorted_indices(dci_out):
    """
    :param dci_out: Output of DCI method
    :return: One array for X coordiantes and another for Y coordinates of sorted p-values.
    """
    ind = np.unravel_index(np.argsort(dci_out, axis=None), dci_out.shape)
    return ind[0], ind[1]


def toNearest(n, res):
    """
    :param n: Number
    :param res: Resolution
    :return: Number rounded up to the closest multiplicate of res
    """
    if n % res == 0:
        return n
    return (n // res + 1) * res


def isChr(s, c):
    """
    :param s: String
    :param c: Chromosome, number or X or Y
    :return: Whether s matches c
    """
    if 'X' == c:
        return 'X' in c
    if 'Y' == c:
        return 'Y' in c
    return str(c) in re.findall("[1-9][0-9]*", s)


def readHiCFile(f, chr, chr2, res, distance, norm_method):
    """
    :param f: .hic file path
    :param chr: Which chromosome to read the file for
    :param res: Resolution to extract information from
    :return: Numpy matrix of contact counts
    """
    if not norm_method:
        result = straw.straw('KR', f, str(chr), str(chr2), 'BP', res)
    else:
        result = straw.straw(str(norm_method), f, str(chr), str(chr2), 'BP', res)
    x = np.array(result[0]) // res
    y = np.array(result[1]) // res
    val = np.array(result[2])
    #n = max(max(x), max(y)) + 1
    n1 = max(x) + 1
    n2 = max(y) + 1
    o = np.zeros((n1, n2))
    o[x, y] = val
    #o[y, x] = val
    np.nan_to_num(o, copy=False, nan=0, posinf=0, neginf=0)
    return o


def readCoolFile(f, chr, chr2, norm_method):
    """
    :param f: .cool file path
    :param chr: Which chromosome to read the file for
    :return: Numpy matrix of contact counts
    """
    clr = cooler.Cooler(f)
    if chr == chr2:
        if not norm_method:
            result = clr.matrix(balance=True).fetch(chr)
        else:
            result = clr.matrix(balance=norm_method).fetch(chr)
    else:
        if not norm_method:
            result = clr.matrix(balance=True).fetch(chr, chr2)
        else:
            result = clr.matrix(balance=norm_method).fetch(chr, chr2)
    result[np.isnan(result)] = 0
    return result


def readMultiCoolFile(f, chr, chr2, res, norm_method):
    """
    :param f: .cool file path
    :param chr: Which chromosome to read the file for
    :param res: Resolution to extract information from
    :return: Numpy matrix of contact counts
    """
    uri = '%s::/resolutions/%s' % (f, res)
    clr = cooler.Cooler(uri)
    if chr == chr2:
        if not norm_method:
            result = clr.matrix(balance=True).fetch(chr)
        else:
            result = clr.matrix(balance=norm_method).fetch(chr)
    else:
        if not norm_method:
            result = clr.matrix(balance=True).fetch(chr, chr2)
        else:
            result = clr.matrix(balance=norm_method).fetch(chr)
    result[np.isnan(result)] = 0
    return result


def readBEDMAT(bed, mat, res, chr, bias):
    """
    :param f: Path for a .bed or .matrix file
    :param res: Resolution of the files
    :param chr: Which chromosome to extract information about
    :param bias: Bias dictionary
    :return: Numpy matrix of contact counts
    """
    d = {}
    with open(bed) as bf:
        for line in bf:
            l = line.strip().split('\t')
            if isChr(l[0], chr):
                d[int(l[3])] = toNearest(int(l[2]), res) // res
    n = max(d.values()) + 1
    o = np.zeros((n, n))
    with open(mat) as mf:
        for line in mf:
            l = line.strip().split('\t')
            if int(l[0]) in d and int(l[1]) in d:
                val = float(l[2])
                if bias:
                    val /= bias.get(int(l[1]), 1)
                    val /= bias.get(int(l[0]), 1)
                o[d[int(l[0])], d[int(l[1])]] = val
                o[d[int(l[1])], d[int(l[0])]] = val
    np.nan_to_num(o, copy=False, nan=0, posinf=0, neginf=0)
    return o


def parseBP(s):
    """
    :param s: string
    :return: string converted to number, taking account for kb or mb
    """
    if not s:
        return False
    if s.isnumeric():
        return int(s)
    s = s.lower()
    if "kb" in s:
        n = s.split("kb")[0]
        if not n.isnumeric():
            return False
        return int(n) * 1000
    elif "mb" in s:
        n = s.split("mb")[0]
        if not n.isnumeric():
            return False
        return int(n) * 1000000
    return False


def main():
    args = parse_args(sys.argv[1:])
    print("\n")

    f1 = args.f1_path
    f2 = args.f2_path
    if args.bed1 and args.mat1:
        f1 = args.mat1
    if args.bed2 and args.mat2:
        f2 = args.mat2

    if not (os.path.exists(f1) and os.path.exists(f2)):
        print("Error: Couldn't find specified contact files")
        return
    res = parseBP(args.resolution)
    if not res:
        print("Error: Invalid resolution")
        return
    distFilter = 0

    biasf1 = False
    if args.biasfile1:
        if os.path.exists(args.biasfile1):
            biasf1 = args.biasfile1
        else:
            print("Error: Couldn't find specified bias1 file")
            return

    biasf2 = False
    if args.biasfile2:
        if os.path.exists(args.biasfile2):
            biasf2 = args.biasfile2
        else:
            print("Error: Couldn't find specified bias2 file")
            return

    if args.distFilter:
        distFilter = parseBP(args.distFilter)
        if not distFilter:
            print("Error: Invalid distance filter")
            return

    tsvout = False
    try:
        if args.tsvout != 0:
            tsvout = float(args.tsvout)
    except:
        pass

    if not args.chromosome2 or args.chromosome2 == 'n':
        args.chromosome2 = args.chromosome

    o, x, y, q_val, changes_array = DCI(f1,
                           f2,
                           bed1=args.bed1,
                           bed2=args.bed2,
                           res=res,
                           sigma0=args.s_z,
                           s=args.s,
                           verbose=args.verbose,
                           mutual_nz=args.mutual_nz,                           
                           distance_filter=distFilter,
                           plot_results=args.plot,
                           bias1=biasf1,
                           bias2=biasf2,
                           changes=args.changedir,
                           low_memory=args.lowmemory,
                           chromosome=args.chromosome,
                           chromosome2=args.chromosome2,
                           norm_method = args.norm_method,
                           st=args.st,
                           tsvout=tsvout)
    

    if tsvout:

        if args.st:
        #indices = np.argwhere(o < tsvout)
        
            with open(args.outdir, 'w') as outfile:
                outfile.write('CHR1\tLOC1_start\tLOC1_end\tCHR2\tLOC2_start\tLOC2_end\tQ_VAL\tLOG_FOLD_CHANGE\n')
                for i in range(len(x)):
                    _x, _y = x[i], y[i]
                    outfile.write(
                        f'{args.chromosome}\t{_x * res}\t{_x * res + res}\t{args.chromosome2}\t{_y * res}\t{_y * res + res}\t{q_val[i]}\t{changes_array[i]}\n')
        else:
            indices = np.argwhere(o < tsvout)
            with open(args.outdir, 'w') as outfile:
                outfile.write('CHR1\tLOC1_start\tLOC1_end\tCHR2\tLOC2_start\tLOC2_end\tQ_VAL\tLOG_FOLD_CHANGE\n')
                for i in indices:
                    _x, _y = i[0], i[1]
                    outfile.write(
                        f'{args.chromosome}\t{_x * res}\t{_x * res + res}\t{args.chromosome2}\t{_y * res}\t{_y * res + res}\t{o[_x,_y]}\t{changes_array[_x,_y]}\n')
    else:
        np.save(args.outdir, o)


def sizeof_fmt(num, suffix='B'):
    ''' By Fred Cirera, after https://stackoverflow.com/a/1094933/1870254'''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


if __name__ == '__main__':
    main()
