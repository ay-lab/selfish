#!/usr/bin/env python3
import argparse
import os
import sys
import re
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
        "-ch",
        "--chromosome",
        dest="chromosome",
        help="REQUIRED: Specify which chromosome to run the program for.",
        default='n',
        required=True)

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

    parser.add_argument("-V", "--version", action="version", version="Selfish {}".format(__version__) \
                        , help="Print version and exit")

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
        o[4] = np.multiply(o[4], factors)
        factors = np.vectorize(bias.get)(o[3], 1)
        o[4] = np.multiply(o[4], factors)

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
            means[i] = mean
            stds[i] = std
        else:
            means[i] = 0
            stds[i] = 1
    return means, stds


def normalize_map(cmap, non_zero):
    """
    :param cmap: contact map
    :param non_zero: same size matrix with map that has True as value for indices to be normalized.
    :return: Changes matrix in place to have each diagonal sum = 0
    """
    means, stds = get_diags(cmap)
    x, y = np.nonzero(non_zero)
    distance = np.abs(x - y)
    m = np.vectorize(means.get)(distance)
    s = np.vectorize(stds.get)(distance)
    cmap[non_zero] -= m
    cmap[non_zero] /= s


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
        distance_filter=5000000,
        bias1=False,
        bias2=False,
        chromosome='n',
        low_memory=False):
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
    if plot_results:
        import seaborn as sns
    scales = [(sigma0 * (2**(i / s))) for i in range(1, s + 3)]

    dt = np.float32 if low_memory else np.float64

    biasDict1 = read_bias(bias1, chromosome, res)
    biasDict2 = read_bias(bias2, chromosome, res)

    if type(f1) == str and type(f2) == str:
        if verbose: print("Reading Contact Map 1")
        if bed1 != '' and f1 != '':
            if chromosome == 'n':
                print("You need to specify a chromosome for matrix files.")
                raise FileNotFoundError
            a = readBEDMAT(bed1, f1, res, chromosome, biasDict1)
        elif ".hic" in f1:
            if chromosome == 'n':
                print("You need to specify a chromosome for hic files.")
                raise FileNotFoundError
            a = readHiCFile(f1, chromosome, res)
        elif f1.endswith('.cool'): #cooler.fileops.is_cooler(f1):
            a = readCoolFile(f1, chromosome)
        elif f1.endswith('.mcool'): #cooler.fileops.is_multires_file(f1):
            a = readMultiCoolFile(f1, chromosome, res)
        else:
            a = read_map_pd(f1, res, biasDict1, dt, chromosome)

        if verbose: print("Reading Contact Map 2")
        if bed2 != '' and f2 != '':
            if chromosome == 'n':
                print("You need to specify a chromosome for matrix files.")
                raise FileNotFoundError
            b = readBEDMAT(bed2, f2, res, chromosome, biasDict2)
        elif ".hic" in f2:
            if chromosome == 'n':
                print("You need to specify a chromosome for hic files.")
                raise FileNotFoundError
            b = readHiCFile(f2, chromosome, res)
        elif f2.endswith('.cool'): #cooler.fileops.is_cooler(f1):
            b = readCoolFile(f2, chromosome)
        elif f2.endswith('.mcool'): #cooler.fileops.is_multires_file(f1):
            b = readMultiCoolFile(f2, chromosome, res)
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
    if verbose: print("Applying distance filter")
    if distance_filter > 0:
        a = np.tril(a, distance_filter // res)
        b = np.tril(b, distance_filter // res)

    tmp_n = max(max(a.shape), max(b.shape))
    tmp = np.zeros((tmp_n, tmp_n))
    tmp[:a.shape[0], :a.shape[1]] = a
    a = tmp.copy()
    tmp = np.zeros((tmp_n, tmp_n))
    tmp[:b.shape[0], :b.shape[1]] = b
    b = tmp.copy()
    tmp = None

    non_zero_indices = np.logical_or(a != 0, b != 0)

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
        sns.heatmap(np.log(p_f))
        plt.title("Contact counts (KR normalised, log scale) " + f2)
        plt.savefig('f2_2_kr.png')

    # np.save("kr_a.dat",a)
    # np.save("kr_b.dat",b)
    if changes != "":
        c_temp = np.zeros_like(a)
        c_temp[non_zero_indices] = np.log2(np.divide(a + 1, b + 1))[non_zero_indices]
        np.save(changes, c_temp)
        c_temp = None

    if verbose: print("Diagonal Normalizing Map 1")
    normalize_map(a, non_zero_indices)

    if verbose: print("Diagonal Normalizing Map 2")
    normalize_map(b, non_zero_indices)
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
        plt.title("Contact difference (Diagonal normalised")
        plt.savefig('f12_diff_diag.png')

    # np.save("diag_a.dat",a)
    # np.save("diag_b.dat",b)
    if verbose: print("Applying gaussians")

    final_p = np.ones(len(a[non_zero_indices]))



    size = a.shape[0]
    b -= a
    #a = None
    diff = b.copy()
    #b = None
    d_pre = gaussian_filter(diff, scales[0])
    count = 1
    for scale in scales[1:]:
        print("Gaussian:", count, "of", len(scales) - 1)
        d_post = gaussian_filter(diff, scale)
        d_diff = d_post - d_pre
        params = norm.fit(d_diff[non_zero_indices])
        p_vals = norm.cdf(d_diff[non_zero_indices],
                          loc=params[0],
                          scale=params[1])
        p_vals[p_vals > 0.5] = 1 - p_vals[p_vals > 0.5]
        p_vals *= 2
        final_p[p_vals < final_p] = p_vals[p_vals < final_p]
        d_pre = d_post.copy()
        count += 1

    d_diff = None
    o = np.ones(diff.shape, dtype=dt)
    diff = None
    _, out_p = smm.multipletests(final_p, method='fdr_bh')[:2]

    o[non_zero_indices] = out_p
    o[o == 0] = 1
    if plot_results:
        plt.clf()
        sns.heatmap(np.abs(np.log10(o)))
        plt.title("Differential analysis")
        plt.savefig("f1f2_selfish.png")
    return o


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


def readHiCFile(f, chr, res):
    """
    :param f: .hic file path
    :param chr: Which chromosome to read the file for
    :param res: Resolution to extract information from
    :return: Numpy matrix of contact counts
    """
    result = straw.straw('KR', f, str(chr), str(chr), 'BP', res)
    x = np.array(result[0]) // res
    y = np.array(result[1]) // res
    val = np.array(result[2])
    n = max(max(x), max(y))
    o = np.zeros((n, n)) + 1
    o[x, y] = val
    o[y, x] = val
    return o

def readCoolFile(f, chr):
    """
    :param f: .cool file path
    :param chr: Which chromosome to read the file for
    :return: Numpy matrix of contact counts
    """
    clr = cooler.Cooler(f)
    result = clr.matrix(balance=True).fetch(chr)
    result[np.isnan(result)] = 0
    return result

def readMultiCoolFile(f, chr, res):
    """
    :param f: .cool file path
    :param chr: Which chromosome to read the file for
    :param res: Resolution to extract information from
    :return: Numpy matrix of contact counts
    """
    uri = '%s::/resolutions/%s' % (f, res)
    clr = cooler.Cooler(uri)
    result = clr.matrix(balance=True).fetch(chr)
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
                    val *= bias.get(int(l[1]), 1)
                    val *= bias.get(int(l[0]), 1)
                o[d[int(l[0])], d[int(l[1])]] = val
                o[d[int(l[1])], d[int(l[0])]] = val
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

    o = DCI(f1,
            f2,
            bed1=args.bed1,
            bed2=args.bed2,
            res=res,
            sigma0=args.s_z,
            s=args.s,
            verbose=args.verbose,
            distance_filter=distFilter,
            plot_results=args.plot,
            bias1=biasf1,
            bias2=biasf2,
            changes=args.changedir,
            low_memory=args.lowmemory,
            chromosome=args.chromosome)

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
