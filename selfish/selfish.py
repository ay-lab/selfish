#!/usr/bin/env python3
import argparse
import os
import pathlib
import sys
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.stats.multitest as smm
from scipy.ndimage import gaussian_filter
from scipy.stats import norm

dir = os.path.dirname(__file__)
version_py = os.path.join(dir, "_version.py")
exec(open(version_py).read())


def parse_args(args):
    parser = argparse.ArgumentParser(description="Check the help flag")

    parser.add_argument("-f1", "--file1", dest="f1_path",
                        help="REQUIRED: Contact map 1", required=True)

    parser.add_argument("-f2", "--file2", dest="f2_path",
                        help="REQUIRED: Contact map 2", required=True)

    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="REQUIRED: where the output files\
                      will be written", required=True)

    parser.add_argument("-r", "--resolution", dest="resolution",
                        help="REQUIRED: Resolution used for the contact maps"
                        , required=True)

    parser.add_argument("-b", "--biases", dest="biasfile", \
                        help="RECOMMENDED: biases calculated by\
                        ICE or KR norm for each locus are read from BIASFILE", \
                        required=False)

    parser.add_argument("-sz", "--sigmaZero", dest="s_z", type=float, default=1.6,
                        help="OPTIONAL: sigma0 value for the method. DEFAULT is 1.6. Experimentally chosen for 5Kb resolution",
                        required=False)

    parser.add_argument("-i", "--iterations", dest="s", default=10,
                        type=int,
                        help="OPTIONAL: iteration count for the method. DEFAULT is 10. Experimentally chosen for 5Kb resolution", \
                        required=False)

    parser.add_argument("-d", "--distanceFilter", dest="distFilter",
                        help="OPTIONAL: If the data is too sparse for distant locations(ie. distance > 5Mb), you can filter them. DEFAULT is None",
                        required=False)

    parser.add_argument("-v", "--verbose", dest="verbose", type=bool,
                        default=True, help="OPTIONAL: Verbosity of the program",
                        required=False)

    parser.add_argument("-p", "--plot", dest="plot",
                        help="OPTIONAL: use this flag for generating plots. \
                      DEFAULT is False.", default=False, required=False)
    parser.add_argument("-lm", "--lowmem", dest="lowmemory",
                        help="OPTIONAL: Uses Float32 instead of Float64, halves the memory usage. Default is False",
                        default=False, required=False)
    parser.add_argument("-ch", "--chromosome", dest="chromosome",
                        help="REQUIRED for BED/Matrix input. Specify which chromosome to run the program for.",
                        default=0, required=False)

    parser.add_argument("-V", "--version", action="version", version="Selfish {}".format(__version__) \
                        , help="Print version and exit")

    return parser.parse_args()


def read_map_pd(path, res, bias, dt):
    sep = get_sep(path)
    o = pd.read_csv(path, sep=sep, header=None)
    o.dropna(inplace=True)
    if np.all(np.mod(o[0], res) == 0):
        o[0] //= res
        o[1] //= res
    n = max(max(o[0]), max(o[1])) + 1
    out = np.zeros((n, n), dtype=dt)
    out[o[0], o[1]] = o[2]
    out[o[1], o[0]] = o[2]
    if bias:
        factors = np.vectorize(bias.get)(o[0], 1)
        o[2] = np.multiply(o[2], factors)
        factors = np.vectorize(bias.get)(o[1], 1)
        o[2] = np.multiply(o[2], factors)

    return np.triu(out)


def get_diags(map):
    """
    Calculates the mean and standard deviation for every
    diagonal of a contact map
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


def normalize_map(map, non_zero):
    """
    :param map: contact map
    :param non_zero: same size matrix with map that has True as value for indices to be normalized.
    :return: Changes matrix in place to have diagonals = 0
    """
    means, stds = get_diags(map)
    x, y = np.nonzero(non_zero)
    distance = np.abs(x - y)
    m = np.vectorize(means.get)(distance)
    s = np.vectorize(stds.get)(distance)
    map[non_zero] -= m
    map[non_zero] /= s


def get_sep(f):
    with open(f) as file:
        for line in file:
            if "\t" in line:
                return '\t'
            elif " " in line.strip():
                return ' '
            elif "," in line:
                return ','
            break
    raise FileNotFoundError


def read_bias(f):
    d = defaultdict(lambda: 1.0)
    if f:
        sep = get_sep(f)
        with open(f) as file:
            for line in file:
                line = line.strip().split(sep)
                val = float(line[1])
                if not np.isnan(val):
                    d[int(float(line[0]))] = val
        return d
    else:
        return False


def DCI(f1,
        f2,
        res=100000,
        sigma0=1.6,
        s=10,
        plot_results=False,
        verbose=True,
        distance_filter=5000000,
        bias=False,
        chromosome=0,
        low_memory=False):
    scales = [(2 * (2 * (sigma0 * (2 ** (i / s)))) + 1) for i in range(s + 2)]

    dt = np.float32 if low_memory else np.float64

    biasDict = read_bias(bias)

    if type(f1) == str and type(f2) == str:
        if verbose: print("Reading Contact Map 1")
        if "bed" not in f1 and "matrix" not in f1:
            a = read_map_pd(f1, res, biasDict, dt)
        else:
            a = readBEDMAT(f1, res, chromosome, bias)

        if verbose: print("Reading Contact Map 2")
        if "bed" not in f2 and "matrix" not in f2:
            b = read_map_pd(f2, res, biasDict, dt)
        else:
            b = readBEDMAT(f2, res, chromosome, bias)

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
        print("Error: inputs should either be file names or square numpy matrices")
        return
    if verbose: print("Applying distance filter")
    if distance_filter > 0:
        a = np.tril(a, distance_filter // res)
        b = np.tril(b, distance_filter // res)

    non_zero_indices = np.logical_and(a != 0, b != 0)

    if plot_results:
        plt.clf()
        p_f = a.copy()
        p_f[p_f == 0] = 1
        sns.heatmap(np.abs(np.log(p_f)))
        plt.title("Contact counts (KR normalised, log scale) " + f1)
        plt.savefig(f1 + "_" + f2 + '_1_kr.png')
        plt.clf()
        p_f = b.copy()
        p_f[p_f == 0] = 1
        sns.heatmap(np.log(p_f))
        plt.title("Contact counts (KR normalised, log scale) " + f2)
        plt.savefig(f1 + "_" + f2 + '_2_kr.png')

    # np.save("kr_a.dat",a)
    # np.save("kr_b.dat",b)

    if verbose: print("Diagonal Normalizing Map 1")
    normalize_map(a, a != 0)

    if verbose: print("Diagonal Normalizing Map 2")
    normalize_map(b, b != 0)
    if plot_results:
        plt.clf()
        sns.heatmap(np.abs(a))
        plt.title("Contact counts (Diagonal normalised) " + f1)
        plt.savefig(f1 + "_" + f2 + '_1_diag.png')
        plt.clf()
        sns.heatmap(np.abs(b))
        plt.title("Contact counts (Diagonal normalised) " + f2)
        plt.savefig(f1 + "_" + f2 + '_2_diag.png')

        plt.clf()
        sns.heatmap(np.abs(np.abs(a) - np.abs(b)))
        plt.title("Contact difference (Diagonal normalised")
        plt.savefig(f1 + "_" + f2 + '_diff_diag.png')

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
        p_vals = norm.cdf(d_diff[non_zero_indices], loc=params[0], scale=params[1])
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
        sns.heatmap(np.abs(np.log2(o)))
        plt.title("Differential analysis")
        plt.savefig(f1 + "_" + f2 + "_selfish.png")
    return a,b,o


def sorted_indices(dci_out):
    ind = np.unravel_index(np.argsort(dci_out, axis=None), dci_out.shape)
    return ind[0], ind[1]

def toNearest(n,res):
    if n % res == 0:
        return n
    return (n//res + 1) * res
def isChr(s, c):
    return str(c) in re.findall("[1-9][0-9]*", s)

def readBEDMAT(f,res, chr, bias):
    bed = f if ".bed" in f else f.replace("matrix", "bed")
    mat = f if ".matrix" in f else f.replace("bed", "matrix")
    if ".bed" != bed[-4:] or ".matrix" != mat[-7:]:
        print("Error: Couldn't find matrix-bed combinations. They must have same names, only difference being extension. (.bed, .matrix)")
        raise FileNotFoundError
    d = {}
    with open(bed) as bf:
        for line in bf:
            l = line.strip().split('\t')
            if isChr(l[0], chr):
                d[int(l[3])] = toNearest(int(l[2]), res) // res
    n = max(d.values()) + 1
    o = np.zeros((n,n))
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
    if not (os.path.exists(f1) and os.path.exists(f2)):
        print("Error: Couldn't find specified contact files")
        return
    res = parseBP(args.resolution)
    if not res:
        print("Error: Invalid resolution")
        return
    distFilter = 0
    biasf = False
    if args.biasfile:
        if os.path.exists(args.biasfile):
            biasf = args.biasfile
        else:
            print("Error: Couldn't find specified bias file")
            return

    if args.distFilter:
        distFilter = parseBP(args.distFilter)
        if not distFilter:
            print("Error: Invalid distance filter")
            return

    o = DCI(f1, f2, res=res, sigma0=args.s_z,
            s=args.s, verbose=args.verbose,
            distance_filter=distFilter,
            plot_results=args.plot,
            bias=biasf,
            low_memory=args.lowmemory,
            chromosome=args.chromosome)

    np.save(str(pathlib.Path(args.outdir).joinpath("selfish.npy")), o)


def sizeof_fmt(num, suffix='B'):
    ''' By Fred Cirera, after https://stackoverflow.com/a/1094933/1870254'''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


if __name__ == '__main__':
    main()
