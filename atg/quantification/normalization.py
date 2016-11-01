import pandas
import numpy as np


def cpm_normalization(count_df):
    normalized_df = count_df.div(count_df.sum()) * 10**6
    return normalized_df


def _tmm_norm_factor_single(obs, ref, log_ratio_trim=0.3, sum_trim=0.05, weighting=True, a_cutoff=-1e10):
    """

    :param obs:
    :param ref:
    :param log_ratio_trim:
    :param sum_trim:
    :param weighting:
    :param a_cutoff:
    :return:
    """

    obs_sum = obs.sum()
    ref_sum = ref.sum()

    # log ration of expression accounting for library size
    lr = np.log2((obs / obs_sum) / (ref / ref_sum))
    # absolute expression
    ae = (np.log2(obs / obs_sum) + np.log2(ref / ref_sum)) / 2
    # estimated asymptotic variance
    v = (obs_sum - obs) / obs_sum / obs + (ref_sum - ref) / ref_sum / ref
    # create mask
    m = np.isfinite(lr) & np.isfinite(ae) & (ae > a_cutoff)
    # drop the masked values
    lr = lr[m]
    ae = ae[m]
    v = v[m]
    assert len(lr) == len(ae) == len(v)

    n = len(lr)
    lo_l = np.floor(n * log_ratio_trim) + 1
    hi_l = n + 1 - lo_l
    lo_s = np.floor(n * sum_trim) + 1
    hi_s = n + 1 - lo_s
    k = ((lr.rank(method="first") >= lo_l) & (lr.rank(method="first") <= hi_l)) \
        & ((ae.rank(method="first") >= lo_s) & (ae.rank(method="first") <= hi_s))

    if weighting:
        return 2 ** (sum(lr[k] / v[k]) / sum(1 / v[k]))
    else:
        return 2 ** (lr[k].mean())


def tmm_norm_factor(count_df):
    """

    :param count_df:
    :return:
    """
    count_nonzero = (count_df.where(count_df > 0)
                     .dropna(how="all").fillna(0))

    # determine reference sample based on expression at 75th quantile
    lib_size = count_nonzero.sum()
    y = (count_nonzero
         .div(lib_size)
         .quantile(0.75))

    sf = y.div(np.exp(np.mean(np.log(y))))
    ref_index = (sf-sf.mean()).abs().idxmin()
    sf_tmm = count_nonzero.apply(_tmm_norm_factor_single, ref=count_nonzero.ix[:, ref_index])

    # adjust size factors to multiply to 1
    final_sf_tmm = sf_tmm / np.exp(np.mean(np.log(sf_tmm)))

    return final_sf_tmm


def tmm_normalization(count_df):
    norm_factor = tmm_norm_factor(count_df)
    normalized_df = count_df.div(count_df.sum() * norm_factor) * 10**6
    return normalized_df


def voom(count_df):
    """
    calculate log2 cpm on TMM-normalized data, as in voom. Further calculation of the dispersion trend would be used
    in linear modeling.

    :param count_df:
    :return: a log2, normalized count dataframe
    """

    lib_size = count_df.sum() * tmm_norm_factor(count_df)
    log_count_df = np.log2(count_df / (lib_size + 1) * 1e6)

    # further calculations for modeling dispersion trend
    # row-wise mean and SD
    # gene_mean = count_df.mean(axis=1)
    # gene_sd = count_df.std(axis=1)
    #
    # # fit lowess trend to sqrt-standard-deviations by log-count-size
    # sx = gene_mean + (np.log2(lib_size+1)-np.log2(1e6)).mean()
    # sy = np.sqrt(gene_sd)
    #
    # lowess_fit = sm.nonparametric.lowess(sy, sx, frac=0.5)

    return log_count_df


def rle_normalization(count_df):
    """
    normalization by relative log expression, as implemented in DESeq2.

    :param count_df:
    :return:
    """

    log_count_df = count_df.apply(np.log)

    size_factor = (log_count_df
                   .subtract(log_count_df.mean(axis=1), axis=0)
                   .replace([np.inf, -np.inf], np.nan)
                   .median()
                   .apply(np.exp))

    return count_df.div(size_factor, axis=1)


def normalize_counts(namespace):
    count_df = pandas.read_csv(namespace.count_filename, index_col=namespace.index, delimiter=namespace.delimiter)
    if namespace.method == 'TMM':
        normalized_df = tmm_normalization(count_df)

    elif namespace.method == 'CPM':
        normalized_df = cpm_normalization(count_df)

    elif namespace.method == 'voom':
        normalized_df = voom(count_df)

    elif namespace.method == 'RLE':
        normalized_df = rle_normalization(count_df)

    else:
        normalized_df = count_df

    normalized_df.to_csv(namespace.output_filename, sep=namespace.delimiter)


def setup_subparsers(subparsers):
    normalization_parser = subparsers.add_parser('normalize', help='Normalize read counts')

    normalization_parser.add_argument('count_filename', help="delimited file containing read counts")
    normalization_parser.add_argument('output_filename', help="")
    normalization_parser.add_argument('-m', '--method', default="TMM", choices=['TMM', 'CPM', 'voom'],
                                      help="TMM: weighted Trimmed Mean of M-values; "
                                           "CPM: Counts Per Million reads"
                                           "voom: log2 count after TMM"
                                           "RLE: Relative Log Expression")
    normalization_parser.add_argument('-i', '--index', default=0, type=int,
                                      help="numerical position of the ID column")
    normalization_parser.add_argument('-d', '--delimiter', default=",", help='Text delimiter, e.g. ","')

    normalization_parser.set_defaults(func=normalize_counts)


# tmm normalized counts = count / (library size * normalization factor)
# In all the downstream code, the lib.size and norm.factors are
# multiplied together to act as the effective library size; this
# (product) would be similar to DESeq's size factor.
