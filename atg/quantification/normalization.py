import pandas
import numpy as np

# TODO: implement DESeq normalization (relative log expression)


def cpm_normalization(count_df):
    normalized_df = count_df.div(count_df.sum()) * 10**6
    return normalized_df


def _tmm_norm_factor_single(obs, ref, log_ratio_trim=0.3, sum_trim=0.05, weighting=True, a_cutoff=-1e10):
    if all(obs == ref):
        return 1

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


def normalize_counts(namespace):
    count_df = pandas.read_csv(namespace.count_filename, index_col=namespace.index, delimiter=namespace.delimiter)
    if namespace.method == 'TMM':
        normalized_df = tmm_normalization(count_df)

    elif namespace.method == 'CPM':
        normalized_df = cpm_normalization(count_df)

    normalized_df.to_csv(namespace.output_filename, sep=namespace.delimiter)


def setup_subparsers(subparsers):
    normalization_parser = subparsers.add_parser('normalize', help='Normalize read counts')

    normalization_parser.add_argument('count_filename', help="delimited file containing read counts")
    normalization_parser.add_argument('output_filename', help="")
    normalization_parser.add_argument('-m', '--method', default="TMM", choices=['TMM', 'CPM'],
                                      help="TMM: weighted Trimmed Mean of M-values; "
                                           "CPM: Counts Per Million reads")
    normalization_parser.add_argument('-i', '--index', default=0, type=int,
                                      help="numerical position of the ID column")
    normalization_parser.add_argument('-d', '--delimiter', default=",", help='Text delimiter, e.g. ","')

    normalization_parser.set_defaults(func=normalize_counts)


# tmm normalized counts = count / (library size * normalization factor)
# In all the downstream code, the lib.size and norm.factors are
# multiplied together to act as the effective library size; this
# (product) would be similar to DESeq's size factor.