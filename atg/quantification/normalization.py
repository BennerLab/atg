import pandas
import numpy as np
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf


DEFAULT_ASYMPTOTIC_DISPERSION = 0.01
DEFAULT_EXTRA_POISSON = 3.0


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

    # ignore zero values so that log can be taken
    nonzero_values = (obs != 0) & (ref != 0)
    obs_nonzero = obs[nonzero_values]
    ref_nonzero = ref[nonzero_values]

    # log ratio of expression accounting for library size
    lr = np.log2((obs_nonzero / obs_sum) / (ref_nonzero / ref_sum))
    # absolute expression
    ae = (np.log2(obs_nonzero / obs_sum) + np.log2(ref_nonzero / ref_sum)) / 2
    # estimated asymptotic variance
    v = (obs_sum - obs_nonzero) / obs_sum / obs_nonzero + (ref_sum - ref_nonzero) / ref_sum / ref_nonzero
    # just keep values that exceed the absolute expression threshold
    m = ae > a_cutoff
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
    log_count_df = np.log2((count_df + 0.5) / (lib_size + 1) * 1e6)

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


def _rle_size_factor(count_df):
    """

    :param count_df:
    :return:
    """

    log_count_df = count_df.apply(np.log)

    size_factor = (log_count_df
                   .subtract(log_count_df.mean(axis=1), axis=0)
                   .replace([np.inf, -np.inf], np.nan)
                   .median()
                   .apply(np.exp))

    return size_factor


def rle_normalization(count_df):
    """
    normalization by relative log expression, as implemented in DESeq2.

    :param count_df:
    :return:
    """

    size_factor = _rle_size_factor(count_df)
    return count_df.div(size_factor, axis=1)


def _dispersion_genewise_mle(normalized_count_df, min_dispersion=1e-8):
    gene_mean = normalized_count_df.mean(axis=1)

    uncorrected_dispersion = ((((normalized_count_df.subtract(gene_mean, axis=0)) ** 2)
                                                    .subtract(gene_mean, axis=0))
                                                    .div(gene_mean ** 2, axis=0)).sum(axis=1)
    corrected_dispersion = uncorrected_dispersion/(normalized_count_df.shape[1]-1)

    return corrected_dispersion.clip(lower=min_dispersion)


def _fit_vst_parameters(count_df, sample_size=1000, min_mean_count=5.0, min_dispersion=1e-8):
    """

    :param count_df:
    :return:
    """

    normalized_count_df = rle_normalization(count_df)

    # ignore genes below the count threshold
    normalized_count_df['baseMean'] = normalized_count_df.mean(axis=1)
    detected_count_df = (normalized_count_df.query('baseMean >= @min_mean_count')
                                            .sort_values('baseMean'))

    # build sub-dataframe from genes spanning a range of genewise average expression
    spaced_index = np.linspace(0, detected_count_df.shape[0]-1, num=sample_size, dtype=int)
    sample_count_df = detected_count_df.iloc[spaced_index, 0:-1]
    dispersion_mle = _dispersion_genewise_mle(sample_count_df, min_dispersion)

    # for parametric fit, ignore genes with very low dispersion
    dispersion_for_fit = dispersion_mle[dispersion_mle > min_dispersion*100]
    dispersion_fit_df = dispersion_for_fit.to_frame(name="dispersion").join(normalized_count_df.ix[:, 'baseMean'])

    formula = 'dispersion ~ I(1/baseMean)'
    model = smf.glm(formula=formula, data=dispersion_fit_df,
                    family=sm.families.Gamma(link=statsmodels.genmod.families.links.identity)).fit()

    return model.params


def vst_transformation(count_df, asymptotic_dispersion=None, extra_poisson=None):
    """
    apply a variance stabilizing transformation to a count matrix, as in DESeq2.

    :param count_df:
    :param asymptotic_dispersion: aka "a0"
    :param extra_poisson: aka "a1"
    :return:
    """

    normalized_count_df = rle_normalization(count_df)

    asymptotic_dispersion_fit, extra_poisson_fit = _fit_vst_parameters(count_df)

    if asymptotic_dispersion and extra_poisson:
        pass
    else:
        if not extra_poisson:
            extra_poisson = extra_poisson_fit
        if not asymptotic_dispersion:
            asymptotic_dispersion = asymptotic_dispersion_fit

    return (np.log((1 + extra_poisson + 2 * asymptotic_dispersion * normalized_count_df +
                    2 * np.sqrt(asymptotic_dispersion * normalized_count_df *
                        (1 + extra_poisson + asymptotic_dispersion * normalized_count_df)))
                    / (4 * asymptotic_dispersion))
            / np.log(2))


def normalize_counts(namespace):
    count_df = pandas.read_csv(namespace.count_filename, index_col=namespace.index, delimiter=namespace.delimiter)
    if namespace.method == 'tmm':
        normalized_df = tmm_normalization(count_df)

    elif namespace.method == 'cpm':
        normalized_df = cpm_normalization(count_df)

    elif namespace.method == 'voom':
        normalized_df = voom(count_df)

    elif namespace.method == 'rle':
        normalized_df = rle_normalization(count_df)

    elif namespace.method == 'vst':
        normalized_df = vst_transformation(count_df)

    else:
        normalized_df = count_df

    normalized_df.to_csv(namespace.output_filename, sep=namespace.delimiter)


def setup_subparsers(subparsers):
    normalization_parser = subparsers.add_parser('normalize', help='Normalize read counts')

    normalization_parser.add_argument('count_filename', help="delimited file containing read counts")
    normalization_parser.add_argument('output_filename', help="")
    normalization_parser.add_argument('-m', '--method', default="TMM", choices=['tmm', 'cpm', 'voom', 'rle', 'vst'],
                                      help="tmm: weighted Trimmed Mean of M-values; "
                                           "cpm: Counts Per Million reads"
                                           "voom: log2 count after TMM"
                                           "rle: Relative Log Expression"
                                           "vst: Variance Stabilizing Transformation")
    normalization_parser.add_argument('-i', '--index', default=0, type=int,
                                      help="numerical position of the ID column")
    normalization_parser.add_argument('-d', '--delimiter', default=",", help='Text delimiter, e.g. ","')

    normalization_parser.set_defaults(func=normalize_counts)
