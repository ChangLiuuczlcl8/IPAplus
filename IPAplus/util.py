import numpy as np
import pandas as pd
from scipy.stats import betabinom


def cal_posterior_3(prior, model_prediction, p):
    bbn_distribution = betabinom.pmf(np.arange(11), 10, p, 1)
    model_prediction_rank = pd.DataFrame(model_prediction).rank(method='min', ascending=False) - 1
    model_prediction_rank = model_prediction_rank[0]
    model_prediction_rank.loc[model_prediction_rank == model_prediction_rank.max()] = 10
    likelihood_3 = bbn_distribution[model_prediction_rank.astype(int)]
    posterior = likelihood_3 * prior
    return posterior / sum(posterior)


def split_every_n_special(s, special="_", n=2):
    parts = s.split(special)  # split by the special character
    result = []
    for i in range(0, len(parts), n):
        group = special.join(parts[i:i+n])
        result.append(group)
    return result
