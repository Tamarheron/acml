#!/usr/bin/env python

import pandas as pd
#  import numpy as np
import statsmodels.formula.api as smf
from matplotlib import pyplot as plt
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import warnings
import matplotlib.patches as mpatches
pandas2ri.activate()

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    r("""
    library(miscTools)
    library(maxLik)
    library(truncreg)
    """)


def coefplot(formula, data, fontsize=5):
    """ Plots coefficients of a regression model.
        formula = patsy-style formula for regression model
        data = pandas dataframe with columns for the variables in the formula
        fontsize = 5 by default
    """
    lm = smf.ols(formula, data=data).fit()
    lm0 = smf.ols(formula + "+ 0", data=data).fit()
    r.assign("data", data)
    r("""
    trunc_reg <- truncreg(%s,
                          data = data,
                          point = 0,
                          direction = 'left')
    summ <- summary(trunc_reg)
    coeffs <- trunc_reg$coefficients
    coeffs <- coeffs[names(coeffs) != "sigma"]
    coeffs_values <- as.vector(coeffs, mode="numeric")
    """ % (formula))
    params = pd.DataFrame(data=r("coeffs_values"), index=r("names(coeffs)"))
    params.index = [":".join(sorted(name.split(":"))) for name in params.index]
    params = params.drop("(Intercept)")
    params.columns = ["truncreg"]
    lm_params = lm.params.drop("Intercept")
    lm_params.index = [":".join(sorted(name.split(":"))) for name in
                       lm_params.index]
    params["lm"] = lm_params
    lm0_params = lm0.params
    lm0_params.index = [":".join(sorted(name.split(":"))) for name in
                        lm0_params.index]
    params["lm0"] = lm0_params
    params = params.sort_values("lm")
    lm_bse = lm.bse
    lm_bse.index = [":".join(sorted(name.split(":"))) for name in
                    lm_bse.index]
    lm0_bse = lm0.bse
    lm0_bse.index = [":".join(sorted(name.split(":"))) for name in
                     lm0_bse.index]
    fig, ax = plt.subplots()
    y = range(len(params.index))
    ax.scatter(list(params["truncreg"]), y, color="g", s=2)
    ax.scatter(list(params["lm"]), y, color="r", s=2)
    ax.scatter(list(params["lm0"]), y,  color="b", s=2)
    for y in range(len(params.index)):
        sub = params.index[y]
        x = params.lm[sub]
        se = lm_bse[sub]
        ax.plot([x - se, x + se], [y, y], color="red")
        x = params.lm0[sub]
        se = lm0_bse[sub]
        ax.plot([x - se, x + se], [y, y], color="blue")
    red_patch = mpatches.Patch(color='red', label='Linear Model')
    blue_patch = mpatches.Patch(color='blue', label='Forced Zero Intercept')
    green_patch = mpatches.Patch(color='green', label='Truncated Regression')
    plt.legend(handles=[red_patch, blue_patch, green_patch], loc=2)
    plt.yticks(range(len(params.index)), params.index)
    ax.set_ylim([-1, len(params)])
    ax.set_yticklabels(["\n".join(name.split(":")) for name in params.index],
                       fontsize=fontsize)
    ax.set_ylabel("Substitutions")
    ax.set_xlabel("Coefficients")
    plt.title("Coefficient plot")
    plt.grid()
    fig.savefig("coefplot.png", dpi=200)
    file = open("lm_summary.txt", "w")
    file.write(str(lm.summary()))
    file.close()
    file = open("lm0_summary.txt", "w")
    file.write(str(lm0.summary()))
    file.close()
    file = open("truncreg_summary.txt", "w")
    file.write(str(r("summ")))
    file.close()


if __name__ == "__main__":
    data = pd.read_csv("culled.csv", sep=" ")
    data = data.drop("Unnamed: 75", axis=1)
    file = open("predictors.csv", "r")
    predictors = file.read().split("\n")[:20]
    file.close()
    coefplot("AGDIST ~ %s" % ("+".join(predictors)), data, fontsize=5)
