#!/usr/bin/env python

#*****************************************************************************80
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 January 2018
#
#  Author:
#
#    Su Yi (suyi@hust.edu.cn)

from nilaplace import numerical_inverse_laplace_transform
from scipy.special import gamma, hyp2f1
import numpy as np
from stirling import stirling2
import math

no = 4

def exp_pdf_laplace(s, mean, lambd=None):
    laplace_exp = 1.0/(1+mean*s)
    return laplace_exp


def exp_mean(mean, lambd=None):
    return mean


def exp_var(mean, lambd=None):
    return mean**2


def determin_pdf_laplace(s, mean, lambd=None):
    laplace_deter = np.exp(-mean*s)
    return laplace_deter


def determin_mean(mean, lambd=None):
    return mean


def determin_var(mean, lambd=None):
    return 0.0


def mm1k_sojourn_pdf_laplace(s, b, lambd):
    rho = lambd * b
    mu = 1 / b
    P0 = 1 / (sum([rho**k for k in range(0, no+1)]))
    # Pn = ((1-rho) * rho**n) / (1 - rho**(n+1))
    Pn = rho**no / (sum([rho**k for k in range(0, no+1)]))
    laplace_disk_lat = ((mu*P0)/(1-Pn)) * ((1 - (lambd/(mu+s))**no) / (mu - lambd + s))
    return laplace_disk_lat


def mm1k_mean(b, lambd):
    rho = lambd * b
    mu = 1 / b
    Nbar = (rho * (1 - (no+1)*rho**no + no*rho**(no+1))) / ((1 - rho) * (1 - rho**(no+1)))
    Pn = rho**no / (sum([rho**k for k in range(0, no+1)]))
    Bbar = Nbar / (lambd*(1 - Pn))
    return Bbar


def mm1_sojourn_pdf_laplace(s, b, lambd):
    rho = lambd * b
    mu = 1.0 / b
    pass


def mm1_mean(b, lambd):
    rho = lambd * b
    mu = 1.0 / b
    return 1.0 / (mu*(1-rho))


def get_cache_pdf_laplace(lambd, theta, mean1, pdf_laplace_func1, mean2, pdf_laplace_func2):
    # theta is the miss ratio
    lambd1 = lambd*(1-theta)
    lambd2 = lambd*theta

    def cache_pdf_laplace(s, mean, lambd):
        laplace_pdf = (1-theta)*pdf_laplace_func1(s, mean1, lambd1) + theta*pdf_laplace_func2(s, mean2, lambd2)
        return laplace_pdf

    return cache_pdf_laplace

def cache_cdf(t, lambd, theta, mean1, pdf_laplace_func1, mean2, pdf_laplace_func2):
    # theta is the miss ratio
    lambd1 = lambd*(1-theta)
    lambd2 = lambd*theta

    def cache_cdf_laplace(s):
        laplace_pdf = (1-theta)*pdf_laplace_func1(s, mean1, lambd1) + theta*pdf_laplace_func2(s, mean2, lambd2)
        return laplace_pdf/s

    p = numerical_inverse_laplace_transform(cache_cdf_laplace, t, n=6)

    return p


def cache_mean(lambd, theta, mean1, mean_func1, mean2, mean_func2):
    return (1-theta)*mean_func1(mean1, lambd*(1-theta)) + theta*mean_func2(mean2, lambd*theta)


def mg1_wait_pdf_laplace(s, mean, lambd, laplace_pdf_func):
    rho = lambd * mean
    laplace_wait = (1-rho) * s / (lambd * laplace_pdf_func(s, mean, lambd) + s - lambd)
    return laplace_wait


def mg1_sojourn_pdf_laplace(s, mean, lambd, laplace_pdf_func):
    laplace_sojourn = mg1_wait_pdf_laplace(s, mean, lambd, laplace_pdf_func) * laplace_pdf_func(s, mean, lambd)
    return laplace_sojourn


def mg1_sojourn_cdf_laplace(s, mean, lambd, laplace_pdf_func):
    return mg1_sojourn_pdf_laplace(s, mean, lambd, laplace_pdf_func) / s


def mg1_sojourn_cdf(t, mean, lambd, laplace_pdf_func):

    def cdf_laplace(s):
        return mg1_sojourn_cdf_laplace(s, mean, lambd, laplace_pdf_func)

    p = numerical_inverse_laplace_transform(cdf_laplace, t, n=6)
    return p


def mg1_queuelen_pgf(z, mean, lambd, laplace_pdf_func):
    # This is the Probability-generating function (pgf), not the Probability mass function (pmf).
    # The pmf is corresponding to pdf of continuous random variable for discrete random variable.
    # pgf(z): G(z) = \sum_{k=0}^{\inf}(p(k)*z^k), p(k) is the pmf(k).
    # Hence:
    # pgf(0) = pmf(0), because only 0^0=1
    # pgf(1) = \sum_{k=0}^{inf}(pmf(k)) = 1
    rho = lambd * mean
    tmp = laplace_pdf_func(lambd-lambd*z, mean, lambd)
    if tmp == z:
        print "zero denominator:", z
        laplace_queuelen = 0
    else:
        laplace_queuelen = tmp * (1-rho) * (1-z)/(tmp - z)
    return laplace_queuelen


def mg1_queuelen_pmf(k, mean, variance, lambd, laplace_pdf_func):
    if k == 0:
        queuelen_p = mg1_queuelen_pgf(k, mean, lambd, laplace_pdf_func)
    else:
        cv = variance**0.5 / mean
        rho = lambd * mean
        q = rho * (cv**0.5 + 1) / (2 + rho*(cv**0.5 - 1))
        queuelen_p = rho * (q**(k-1)) * (1-q)
    # print queuelen_p, "#", k, mean, variance, lambd, "#"
    return queuelen_p


def mm1_queuelen_pmf(k, mean, lambd):
    rho = mean * lambd
    return (1-rho)*(rho**k)


def get_queuelen_pmf(mean, variance, lambd, laplace_pdf_func):
    
    def pmf_func(k):
        if mean**2 == variance:
            # print "M/M/1"
            # p_k = mm1_queuelen_pmf(k, mean, lambd)
            p_k = mg1_queuelen_pmf(k, mean, variance, lambd, laplace_pdf_func)
        else:
            p_k = mg1_queuelen_pmf(k, mean, variance, lambd, laplace_pdf_func)
        # print k, p_k
        return p_k

    return pmf_func


def n_not_blocked(n, m, theta):
    m = float(m)
    n_not_blocked_p = 1 - ((math.factorial(m)*stirling2(n, m)[n-1, m-1])/(m**n))
    return n_not_blocked_p


def get_probability_hit_not_blocked(m, theta, pmf_func):
    hit_not_blocked_p = 0.0
    for k in range(0, m):
        hit_not_blocked_p += (1-theta) * pmf_func(k)
        # print k, hit_not_blocked_p
    k = m
    while pmf_func(k) > 0.00001:
        hit_not_blocked_p += (1-theta) * pmf_func(k) * n_not_blocked(k, m, theta)
        # print k, hit_not_blocked_p
        k += 1
    return hit_not_blocked_p


def cosmodel_reqbased_backend(t, lambd, theta, no, mean1, var1, laplace_pdf_func1, mean2, var2, laplace_pdf_func2):
    lambd_m = lambd * theta
    rho_dev = lambd_m * mean2
    pmf_func = get_queuelen_pmf(mean2, var2, lambd_m, laplace_pdf_func2)
    p_star = get_probability_hit_not_blocked(no, theta, pmf_func)
    # print "p_star :", p_star

    lambd_ag = lambd * (1 - p_star)
    theta_ag = theta / (1 - p_star)
    mean_ag = theta_ag * mean2 + (1 - theta_ag) * mean1
    rho_ag = mean_ag * lambd_ag

    # print "lambd_ag :", lambd_ag, "theta_ag :", theta_ag, "mean_ag :", mean_ag, "rho_ag :", rho_ag
    # print "lambd_m :", lambd_m, "lambd_ag - lambd_m", lambd_ag - lambd_m

    def ag_pdf_laplace(s, mean_param, lambd_param):
        laplace_p = theta_ag * laplace_pdf_func2(s, mean2, lambd_m) + (1 - theta_ag) * laplace_pdf_func1(s, mean1, lambd_ag - lambd_m)
        return laplace_p

    p_ag = (1 - p_star) * mg1_sojourn_cdf(t, mean_ag, lambd_ag, ag_pdf_laplace)

    return (p_star, p_ag)


if __name__ == '__main__':
    # mean = 5.0
    # variance = 0
    # lambd = 1.0/mean * 0.90

    # p = 0.0
    # k = 0
    # while k < 10:
    #     pp = mg1_queuelen_pmf(k, mean, variance, lambd, determin_pdf_laplace)
    #     p += pp
    #     print k, pp, p
    #     k += 1
    no = 10
    srate = 100.0
    cache_hit_ratio = 0.65
    arate_max = srate/100.0/(1-cache_hit_ratio)
    arate = arate_max*0.7
    theta = 1 - 0.675
    lambd = arate
    mean1 = 1.0/srate
    mean2 = 1.0/(srate/100.0)

    pdf_laplace_func1 = exp_pdf_laplace
    mean_func1 = exp_mean
    var1 = exp_var(mean1)
    pdf_laplace_func2 = exp_pdf_laplace
    var2 = exp_var(mean2)
    lats = [0.0035, 0.009, 0.023, 1.422, 3.436, 5.491, 10.373]

    for t in lats:
        print t, cosmodel_reqbased_backend(t, lambd, theta, no, mean1, var1, pdf_laplace_func1, mean2, var2, pdf_laplace_func2)