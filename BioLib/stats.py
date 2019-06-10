#!/usr/bin/env python

'''
for python3
'''

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

import sys
import os
import math
from copy import deepcopy


def percentile(N=None, percent=None, key=lambda x:x):
    '''
    Find the percentile of a list of values.
    Taken from http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    '''
    if not N:
        return None

    N = sorted(N)

    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)

    return d0+d1


def mean(data=None):
    """Return the sample arithmetic mean of data.
    http://stackoverflow.com/a/27758326/632242
    """
    n = len(data)
    #if n < 1:
    #    raise ValueError('mean requires at least one data point')
    return sum(data)/float(n)


def _ss(data=None):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def pstdev(data=None):
    """Calculates the population standard deviation."""
    n = len(data)
    #if n < 2:
    #    raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5


def median(N=None, key=lambda x:x):
    """median is 50th percentile."""
    return percentile(N, 0.5, key)


def multiple_testing_correction(pvalues=None, correction_type="Benjamini-Hochberg"):
    """
    Copyright 2017 Francisco Pina Martins <f.pinamartins@gmail.com>
    Taken from https://stackoverflow.com/a/21739593/3091595, remove numpy dependence

    @parameter pvalues - a list of pvalues
    @parameter correction_type - pvalue correction method

    @return qvalues - a list of qvalues
    """
    n = len(pvalues)
    qvalues = [0]*n
    if correction_type == "Bonferroni":
        qvalues = n * pvalues

    elif correction_type == "Bonferroni-Holm":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            qvalues[i] = (n-rank) * pvalue

    elif correction_type == "Benjamini-Hochberg":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in xrange(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            qvalues[index] = new_values[i]

    return qvalues



