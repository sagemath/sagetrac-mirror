# -*- coding: utf-8 -*-
r"""
Time series network (correlation graph)

This module defines an object class for creating a graph with edges defined by 
correlations between elements of a set of fixed-length vectors.  While the 
motivating context is time series, any real-valued vectors will work.

{INDEX_OF_METHODS}

## Design notes

### Inputs
** `x` A 2-D array containing multiple variables and observations. Each row 
    of `x` contains a variable, and each column a single observation of those 
    variables
** `labels` A list of descriptive labels for the variables; must be of the same
    length as the number of rows in `x`. These labels will be shown on the 
    vertices of the graph
** `alpha=0.1` The significance level the algorithm should use to determine
    which vertices should be connected by edges. Only pairs of vertices whose
    observations are correlated at a significance level of `alpha` will be
    connected

### Constraints
There are currently no constraints on the size of `x` or the values it contains.

### Implementation
The class's `__init__` method computes a correlation matrix from `x` using
the `numpy.corrcoef` method. It then generates the graph using the specified
labels and `alpha` parameter.

The `plot` method updates the graph (to incorporate any changes in `alpha`) 
and then generates a plot using the Sage `Graph` plot method. If there are 
more than two vertices, the plot is generated using a circular layout; however
with two vertices edges are not displayed in circular layout, so no layout 
is specified.

The correlation matrix obtained from `x` is stored as a property called `cm`. 
It is a 2-D numpy array.

AUTHORS:

-  Robert J. Kingan (2018-07-26): initial version

"""

# ****************************************************************************
#       Copyright (C) 2018 Robert Kingan (rjk2008@yahoo.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import numpy as np 
from scipy.stats import t 
from graph import Graph

def _is_sig(r, n, alpha):
    tscore = r * np.sqrt(n-2) / np.sqrt(1-r*r)
    tp = t.pdf(tscore, n-2)
    result =  not(alpha < tp < 1-alpha)
    return result

class TimeSeriesNetwork(object):
    def __init__(self, x, labels, alpha=0.1):
        self.cm = np.corrcoef(x)
        self.k = x.shape[1]
        self.alpha = alpha
        self.labels = labels
        self._make_graph()
        
    def _make_graph(self):
        self.graph = Graph()
        self.graph.add_vertices(self.labels)
        for i in range(len(self.cm)):
            for j in range(i+1, len(self.cm)):
                if _is_sig(self.cm[i][j], self.k, self.alpha):
                    self.graph.add_edge(self.graph.vertices()[i], self.graph.vertices()[j])
    def plot(self):
        self._make_graph()
        if len(self.cm) > 2:
            return self.graph.plot(layout='circular')
        return self.graph.plot()
