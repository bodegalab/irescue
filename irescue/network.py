#!/usr/bin/env python

# NB: This module include partly modified third-party code distributed under the
# license below.

##############################################################################
# The MIT License (MIT)

# Copyright (c) 2015 CGAT

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##############################################################################

from collections import defaultdict

def get_substr_slices(umi_length, idx_size):
    '''
    Create slices to split a UMI into approximately equal size substrings
    Returns a list of tuples that can be passed to slice function
    '''
    cs, r = divmod(umi_length, idx_size)
    sub_sizes = [cs + 1] * r + [cs] * (idx_size - r)
    offset = 0
    slices = []
    for s in sub_sizes:
        slices.append((offset, offset + s))
        offset += s
    return slices

def build_substr_idx(equivalence_classes, length, threshold):
    '''
    Group equivalence classes into subgroups having a common substring
    '''
    slices = get_substr_slices(length, threshold+1)
    substr_idx = {k: defaultdict(set) for k in slices}
    for idx in slices:
        for ec in equivalence_classes:
            sub = ec.umi[slice(*idx)]
            substr_idx[idx][sub].add(ec)
    return substr_idx

def gen_ec_pairs(equivalence_classes, substr_idx):
    '''
    Yields equivalence classes pairs from build_substr_idx()
    '''
    for i, ec in enumerate(equivalence_classes, start=1):
        neighbours = set()
        for idx, substr_map in substr_idx.items():
            sub = ec.umi[slice(*idx)]
            neighbours = neighbours.union(substr_map[sub])
        neighbours.difference_update(equivalence_classes[:i])
        for nbr in neighbours:
            yield ec, nbr