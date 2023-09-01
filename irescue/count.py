#!/usr/bin/env python

from collections import Counter
import numpy as np
import networkx as nx
from irescue.misc import get_ranges, getlen, writerr, flatten, run_shell_cmd, iupac_nt_code
from irescue.em import run_em
import gzip
import os

def hdist(x, y):
    """
    Calculate hamming distance between sequences of the same length
    """
    if len(x) != len(y):
        writerr('Error: all UMI lengths must be the same.', error=True)
    return sum(1 for i in range(len(x)) if x[i] != y[i])

def connect_umis(x, y, threshold=1):
    """
    Check if UMI x connects to y.

    x, y : iterable
        (UMI, TE, count) : (str, set, int)
    """
    return (hdist(x[0], y[0]) <= threshold
            and x[2] >= (2 * y[2]) - 1
            and x[1].intersection(y[1]))

def pathfinder(graph, start_node, path=[], start_ft=None):
    """
    Finds first valid path of UMIs with compatible equivalence class given
    a starting node. Can be used iteratively to find all possible paths.
    """
    if not start_ft:
        start_ft = graph.nodes[start_node]['ft']
    else:
        start_ft = start_ft.intersection(graph.nodes[start_node]['ft'])
        if not start_ft:
            return path
    path = path + [start_node]
    for node in graph.successors(start_node):
        if start_ft.intersection(graph.nodes[node]['ft']) and node not in path:
            extended_path = pathfinder(graph, node, path, start_ft)
            if extended_path:
                return extended_path
    return path

def index_features(features_file):
    idx = {}
    with gzip.open(features_file, 'rb') as f:
        for i, line in enumerate(f, start=1):
            ft = line.strip().split(b'\t')[0]
        #    idx[i] = ft
            idx[ft] = i
    return idx

def parse_maps(maps_file, feature_index):
    """
    maps_file : str
        Content: "CB UMI FEATs count"
    out : str, list
        CB,
        [(UMI <str>, {FT <int>, ...} <set>, count <int>) <tuple>, ...]
    """
    with gzip.open(maps_file, 'rb') as f:
        eqcl = []
        it = f.readline().strip().split(b'\t')[0]
        f.seek(0)
        for line in f:
            cb, umi, feat, count = line.strip().split(b'\t')
            count = int(count)
            feat = {feature_index[ft] for ft in feat.split(b',')}
            if cb == it:
                eqcl.append((umi, feat, count))
            else:
                yield it, eqcl
                it = cb
                eqcl = [(umi, feat, count)]
        yield it, eqcl

def compute_cell_counts(equivalence_classes, number_of_features):
    """
    Calculate TE counts of a single cell, given a list of equivalence classes.

    Parameters
    ----------
    equivalence_classes : list
        (UMI_sequence, {TE_index}, read_count) : (str, set, int)
        Tuples containing UMI-TEs equivalence class infos.

    Returns
    -------
    out : dict
        feature <int>: count <float> dictionary.
    """
    # initialize TE counts
    counts = Counter()
    # build cell-wide UMI deduplication graph
    graph = nx.DiGraph()
    for i, eqc1 in enumerate(equivalence_classes):
        graph.add_node(i)
        graph.nodes[i]['ft'] = eqc1[1]
        graph.nodes[i]['count'] = float(eqc1[2])
        for j, eqc2 in enumerate(equivalence_classes):
            if i != j and connect_umis(eqc1, eqc2):
                graph.add_edge(i, j)
    # split cell-wide graph into subgraphs of connected nodes
    subgraphs = [graph.subgraph(x) for x in
                 nx.connected_components(graph.to_undirected())]
    # put aside networks that will be solved with EM
    em_array = []
    for subg in subgraphs:
        # find all parent nodes in graph
        parents = [x for x in subg if not list(subg.predecessors(x))]
        if not parents:
            # if no parents are found due to bidirected edges, take all nodes
            parents=list(subg.nodes)
        # initialize dict of possible paths configurations, starting from
        # each parent node.
        paths = {x: [] for x in parents}
        # find paths starting from each parent node
        for parent in parents:
            # populate this list with nodes utilized in paths
            blacklist = []
            # find paths in list of nodes starting from parent
            path = []
            subg_copy = subg.copy()
            nodes = [parent] + [x for x in subg_copy if x != parent]
            for node in nodes:
                # make a copy of subgraph and remove nodes already used
                # in a path
                if node not in blacklist:
                    path = pathfinder(subg_copy, node)
                    [blacklist.append(x) for x in path]
                    paths[parent].append(path)
                    [subg_copy.remove_node(x) for x in path]
        # find the path configuration leading to the minimum number of
        # deduplicated UMIs
        path_config = [
            paths[k] for k, v in paths.items()
            if len(v) == min([len(x) for x in paths.values()])
        ][0]
        # assign UMI count to features
        for path in path_config:
            # get list of nodes' features sets
            features_sets = [subg.nodes[x]['ft'] for x in path]
            # find the features common to all nodes
            feat = features_sets[0]
            for x in features_sets:
                feat = feat.intersection(x)
            feat = list(feat)
            if len(feat) == 1:
                # unambiguous assignment to feature
                counts[feat[0]] += 1
            elif len(feat) > 1:
                # add UMI-TE compatibility matrix to em_array
                row = [1 if x in feat else 0
                       for x in range(1, number_of_features+1)]
                em_array.append(row)
                # run EM to distribute count across features
                #array = []
                #for node in path:
                #    row = [0.0 for _ in feat]
                #    for i, f in enumerate(feat):
                #        if f in subg.nodes[node]['ft']:
                #            row[i] = subg.nodes[node]['count']
                #    array.append(row)
                #array = np.array(array)
                #abundances = run_em(array, cycles=100)
                #for i, f in enumerate(feat):
                #    counts[f] += abundances[i]
            else:
                print(nx.to_dict_of_lists(subg))
                print([subg.nodes[x]['ft'] for x in subg.nodes])
                print([subg.nodes[x]['count'] for x in subg.nodes])
                print(path_config)
                print(path)
                print(features_sets)
                print(feat)
                writerr("Error: no common features detected in subgraph's"
                        " path.", error=True)
    # run EM to optimize the assignment of UMI from multimapping reads
    em_array = np.array(em_array)
    em_counts = run_em(em_array, cycles=100)
    em_counts *= em_array.shape[1]
    for i, c in enumerate(em_counts, start=1):
        counts[i] += c
    return dict(counts)

def split_barcodes(barcodes_file, n):
    """
    barcode_file : iterable
    n : int
    ----------
    out : int, dict
    """
    nBarcodes = getlen(barcodes_file)
    with gzip.open(barcodes_file, 'rb') as f:
        for i, chunk in enumerate(get_ranges(nBarcodes, n)):
            yield i, {next(f).strip(): x+1 for x in chunk}

def run_count(maps_file, feature_index, tmpdir, verbose, barcodes_set):
    taskn, barcodes = barcodes_set
    matrix_file = os.path.join(tmpdir, f'{taskn}_matrix.mtx.gz')
    with gzip.open(matrix_file, 'wb') as f:
        for cellbarcode, cellmaps in parse_maps(maps_file, feature_index):
            if cellbarcode not in barcodes:
                continue
            cellidx = barcodes[cellbarcode]
            cellcounts = compute_cell_counts(
                equivalence_classes=cellmaps,
                number_of_features=len(feature_index)
            )
            lines = [f'{feature} {cellidx} {count}\n'.encode()
                     for feature, count in cellcounts.items()]
            f.writelines(lines)
    return matrix_file

def formatMM(matrix_files, feature_index, barcodes_chunks, outdir):
    if type(matrix_files) is str:
        matrix_files = [matrix_files]
    matrix_out = os.path.join(outdir, 'matrix.mtx.gz')
    features_count = sum(1 for _ in feature_index if type(_) is int)
    barcodes_count = sum(len(x) for _, x in barcodes_chunks)
    mmsize = sum(getlen(f) for f in matrix_files)
    mmheader = b'%%MatrixMarket matrix coordinate real general\n'
    mmtotal = f'{features_count} {barcodes_count} {mmsize}\n'.encode()
    with gzip.GzipFile(matrix_out, 'wb', mtime=0) as mmout:
        mmout.write(mmheader)
        mmout.write(mmtotal)
    mtxstr = ' '.join(matrix_files)
    cmd = f'zcat {mtxstr} | LC_ALL=C sort -k2,2n -k1,1n | gzip >> {matrix_out}'
    run_shell_cmd(cmd)
    return(matrix_out)
