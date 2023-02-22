#!/usr/bin/env python

import numpy as np
from irescue.misc import getlen, writerr, flatten, run_shell_cmd
import gzip
import os

# calculate mismatches between sequences of same length
def find_mm(x, y):
    if len(x) != len(y):
        return -1
    mm = 0
    for i in range(len(x)):
        if x[i] != y[i]:
            mm += 1
    return mm

# collapse UMI graph
def collapse_networks(graph):
    out = dict()
    for key, value in graph.items():
        out[key] = []
        if len(value) == 1:
            # check if it's a single node, then add to the output and go to the next node
            out[key].append(value)
            continue
        for val in graph.values():
            # check if there are other nodes that contains all the values of the current one
            if all(i in value for i in val):
                out[key].append(val)
        if len(out[key]) <= 1:
            out.popitem()
    return out

# calculate counts of a cell from mappings dictionary
def cellCount(maps, intcount = False):

    # get and index equivalence classes from maps
    eclist = list()
    for v in maps.values():
        eclist.append(tuple(sorted(v.keys())))
    eclist = list(set(eclist))

    # make a simple mapping dict (index number in place of families) and its reverse
    smaps = dict([(i,eclist.index(tuple(sorted(j.keys())))) for i,j in maps.items()])
    rsmaps = dict()
    for key, value in smaps.items():
        rsmaps.setdefault(value, list()).append(key)

    # compute the count of each equivalence class in the cell barcode
    counts = dict()
    for ec in rsmaps:
        # list of UMIs associated to EC
        umis = rsmaps[ec]
        ### compute the total count of the equivalence class
        if len(umis) > 1:
            ### Find and collapse duplicated UMIs ###
            # Make an NxN array of number of mismatches between N UMIs
            mm_arr = np.array([[find_mm(ux,i) for i in umis] for ux in umis])
            # Find UMI pairs with up to 1 mismatch, where UMIs are representad
            # by integers: [[i, j], [i, k], [k, m]]
            mm_check = np.argwhere(mm_arr <= 1)
            # Make a graph that connects UMIs with <=1 mismatches
            # {NODE: EDGES} or {UMI: [CONNECTED_UMIS]}
            graph = dict()
            for key, value in mm_check:
                graph.setdefault(key, set()).add(value)
            # Check if all nodes are connected (i.e. complete graph)
            if all([x == set(graph.keys()) for x in graph.values()]):
                # Set EC final count to 1
                ec_count = 1
            else:            
                # Collapse networks based on UMI similarity: {HUB: [UMI_GRAPHS]}
                coll_nets = collapse_networks(graph)
                # Get EC final count after collapsing
                ec_count = len(coll_nets)
        else:
            # If only one umi, skip collapsing and assign 1 to the final count
            ec_count = 1

        ### find the predominant TE family in the equivalence class
        # make count matrix from mappings (row = UMI, column = TE)
        ec_counts = np.array([list(j.values()) for i,j in maps.items() if i in rsmaps[ec]])
        # sum counts by TE
        ec_sum = ec_counts.sum(axis = 0)
        # find the index of the highest count
        ec_max = np.argwhere(ec_sum == ec_sum.max()).flatten()

        # retrieve the TEs with highest count
        te_max = list()
        for i in ec_max:
            te_max.append(eclist[ec][i])

        # add count
        for te in te_max:
            # initialize the feature in the cell barcode dictionary
            if te not in counts:
                counts[te] = 0
            # get the normalized count by dividing the raw count by the number of predominant TEs
            norm_count = ec_count / len(te_max)
            # if integers are needed, round the normalized count
            if intcount:
                norm_count = round(norm_count)
            # add count to dictionary
            counts[te] += norm_count

    return counts

def parse_features(features_file):
    '''Generator for (index,feature) tuples'''
    with gzip.open(features_file, 'rb') as f:
        for i, line in enumerate(f):
            l = line.decode('utf-8').strip().split('\t')
            yield (l[0], i+1)

def split_int(num,div):
    '''Splits an integer x into n integers whose sum is equal to x'''
    split = int(num/div)
    for i in range(0, num, split):
        j = i + split
        if j > num-split:
            j = num
            yield range(i, j)
            break
        yield range(i, j)

def split_bc(barcode_file, n):
    '''Yields barcodes (index,sequence) tuples in n chunks'''
    bclen = getlen(barcode_file)
    #split = round(bclen/n)
    with gzip.open(barcode_file, 'rb') as f:
        c=0
        for chunk in split_int(bclen, n):
            yield (c,[(next(f).decode('utf-8').strip(),x+1) for x in chunk])
            c+=1

def count(mappings_file, outdir, tmpdir, features, intcount, verbose, bc_split):
    '''Runs cellCount for a set of barcodes'''
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)

    # set temporary matrix name prefix as chunk number
    chunkn = bc_split[0]
    matrix_file = tmpdir + f'/{chunkn}_matrix.mtx.gz'

    # parse barcodes in a SEQUENCE:INDEX dictionary
    barcodes = dict(bc_split[1])
    writerr(f'Processing {len(barcodes)} barcodes from chunk {chunkn}', send=verbose)

    # get number of lines in mappings_file
    nlines = getlen(mappings_file)

    # initialize mappings dictionary {UMI: {FEATURE: COUNT}}
    maps = dict()

    # cell barcode placeholder
    cell = None

    with gzip.open(mappings_file, 'rb') as data, \
    gzip.open(matrix_file, 'wb') as mtxFile:

        for line in enumerate(data, start=1):
            # gather barcode, umi and feature from mappings file
            cx, ux, te = line[1].decode('utf-8').strip().split('\t')
            if '~' in te:
                te = te[:te.index('~')]

            if len(barcodes)==0:
                # interrupt loop when reaching the end of the barcodes chunk
                break

            if not cell:
                # skip to the first cell barcode contained in the current barcodes chunk
                if cx not in barcodes:
                    continue
                else:
                    cell = cx

            # if cell barcode changes, compute counts from previous cell's mappings
            if cx != cell and cell in barcodes:
                cellidx = barcodes.pop(cell)
                writerr(f'[{chunkn}] Computing counts for cell barcode {cellidx} ({cell})', send=verbose)
                # compute final counts of the cell
                counts = cellCount(maps, intcount=intcount)
                # arrange counts in a data frame and write to text file
                lines = [ f'{str(k)} {str(cellidx)} {str(v)}\n'.encode() \
                        for k, v in counts.items() ]
                mtxFile.writelines(lines)
                # re-initialize mappings dict
                maps = dict()
            
            # reassign cell to current barcode
            cell = cx

            # add features count to mappings dict
            if cx in barcodes:
                teidx = features[te]
                if ux not in maps:
                    # initialize UMI if not in mappings dict
                    maps[ux] = dict()
                if teidx in maps[ux]:
                    # initialize feature count for UMI
                    maps[ux][teidx]+=1
                else:
                    # add count to existing feature in UMI
                    maps[ux][teidx]=1

            # if end of file is reached, compute counts from current cell's mappings
            if line[0] == nlines and cell in barcodes:
                cellidx = barcodes.pop(cell)
                writerr(f'[{chunkn}] [file_end] Computing counts for cell barcode {cellidx} ({cell})', send=verbose)
                # compute final counts of the cell
                counts = cellCount(maps, intcount=intcount)
                # arrange counts in a data frame and write to text file
                lines = [ f'{str(k)} {str(cellidx)} {str(v)}\n'.encode() \
                        for k, v in counts.items() ]
                mtxFile.writelines(lines)

    writerr(f'Barcodes chunk {chunkn} written to {matrix_file}', send=verbose)
    return matrix_file

# Concatenate matrices in a single MatrixMarket file with proper header
def formatMM(matrix_files, outdir, features, barcodes):
    if type(matrix_files) is str:
        matrix_files = [matrix_files]
    matrix_out = outdir + '/matrix.mtx.gz'
    features_count = len(features)
    barcodes_count = len(flatten([j for i,j in barcodes]))
    mmsize = sum(getlen(f) for f in matrix_files)
    mmheader = '%%MatrixMarket matrix coordinate real general\n'
    mmtotal = f'{features_count} {barcodes_count} {mmsize}\n'
    with gzip.GzipFile(matrix_out, 'wb', mtime=0) as mmout:
        mmout.write(mmheader.encode())
        mmout.write(mmtotal.encode())
    mtxstr = ' '.join(matrix_files)
    cmd = f'zcat {mtxstr} | LC_ALL=C sort -k2,2n -k1,1n | gzip >> {matrix_out}'
    run_shell_cmd(cmd)
    return(matrix_out)
