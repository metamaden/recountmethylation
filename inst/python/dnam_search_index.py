#!/usr/bin/env python3

# Author: Sean Maden
#
# Manage DNAm array search indexes. Options to make an HNSW-based 
# search index from hashed features tables of DNAm array data.
# 

import mmh3
import numpy as np
import pandas as pd
import hnswlib, sys, os, re, pickle, random
from time import time
import faulthandler
faulthandler.enable()

def feature_hash(arr, target_dim=10000):
    """ Perform feature hashing on an array of data
    
    Perform feature hashing on the data in arr, into a vector of target_dim 
    total hashed features.

    Arguments:
        * arr: An array of values to be hashed.
        * target_dim: The target number of hashed values.
    Returns:
        * low_d_rep, or an array of hashed values of len == target_dim

    """ 
    low_d_rep = [0.0 for _ in range(target_dim)]
    for i, el in enumerate(arr):
        hashed = mmh3.hash(str(i))
        if hashed > 0.0:
            low_d_rep[hashed % target_dim] += arr[i]
        else:
            low_d_rep[hashed % target_dim] -= arr[i]
    return low_d_rep

def make_fhmatrix_autolabel(wf_name, of_name, lnotfloat = ['','NA','NaN'], 
    ndim = 10000, lstart = 0):
    """
    
    Get the hashed features table from an input data table. This function 
    automatically sets row labels as the first column of data in the input data
    table at `of_name.`

    Arguments:
        * wf_name: Name/path of hashed features table to write (required, 
            string, rows = samples, cols = hashed features).
        * of_name: Name/path of table to hash (required, tring, rows =  samples, 
            cols = probes).
        * lnotfloat: List of expected missing value symbols (required, ['','NA',
            'NaN']). These are replaced by the row-wise meidans of non-missing 
            values.
        * ndim: Number of hashed features (integer, 1000).
        * lstart: Line to start reading (required, int, 0).
    Returns:
        * None, saves new hashed features table to `wf_name`.

    """
    with open(of_name, "r") as fr:
        with open(wf_name, "w") as fw:
            for li, line in enumerate(fr):
                line_format = line.replace('\n', '').split(',')
                newrow = line_format[0]; lli = line_format[1::]
                if li >= lstart:
                    # replace NAs with median values
                    lli_median = np.median([float(ii) for ii in lli 
                        if not ii in lnotfloat])
                    lli_format = [float(ii) if not ii in lnotfloat
                                    else lli_median for ii in lli]
                    lli_fh = feature_hash(lli_format, target_dim = ndim)
                    newrow = newrow + ',' + ','.join([str(ii) for ii in lli_fh])
                    newrow = newrow + '\n'
                    print('Found new sample: '+newrow[0:100])
                    fw.write(newrow)    
                print("Finished with line number "+str(li))
    return None

def make_fhmatrix_specifylabels(labels_list, wf_name, of_name, 
    lnotfloat = ['','NA','NaN'], ndim = 10000):
    """ Make the hashed feaures matrix a Beta-values table

    Saves a .csv table of ndim total hashed features (columns) by elements, 
    where the number of elements is equal to the row count in the of_name table.
    Missing/NA values are replaced by the row median. New rows are processively
    written to the file wf_name.

    Arguments:
        * labels_list : Column names of the output table
        * wf_name : Name of the output table
        * of_name : Name of table to be hashed. Rows will be hashed. First line 
                    is colnames and thus skipped.
    Returns:
        * None, produces a new hashed features .csv table of dim nrow_of_name by 
        ndim+1 (columns), where first column has feature labels


    """
    with open(of_name, "r") as fr:
        with open(wf_name, "w") as fw:
            for li, line in enumerate(fr):
                if li > 0:
                    lli = line.replace('\n', '').split(',')[1::]
                    newrow = labels_list[li] # append label to new row
                    # replace NAs with median values
                    lli_median = np.median([float(ii) for ii in lli 
                        if not ii in lnotfloat])
                    lli_format = [float(ii) if not ii in lnotfloat
                                    else lli_median for ii in lli]
                    lli_fh = feature_hash(lli_format, target_dim = ndim)
                    newrow = newrow + ',' + ','.join([str(ii) for ii in lli_fh])
                    newrow = newrow + '\n'
                    print('Writing new row: '+newrow[0:100])
                    fw.write(newrow)
                    print("Finished with line number "+str(li))
    return None

def make_hnsw_si(fname, index_name, dindex_name, space_val = 'l2', threads = 10, 
    efc_val = 2000, m_val = 1000, ef_val = 2000):
    """ Makes a new search index from a csv table

    Arguments:
        * fname: Name of CSV file containing index data.
        * index_name: Name of new search index file to save.
        * dindex_name: Name of new index dict, with str labels, to save.
        * space_val: Space value for the index.
        * threads: Number of threads for processing the index.
        * efc_val: EFC value for the index.
        * m_val: M value for the index.
        * ef_val: EF value for the index.
    Returns
        * None, saves the new index and index dict

    """
    print("Loading the dataset...")
    df=pd.read_csv(fname, sep=',',header=None)
    df_str_labels = df.iloc[0::,0]
    df = df.iloc[0::,1::] # subset to exclude string labels
    num_elements = df.shape[0]
    dim = df.shape[1]
    print("Making the new index...")
    p_iter = hnswlib.Index(space = space_val, dim = dim)
    p_iter.set_num_threads(threads)
    df_labels = np.arange(num_elements)
    p_iter.init_index(max_elements = num_elements, 
        ef_construction = efc_val, M = m_val)
    p_iter.set_ef(ef_val)
    print("Adding new data to the search index...")
    p_iter.add_items(df, df_labels)
    print("Saving the search index...")
    p_iter.save_index(index_name)
    print("Saving the index dict with str labels, as pickle object...")
    dindex = {'index' : p_iter, 'strlabels' : df_str_labels}
    file_open = open(dindex_name, "wb")
    pickle.dump(obj = dindex, file = file_open)
    return None

def query_si(query_data, si_fname, si_labels = [], kval = 2):
    """ Perform a lookup on an hnswlib-saved search index

    Arguments:
        * query_data : Vector of feature hashed data, of dims R x C, where C is the 
                same as in siobject/search index object.
        * si_fname : Name/path to the queried search index.
        * si_labels : Sample ID labels corresponding to index positions in 
            search index at `si_fname`.
        * kval : K number of nearest neighbors to return from siobject lookup.

    Returns:
        * Vector of elements, where elements are indices if si_labels = False, 
                or else element labels if si_labels is provided.

    """
    time_start = time()
    print("Loading search index...")
    sidict = pickle.load(open(si_fname, "rb"))
    siobject = sidict["index"]
    print("Querying "+str(query_data.shape[0])+" elements in data with k = "+
        str(kval)+" nearest neighbors...")
    kval_labels, kval_distances = siobject.knn_query(query_data, k = kval)
    time_str = str(time()-time_start);print("Query completed, time: "+time_str)
    if len(si_labels) > 0:
        print("Applying labels to query results..."); kval_str_labels = []
        for ll in kval_labels:
            kval_str_labels.append([si_labels[ii] for ii in ll])
        print("Returning data (sample id, k index, and distance)...")
        return kval_str_labels, kval_labels, kval_distances
    else:
        print("Returning data (k index and distance)...")
        return kval_labels, kval_distances
    return NULL

def make_dfk_sampleid(sample_idv, lk = [1,2], fh_csv_fname = "bval_100_fh10.csv",
    index_dict_fname = "new_index_dict.pickle"):
    """ make_dfk_sampleid

    Get the sample labels for the k nearest neighbors on a series of queries.
    Use a vector gsmv to identify samples for the query, from metadata.
    
    Arguments:
        * sample_idv : Vector or list of sample ID strings, corresponding to 
            sample ID labels in the rownames/first column of the hashed features 
            table at `fh_csv_fname`, which can correspond to sample names in the 
            queried search index (requried, list/vector of strings).
        * lk : List of k nearest neighbors to query (required, list of int 
            values, 1000).
        * index_dict_fname: Name/path of the search index file (required, 
            string).
        * fh_csv_fname: Name/path of the hashed features csv to read (required, 
            string, )
    Returns:
        * dfk_final

    """
    print("Getting hashed features data for samples...")
    fhdict = {}
    with open(fh_csv_fname, "r") as of:
        for line in of:
            lline = line.split(",")
            sample_id = lline[0].replace('"', '').split(".")[0]
            if sample_id in sample_idv:
                print("Getting index data for sample: '{0}'".format(sample_id))
                fhdict[sample_id] = [float(fhi.replace("\n", "")) 
                    for fhi in lline[1::]
                ]
    dfi = pd.DataFrame.from_dict(fhdict).T
    si_dict = pickle.load(open(index_dict_fname, "rb"))
    dim_si_keylabv = len([ii for ii in si_dict['strlabels']])
    dfk = pd.DataFrame()
    dfk["sample_id"] = fhdict.keys()
    for ii, ki in enumerate(lk):
        print("Querying ",ki," neighbors from lk...")
        if ki <= dim_si_keylabv:
            kval_str_labels, kval_labels, kval_distances = query_si(
                query_data = dfi, si_fname = index_dict_fname, 
                si_labels = si_dict["strlabels"], kval = ki
            )
            kli_format = [";".join(ii) for ii in kval_str_labels]
            dfki = pd.DataFrame(kli_format)
            dfk['k=' + str(ki)] = [ii for ii in dfki[0]]
        else:
            print("Provided k '{0}' > n si samples, skipping...".format(ki))
    print("Returning query results...")
    return dfk

def main():
    """
    """