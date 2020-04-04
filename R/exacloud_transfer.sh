#!/usr/bin/bash

# get onto exacloud
ssh maden@acc.ohsu.edu
ssh exacloud.ohsu.edu
cd /home/exacloud/lustre1/CompBio/recount_meth


# transfer calls
scp -r -P 21747 metamaden@96.79.111.125:/home/metamaden/recount-methylation/recount-methylation-analysis/files/mdata/compilations/ ./


scp -r -P 21747 metamaden@96.79.111.125:/home/metamaden/recount-methylation/recountmethylation_manuscript_supplement ./

scp -r -P 21747 metamaden@96.79.111.125:/home/metamaden/recount-methylation/recount-methylation-analysis/files/mdata/compilations/remethdb_h5se_rg_00-00-01_1583780004 ./
scp -r -P 21747 metamaden@96.79.111.125:/home/metamaden/recount-methylation/recount-methylation-analysis/files/mdata/compilations/remethdb_h5se_gr_00-00-01_1583780004 ./
scp -r -P 21747 metamaden@96.79.111.125:/home/metamaden/recount-methylation/recount-methylation-analysis/files/mdata/compilations/remethdb_h5se_gm_00-00-01_1583780004 ./
