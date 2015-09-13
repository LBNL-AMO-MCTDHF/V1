#!/bin/bash

export OMP_NUM_THREADS=1

 for ext in doub quart; do for sym in sig pi; do for ger in u g; do chmctdhf_diatom Inp=Input.Inp.Cation.$ext.${sym}_${ger} |tee Outs/Out.cation.$ext.${sym}_${ger} & done; done; done

 wait

