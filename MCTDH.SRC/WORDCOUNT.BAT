
 cat *.f90 COREPROJECT/*90 H2PROJECT/*90 HEPROJECT/*90 SINCDVR/*90 |sed '/^\s*!/d;/^s*$/d'|wc
