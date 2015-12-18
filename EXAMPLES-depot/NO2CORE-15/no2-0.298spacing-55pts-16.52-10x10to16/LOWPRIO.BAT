perl -pi -e "s/condo_axl/lr_lowprio/" run.axl.pulse.*
perl -pi -e "s/nodes=7/nodes=11/" run.axl.pulse.*
perl -pi -e "s/504/72/" run.axl.pulse.*
perl -pi -e "s/lr3 --constraint=lr3_c20/lr2/" run.axl.pulse.*

perl -pi -e "s/NPN=18/NPN=11/" Pulse.*Bat*
perl -pi -e "s/mctdhf_sinc/mctdhf_sinc.lr2/" Pulse.*Bat*

