#!/bin/sh
csplit -z -n 1 $1 "/Copyright/"-1 {*}
mv xx0 btquotient.py
mv xx1 pautomorphicform.py
mv xx2 ocmodule.py
mv xx3 utility.py

#printf '0a\nfrom sage.modular.btquotients.utility import our_log,our_exp\n.\nw\n' | ed ocmodule.py
#printf '0a\nfrom sage.modular.btquotients.utility import our_sqrt\n.\nw\n' | ed btquotient.py
printf '0a\nfrom sage.modular.btquotients.btquotient import *\nfrom sage.modular.btquotients.ocmodule import *\n.\nw\n' | ed pautomorphicform.py
