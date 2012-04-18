# 1. Make a clone of the git repository anywhere you want: 
     git clone  git@github.com:haikona/OMS.git

# 2. Apply a patch to the Sage library:
     sage -sh
     cd $SAGE_ROOT/devel/sage
     hq qimport changes_to_sagelib.patch
     hg qpush

# 3. Make symlinks:
     sage -sh
     cd $SAGE_ROOT/devel/sage-main/sage/modular/
     ln -s /path/to/OMS/sage/modular/btquotients .
     ln -s /path/to/OMS/sage/modular/pollack_stevens .
     cd $SAGE_ROOT/devel/sage-main/sage/modular/overconvergent/
     ln -s /path/to/OMS/sage/modular/overconvergent/pollack .

# 4. Synchronize

     git pull
     git push

If changes_to_sagelib.patch changes, do this:

     sage -sh
     cd $SAGE_ROOT/devel/sage/
     hg qpop
     hg qrm   changes_to_sagelib.patch
     hq qimport /path/to/changes_to_sagelib.patch
     hg qpush
     sage -br
