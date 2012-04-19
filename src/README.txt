# 1. Make a clone of the git repository anywhere you want: 
     git clone  git@github.com:haikona/OMS.git

# 2. Apply a patch to the Sage library:
     sage -sh
     cd $SAGE_ROOT/devel/sage
     hg qimport /path/to/changes_to_sagelib.patch
     hg qpush

# 3. Edit $SAGE_ROOT/local/include/zn_poly/zn_poly.h
     sage -sh
     cd $SAGE_ROOT/local/include/zn_poly
     emacs zn_poly.h

     DELETE line 72:  typedef unsigned long ulong
     ADD three lines: #ifndef ulong
                      #define ulong unsigned long
		      #endif

## The following will not work, since dist.pyx uses a relative link.
# 4. Make symlinks:
     sage -sh
     cd $SAGE_ROOT/devel/sage-main/sage/modular/
     ln -s /path/to/OMS/sage/modular/btquotients .
     ln -s /path/to/OMS/sage/modular/pollack_stevens .
     cd $SAGE_ROOT/devel/sage-main/sage/modular/overconvergent/
     ln -s /path/to/OMS/sage/modular/overconvergent/pollack .

## Instead, you can do the following and develop in place
# 4'. Copy files from your clone

      sage -sh
      cd $SAGE_ROOT/devel/sage-main/
      cp -r /path/to/OMS/.git .git
      rm README.txt
      cp /path/to/OMS/README.txt .
      cp /path/to/OMS/changes_to_sagelib.patch .
      cd sage/modular/
      cp -r /path/to/OMS/sage/modular/btquotients btquotients
      cp -r /path/to/OMS/sage/modular/pollack_stevens pollack_stevens
      cd overconvergent/
      cp -r /path/to/OMS/sage/modular/overconvergent/pollack pollack

## You should run git commands in /path/to/OMS if you went with option 4, or $SAGE_ROOT/devel/sage-main in option 4'.
# 5. Synchronize

     git pull
     git push

If changes_to_sagelib.patch changes, do this:

     sage -sh
     cd $SAGE_ROOT/devel/sage/
     hg qpop
     hg qrm   changes_to_sagelib.patch
     hg qimport /path/to/changes_to_sagelib.patch
     hg qpush
     sage -br
