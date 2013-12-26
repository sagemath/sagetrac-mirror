This directory contains initialization files for the shell to open if
you run "sage -sh". Note that we want the following to be true:

* User settings (from ~/.bashrc, for example) are available in the
  Sage shell.

* Anything that conflicts with Sage is overwritten. I particular, it
  is your responsibility to source sage-env to set up the Sage
  environment variables. There is an environment variable
  _SAGE_ENV_SCRIPT exported to help you find it. 

The best thought-out example is probably the bash initialization
file. Note that it sources the user's .bashrc, and then modifies it as
necessary. You should follow the same pattern with other shells.

Some shells don't allow us to change the startup file. In that case we
can't help you if you overwrite the PATH in your startup files.



