# -*- shell-script -*-
#
# The Sage zsh Startup File
#
# See README.txt for some background considerations.

# Sanity check
if [ -z "$_SAGE_ENV_SCRIPT" ] ; then
    echo >&2 "Error, can only be used from within the main sage script."
fi

# Set _SAGE_ENV_SCRIPT readonly so .zshrc cannot modify it.
readonly _SAGE_ENV_SCRIPT

# Source user definitions
if [ -f ~/.zshrc ]; then
    . ~/.zshrc
fi

# Now overwrite settings with Sage environment variables
source "$_SAGE_ENV_SCRIPT"

# Prompt setup
PS1="%S(sage-sh)%s %n@%m:%~$ "
export PS1
