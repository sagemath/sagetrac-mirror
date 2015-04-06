from sage.interfaces.cmdline.tool import (
    Tool, ToolAdapter, CompilerTool, ToolMissingException
)

# Import all submodules to make all instances known
import sage.interfaces.cmdline.posix
import sage.interfaces.cmdline.pdf
import sage.interfaces.cmdline.gcc
import sage.interfaces.cmdline.latex

assert not sage.misc.lazy_import.is_during_startup(), \
    'the command line interfaces must not be imported during startup'

Tool._disallow_further_instances()
