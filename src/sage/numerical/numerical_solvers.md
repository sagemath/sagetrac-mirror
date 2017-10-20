# What we've done

## Basic polynomial homotopy continuation classes (polynomial systems, numerical points, etc.)
## Abstract engine classes (PHCpackEngine, HomotopyContinuationEngine)


# Big picture to do

## Expand coverage/everything needs a doctest
## Good documentation (examples, etc.)
## Thematic tutorial (distant future)
## TestEngine: add in a symbolic solver, whch will work even if phcpy isn't installed. This will help drive out bugs.
## Add more features to solvers.
   -PHCpackEngine still needs:
       -path tracking (taking in a Homotopy object)
       -anything positive dimensional
## Better integration between engines and homotopy types
## Necessary function that lives in the homotopy engine file: take in a string, return an instance of the desired class.


# Small picture to do

## Sort out circular imports
## Get solvers + engine to the point where they can replicate examples in phcpy tutorial
## Make many of the polynomial_homotopy types' data attributes to private + add methods for access
## WitnessSet -- isIrreducible method
## I think some __repr__ functions in polynomial_homotopy probably need to be changed
## Get/maintain good coverage/pylint scores, make sure test and build work
## storing fibres that have been solved on parametrized polynomial systems
