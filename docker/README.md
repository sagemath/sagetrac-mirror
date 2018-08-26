[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://github.com/sagemath/sage/COPYING.txt) [![Maintained](https://img.shields.io/maintenance/yes/2018.svg)](https://github.com/sagemath/sage/commits/master)

# Supported tags

* `latest` — the stable `master` branch [![GitHub last commit (branch)](https://img.shields.io/github/last-commit/sagemath/sage/master.svg)](https://github.com/sagemath/sage/commits/master) [![CircleCI branch](https://img.shields.io/circleci/project/github/sagemath/sage/master.svg)](https://circleci.com/gh/sagemath/sage/tree/master) [![GitLab CI](https://gitlab.com/sagemath/sage/badges/master/pipeline.svg)](https://gitlab.com/sagemath/sage/commits/master)
* `x.x` — all stable releases of Sage are tagged with their version number.
* `x.x.{beta,rc}x` - betas and release candidates of Sage as [tagged in our git repository](https://github.com/sagemath/sage/tags).
* `develop` — the current development version of Sage which gets merged into the `master` branch when a new version of Sage is released [![GitHub last commit (branch)](https://img.shields.io/github/last-commit/sagemath/sage/develop.svg)](https://github.com/sagemath/sage/commits/develop) [![CircleCI branch](https://img.shields.io/circleci/project/github/sagemath/sage/master.svg)](https://circleci.com/gh/sagemath/sage/tree/master) [![GitLab CI](https://gitlab.com/sagemath/sage/badges/develop/pipeline.svg)](https://gitlab.com/sagemath/sage/commits/develop)


# What is SageMath

SageMath is a free open-source mathematics software system licensed under the GPL. It builds on top of many existing open-source packages: NumPy, SciPy, matplotlib, Sympy, Maxima, GAP, FLINT, R and many more. Access their combined power through a common, Python-based language or directly via interfaces or wrappers. 

**Mission**: *Creating a viable free open source alternative to Magma, Maple, Mathematica and Matlab.*

# What's in this image

There are several flavours of this image.

* [`sagemath/sagemath`![image size](https://img.shields.io/microbadger/image-size/sagemath/sagemath/latest.svg)](https://hub.docker.com/r/sagemath/sagemath) contains everything necessary to run Sage on the command line. Run it with:
    ```
    docker run -it sagemath/sagemath:latest
    ```
    You can start a graphical [Jupyter Notebook](https://jupyter.org) at http://localhost:8888 instead. To use the notebook, follow the instructions printed when you run:
    ```
    docker run -p8888:8888 sagemath/sagemath:latest sage-jupyter
    ```
* [`sagemath/sagemath-dev`![image size](https://img.shields.io/microbadger/image-size/sagemath/sagemath-dev.svg)](https://hub.docker.com/r/sagemath/sagemath-dev) contains all the build artifacts to rebuild Sage quickly. This version is probably only relevant for Sage developers. Run this image with:
    ```
    docker run -it sagemath/sagemath-dev:develop
    ```
    This triggers a rebuild and drops you in a shell afterwards. Note that the git repository has been emptied to save space. If you want to use git, fetch from your git repository with `git fetch trac` and go to the commit that was used to create this image with
    ```
    git reset $(cat docker/.commit)
    ```

# How to build your own SageMath images

Run `docker build -f docker/Dockerfile --build-arg ARTIFACT_BASE=sagemath/sagemath-dev:develop --target TARGET .` in the Sage repository with `TARGET` one of `sagemath` or `sagemath-dev`.

# How these images get updated

Every push to our [github repository](https://github.com/sagemath/sage) triggers a build in [CircleCI](https://circleci.com) which builds and pushes the docker images.
A push to master also triggers a "build" on our [Docker Hub](https://hub.docker.com) repositories. The latter build is mostly disabled by the `hooks/` and only updates the `README.md`.

Every push to our [GitLab repository](https://gitlab.com/sagemath/sage) triggers a pipeline in GitLab CI. This build also pushes images to Docker Hub.

Have a look at `.circleci/` and `.gitlab-ci.yml` if you want to setup either continuous integration service for your own fork of the SageMath repository.

# Report bugs and issues

Please tell us of any bugs or omissions at our [Issue Tracker](https://trac.sagemath.org) or contact us through the [sage-support](https://groups.google.com/forum/#!forum/sage-support) or the [sage-devel](https://groups.google.com/forum/#!forum/sage-devel) mailing lists.

# License

The whole Sage software distribution is licensed under the General Public License, version 3. More details can be found in our [COPYING.txt](https://github.com/sagemath/sage/blob/master/COPYING.txt)

[//]: # (Please don't break long lines in this files as dockerhub then gets the formatting of this file wrong.)
