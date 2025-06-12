#!/bin/bash

# Load CMS grid environment - only on DESY NAF
if [ -f "/cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh" ]; then
    source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh > /dev/null
fi
# Set custom VOMS proxy path. This needs to be accessible from Condor
# Please do not forget to run voms-proxy-init --voms cms --out $X509_USER_PROXY
export X509_USER_PROXY=~/.globus/x509up
# Load LCG
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107a x86_64-el9-gcc14-opt

# Load your virtual environment - comment this out if you use one!
# source /INSERT/ABSOLUTE/PATH/TO/YOUR/ENVIRONMENT/bin/activate

if test -n "$VIRTUAL_ENV"; then
    # Make sure python libs installed in virtual environment are prefered over system-wide ones
    VENV_SITE_PACKAGES="$(realpath ${VIRTUAL_ENV}/lib/python3.*/site-packages)"
    export PYTHONPATH=$VENV_SITE_PACKAGES:$PYTHONPATH
else
    # Make sure python libs installed in the user directory are prefered over system-wide ones
    export PYTHONPATH=`python3 -c 'import site; print(site.getusersitepackages())'`:$PYTHONPATH
fi
# On DESY NAF, old HDF5 plugins that are installed system-wide break HDF5 functionality. Disable
unset HDF5_PLUGIN_PATH
# Use this script also as environment script when running Pepper on HTCondor
if test -n "$BASH_VERSION"; then
    export PEPPER_CONDOR_ENV="$(realpath $BASH_SOURCE)"
elif test -n "$ZSH_VERSION"; then
    export PEPPER_CONDOR_ENV="$(realpath ${(%):-%N})"
fi
