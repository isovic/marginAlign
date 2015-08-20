#!/bin/bash

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

#Get directory containing shell script
BASEDIR=$(dirname $0)/../

# Check if virtualenv has been created, and if it exists, activate it
if [ -d ${BASEDIR}/env/bin/ ]; then
# Control will enter here if $DIRECTORY exists.
source ${BASEDIR}/env/bin/activate
# set a virtenv flag so it can be deactivated once done
virtenv=true
fi

#Set the python path to just the relevant directories
export PYTHONPATH=${BASEDIR}/src:${BASEDIR}/submodules

#Set the path to just the relevant directories
export PATH=${BASEDIR}/submodules/sonLib/bin:${BASEDIR}/submodules/jobTree/bin:${BASEDIR}/submodules/bwa/:${BASEDIR}/submodules/last/src/:${BASEDIR}/submodules/last/scripts/:${BASEDIR}/submodules/graphmap/bin/Linux-x64/:${PATH}

#Run described script
python $@

# Deactivate virtualenv if it was activated
if [ $virtenv ]; then
# deactivate virtualenv
deactivate
fi
