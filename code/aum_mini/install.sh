#!/bin/bash

rm -rf build/ cosmic/
python setup.py build_ext
python setup.py install --root=`pwd`/cosmic

PPATH=`find . -iname aum-1.01.egg-info | sed s/.// | sed s/aum-1.01.egg-info//`
echo $PPATH > configfile
echo "NOTE: Check this output "
echo "***************************"
python test.py
echo "***************************"
