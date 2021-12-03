#!/usr/bin/env python
import sys;
import os;

pwd=os.getcwd()

f=open("configfile");
for line in f:
    linecont=line.strip("\n")[1:];
    sys.path.append(linecont);

import cosmology as c;

aa=c.cosmology();

print aa.nofm(1E10,0.1);
print aa.nofm(1E10,0.5);

print "NOTE: If the code printed out: "
print "1.94184058513e-11"
print "2.02622771397e-11"
print "The code is properly setup"

print "Checking that both comoving distances are equal"
print "Chi(3.035):",aa.Chiofz(3.035),"Chi_num(3.035):",aa.Chiofz_num(3.035);
print "Chi(5.012):",aa.Chiofz(5.012),"Chi_num(5.012):",aa.Chiofz_num(5.012);
print "\n"
print "#######################################################"
print "USE THE FOLLOWING LINES TO REPLACE THE INPUT*.PY FILES IN CLUS_LENS and GAL_LENS"
print "################"

print "ff=open(\""+pwd+"/configfile\");"
print "sys.path.append(\""+pwd+"/\"+linecont);"

print "#######################################################"
