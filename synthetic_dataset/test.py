#!/usr/bin/env python

import os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
os.sys.path.insert(0,parentdir) 
print currentdir
print parentdir
from lib.foo import print_foo
print_foo()
