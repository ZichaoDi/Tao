import os
import sys

topdir = r'./data'
for root, dirnames, filenames in os.walk(topdir):
    for filename in filenames:
           sys.path.append(os.path.abspath(root))
