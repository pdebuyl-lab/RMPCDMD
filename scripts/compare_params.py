#!/usr/bin/env python
r"""
Compare the arguments listed by Doxygen's \param directive and those read using
PTread_* functions. This program is use to check the consistency of the
documentation of the programs.
"""
from __future__ import print_function, division

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file')
parser.add_argument('--print-line', action='store_true')
args = parser.parse_args()

import re

doxygen_parameters = []
PT_parameters = []
with open(args.file, 'r') as f:
    for l in f.readlines():
        doxy_match = re.search(r'!! \\param (\w+) ', l)
        if doxy_match:
            if args.print_line: print(l, end='')
            doxygen_parameters.append(doxy_match.group(1))
        PT_match = re.search('PTread_\w+\(\s*config\s*,\s*\'(\w+)\'[,\)]', l)
        if PT_match:
            if args.print_line: print(l, end='')
            PT_parameters.append(PT_match.group(1))

doxygen_parameters.sort()
PT_parameters.sort()

print(doxygen_parameters)
print(PT_parameters)
print(doxygen_parameters == PT_parameters)

