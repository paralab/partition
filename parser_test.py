#!/usr/bin/env python 

import sys
import gmshparser

if len(sys.argv) < 2:
  print("Usage: " + sys.argv[0] + " file")
  exit

mesh = gmshparser.parse(sys.argv[1])


print(mesh)