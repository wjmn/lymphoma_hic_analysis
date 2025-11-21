#!/usr/bin/env python
# Creates metadata for use with hitad from TADLib

import sys

input = sys.argv[1]
output = sys.argv[2]
res = sys.argv[3]


content = f"""
res:{res}
  rep1:{input}
"""

with open(output, "w") as outfile:
    outfile.write(content)
