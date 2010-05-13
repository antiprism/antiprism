#!/usr/bin/python

# Antiprism Resource File - http://www.antiprism.com
# This file may be copied, modified and redistributed
#
# Make a helix by bonding face 'f0' of 'poly_file' to face 'f1' of a copy
# turning the face 'offset' vertices, and repeating 'length' times
# Adrian Rossiter <adrian@antiprism.com>

import os
import sys

if len(sys.argv) < 2:
   sys.stderr.write("usage: helix.py poly_file f0 f1 offset length\n")
   sys.exit(0);

poly = sys.argv[1]
f0 = int(sys.argv[2])

f1 = 0
if(len(sys.argv) > 3):
   f1 = int(sys.argv[3])

offset = 0
if(len(sys.argv) > 4):
   offset = int(sys.argv[4])

length = 10
if(len(sys.argv) > 5):
   length = int(sys.argv[5])

os.system("off_util %s > tmp_helix.off" % poly)
for i in range(length-1):
   os.system("off_align -o tmp_helix.off -b tmp_helix.off -F %d,%d,%d %s " % (f1, f0, offset, poly))
os.system("off_util tmp_helix.off")
   

   

