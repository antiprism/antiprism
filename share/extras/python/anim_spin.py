#!/usr/bin/python

# Antiprism Resource File - http://www.antiprism.com
# This file may be copied, modified and redistributed
#
# Example animation of a spinning polyhedron
# Adrian Rossiter <adrian@antiprism.com>

import os
import math

# command line to render a file
povray_cmd = "povray +a +H100 +W100"                     # unix/linux povray
#povray_cmd = "pvengine /exit /render +a +H300 +W400"     # windows povray
#povray_cmd = "povray +a +H300 +W400"                     # windows megapov

length = 15

for i in range(0, length):
   t = float(i)/length
   rot_ang = 72*t
   os.system("polygon anti 5/3 | off_color -f P -m map_darkblue:ivory | off_trans -R 0,0,%g -R -60,0,0 -o tmp.off" % (rot_ang))
   os.system("off2pov -v 0.0 -E 0.5,0.5,0.7 -e 0.013 -B white -D 1.4 -C 0,-0.03,0 -o panim_%03d.pov tmp.off" % i)
   os.system("%s declare=AspectRatio=1 panim_%03d.pov" % (povray_cmd, i))


