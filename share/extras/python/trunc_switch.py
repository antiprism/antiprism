#!/usr/bin/python

# Antiprism Resource File - http://www.antiprism.com
# This file may be copied, modified and redistributed
#
# Transformation by truncation between a regular polyhedron and dual
# Adrian Rossiter <adrian@antiprism.com>

import os
import math

# command line to render a file
povray_cmd = "povray +a +H300 +W400"                     # unix/linux povray
#povray_cmd = "pvengine /exit /render +a +H300 +W400"     # windows povray
#povray_cmd = "povray +a +H300 +W400"                     # windows megapov

poly = 4                               # index of poly in polys list below
num_frames = 100                       # number of frames to render

polys = [
   [ "U1", 1, 1.2, "10,-15,180"],
   [ "U5", math.sqrt(2)/2, 1.6, "20,20,0"],
   [ "U23", (math.sqrt(5)+1)/2, 3, "20,15,0"],
   [ "U35", 1+(math.sqrt(5)+1)/2, 3, "20,15,0"],
   [ "U53", (math.sqrt(5)+1)/2, 1.8, "20,15,0"]
]

# set up a colour map file
os.system("echo \"3 = 1 0.5 0 0.7\" >  cmap.txt")
os.system("echo \"6 = 1 0.5 0 0.7\" >> cmap.txt")
os.system("echo \"4 = 0 1 0.5 0.7\" >> cmap.txt")
os.system("echo \"8 = 0 1 0.5 0.7\" >> cmap.txt")
os.system("echo \"5 = 0 0.5 1 0.7\" >> cmap.txt")
os.system("echo \"10 = 0 0.5 1 0.7\" >> cmap.txt")

for i in range(0, num_frames):
   os.system("unipoly %s | off_trans -s e > t1.off" % polys[poly][0])
   div_rat = polys[poly][1]
   t = (1+div_rat) * float(i)/num_frames
   if(t<1):
      rat = math.fabs(0.5 - t)
   else:
      rat = math.fabs((1 + 0.5*div_rat - t)/div_rat)
      os.system("pol_recip -c 0,0,0 -r E -o t1.off t1.off")
   os.system("off_util -T %g t1.off | off_color -f n -m cmap.txt > tmp.off" %(rat))
   os.system("off2pov -B black -E white -F 0.0,0.5,1,0.7 -S 2 -v 0.0125 -e 0.01 -D %g -C 0,0,0 -R %s -o panim_%03d.pov tmp.off" % (polys[poly][2], polys[poly][3], i))
   os.system("%s panim_%03d.pov" % (povray_cmd, i) )

   

