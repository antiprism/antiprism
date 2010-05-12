#!/usr/bin/python

# Antiprism Resource File - http://www.antiprism.com
# This file may be copied, modified and redistributed
#
# Idealised truncation between +infinity and -infinity
# Adrian Rossiter <adrian@antiprism.com>

import os
import math

# command line to render a file
povray_cmd = "povray +a +H300 +W400"                     # unix/linux povray
#povray_cmd = "pvengine /exit /render +a +H300 +W400"     # windows povray
#povray_cmd = "povray +a +H300 +W400"                     # windows megapov

num_frames = 100                       # number of frames to render
poly = "U53"                           # uniform poly to truncate

# set up a colour map file
os.system("echo \"3 = 1 0.5 0 0.7\" >  cmap.txt")
os.system("echo \"6 = 1 0.5 0 0.7\" >> cmap.txt")
os.system("echo \"4 = 0 1 0.5 0.7\" >> cmap.txt")
os.system("echo \"8 = 0 1 0.5 0.7\" >> cmap.txt")
os.system("echo \"5 = 0 0.5 1 0.7\" >> cmap.txt")
os.system("echo \"10 = 0 0.5 1 0.7\" >> cmap.txt")


for i in range(0, num_frames):
   t = float(i)/num_frames;
   third = 1.0/3
   if(t<third):
      a = (third-t)/third;
      fact = -math.pow(2, 10*a*a)+1
   elif(t<(2*third)):
      fact = (t-third)/third
   else:
      a = (t-2*third)/third;
      fact = math.pow(2, 10*a*a)
      
   os.system("unipoly %s | off_util -T %g | off_util -S | off_color -f n -m cmap.txt > tmp.off" % (poly, fact))
   os.system("off2pov -S 2 -B black -E white -v 0.017 -e 0.015 -D 2.4 -C 0,0,0 -R -20,30,0 -o panim_%03d.pov tmp.off" % i)
   os.system("%s panim_%03d.pov" % (povray_cmd, i))
   #print "third=%g, a=%g, t=%g, fact=%g\n" %(third, a, t, fact)

   

