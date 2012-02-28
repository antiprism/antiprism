# Antiprism Example File - http://www.antiprism.com
# This file may be copied, modified and redistributed
#
# Display the demonstration models - Jean-Pierre Demailly

if test "$1" != ""
then
   GEOMETRY="-geometry $1"
else
   GEOMETRY="-geometry 600x500+20+20"
fi

off_color -f P -m map_orange:purple U75 -e white | antiview $GEOMETRY
geodesic -M p -f 8 ico | antiview -v b -x ef $GEOMETRY
poly_weave ico | off_util -M a | off_color -f U -e F | antiview -v 0.08 -x f $GEOMETRY
off_color -v 0.8,0.0,0.0 -f white -e grey30 geo_6 | antiview -v 0.04 -e 0.02 $GEOMETRY
string_art -l 20 -T 0,0,0.707 -R 0,0,45 -l 20 -T 0,0,-0.707 -R 0,0,-45 | off_color -e 1.0,1.0,0.7 -v 1.0,0.2,0.0 | poly_kscope -s O | antiview -v 0.03 -B black $GEOMETRY
zono -m v tr_ico | off_color -f 0.4,0.4,0.9 -v 0,1.0,1.0 | antiview -v 0.5 $GEOMETRY
zono -m v tr_ico |pol_recip| off_color -f P  -m map_0/0.5/0.5:darkred | antiview -v 0.5 -E white $GEOMETRY
unitile2d -s t -l 90 -w 30 11 | off_color -f N -v 1.0,0.0,0.0 -m spread+8  | antiview -v 0.03 -e 0.025 -E grey80 -t no_tri $GEOMETRY

