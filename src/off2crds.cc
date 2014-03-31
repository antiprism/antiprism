/*
   Copyright (c) 2006-2009, Adrian Rossiter

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

      The above copyright notice and this permission notice shall be included
      in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
  IN THE SOFTWARE.
*/

/*
   Name: off2crds.cc
   Description: extract coordinates from an OFF file 
   Project: Antiprism - http://www.antiprism.com
*/



#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"


using std::string;
using std::vector;


class o2c_opts: public prog_opts {
   public:
      string sep;
      int sig_digits;
      string ifile;
      string ofile;

      o2c_opts(): prog_opts("off2crds"), sep(" "), sig_digits(17) {}
      void process_command_line(int argc, char **argv);
      void usage();
};



void o2c_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Extract coordinates from a file in OFF format. If input_file is not given\n"
"then input is read from standard input.\n"
"\n"
"Options\n"
"%s"
"  -s <sep>  string to separate coordinates (default \" \")\n"
"  -d <dgts> number of significant digits (default 17) or if negative\n"
"            then the number of digits after the decimal point\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


void o2c_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hs:o:d:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 's':
            sep = optarg;
            break;
            
         case 'd':
            if(!read_int(optarg, &sig_digits, errmsg))
               error(errmsg, c);
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }
   
   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];
   
}


int main(int argc, char *argv[])
{
   o2c_opts opts;
   opts.process_command_line(argc, argv);

   col_geom_v geom;
   geom_read_or_error(geom, opts.ifile, opts);

   char errmsg[MSG_SZ];
   if(!geom.write_crds(opts.ofile, errmsg, opts.sep.c_str(), opts.sig_digits))
      opts.error(errmsg);

   return 0;
}


