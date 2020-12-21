/*
   Copyright (c) 2007-2020, Roger Kaufman
   Includes ideas and algorithms by George W. Hart, http://www.georgehart.com

   Antiprism - http://www.antiprism.com

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
   Name: conway.cc
   Description: Conway Notation
                Implementation of George Hart's Conway Notation
                http://www.georgehart.com/virtual-polyhedra/conway_notation.html
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

using namespace anti;

#define CN_ONE_THIRD 1 / 3.0
#define CN_ONE_HALF 0.5

struct ConwayOperator {
  string operator_short;
  string operator_name;
  int sub1;
  int sub2;
  int sides;
  bool hart_operator;
};

// clang-format off
ConwayOperator conway_operator_list[]{
    {"a",  "ambo",            -1, -1, -1, true  },
    {"B",  "bowtie",          -1, -1, -1, false },
    {"b",  "bevel",            0, -1, -1, false }, // subscript >= 0
    {"c",  "chamfer",         -1, -1, -1, true  },
    {"d",  "dual",            -1, -1, -1, true  },
    {"E",  "ethyl",           -1, -1, -1, false },
    {"e",  "expand",           0, -1, -1, false }, // subscript >= 0
    {"G",  "opposite-lace",   -1, -1, -1, false },
    {"g",  "gyro",             1, -1, -1, true  }, // subscript >= 1
    {"J",  "joined-medial",   -1, -1, -1, false }, // replaces wiki M0
    {"j",  "join",            -1, -1, -1, false },
    {"K",  "stake",           -1, -1,  3, false }, // face sides >= 3
    {"k",  "kis",             -1, -1,  3, true  }, // vertex sides >= 3
    {"L",  "lace",             0, -1,  3, false }, // subscript >= 0, face sides >= 3
    {"l",  "loft",             0, -1,  3, false }, // subscript >= 0, face sides >= 3
    {"M",  "medial",           0, -1, -1, false }, // subscript >= 0
    {"m",  "meta",             0, -1, -1, false }, // subscript >= 0
    {"n",  "needle",          -1, -1, -1, false },
    {"o",  "ortho",            0, -1, -1, false }, // subscript >= 0
    {"p",  "propeller",       -1, -1, -1, true  },
    {"q",  "quinto",          -1, -1, -1, false },
    {"r",  "reflect",         -1, -1, -1, false },
    {"S",  "seed",            -1, -1, -1, false },
    {"s",  "snub",             1, -1, -1, false }, // subscript >= 1
    {"t",  "truncate",        -1, -1,  2, false }, // vertex sides >= 2
    {"u",  "subdivide",        1,  1, -1, false }, // two subscripts >= 1
    {"W",  "waffle",          -1, -1, -1, false },
    {"w",  "whirl",           -1, -1, -1, true  },
    {"X",  "cross",           -1, -1, -1, false },
    {"z",  "zip",             -1, -1, -1, false },
    {"+",  "orient positive", -1, -1, -1, false },
    {"-",  "orient negative", -1, -1, -1, false },
};

struct ConwaySeed {
  string seed;
  string seed_name;
  int sides;
};

ConwaySeed conway_seed_list[]{
    {"T",  "tetrahedron",     -1 },
    {"C",  "cube",            -1 },
    {"O",  "octahedron",      -1 },
    {"I",  "icosahedron",     -1 },
    {"D",  "dodecahedron",    -1 },
    {"P",  "prism",            3 }, // sides >= 3 required
    {"A",  "antiprism",        3 }, // sides >= 3 required
    {"Y",  "pyramid",          3 }, // sides >= 3 required
    {"Z",  "polygon",          3 }, // sides >= 3 required
    {"R",  "random",           3 }, // sides >= 3 required
};
// clang-format on

bool find_operator(const char &op)
{
  bool found = false;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short[0]) {
      found = true;
      break;
    }
  }
  return found;
}

int find_sub1_allowed(const char &op)
{
  int found = -1;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short[0]) {
      if (conway_operator_list[i].sub1 > -1)
        found = conway_operator_list[i].sub1;
      break;
    }
  }
  return found;
}

int find_sub2_allowed(const char &op)
{
  int found = -1;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short[0]) {
      if (conway_operator_list[i].sub2 > -1)
        found = conway_operator_list[i].sub2;
      break;
    }
  }
  return found;
}

int find_sides_allowed(const char &op)
{
  int found = -1;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short[0]) {
      if (conway_operator_list[i].sides > -1)
        found = conway_operator_list[i].sides;
      break;
    }
  }
  return found;
}

bool find_seed(const char &seed)
{
  bool found = false;
  unsigned int last_op = sizeof(conway_seed_list) / sizeof(conway_seed_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (seed == conway_seed_list[i].seed[0]) {
      found = true;
      break;
    }
  }
  return found;
}

int find_seed_sides(const char &seed)
{
  int found = -1;
  unsigned int last_op = sizeof(conway_seed_list) / sizeof(conway_seed_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (seed == conway_seed_list[i].seed[0]) {
      if (conway_seed_list[i].sides > -1)
        found = conway_seed_list[i].sides;
      break;
    }
  }
  return found;
}

class ops {
public:
  int op_pos;
  char op;
  int sub1;
  int sub2;
  int sides;
  ops(int n, char o, int s1, int s2, int si)
      : op_pos(n), op(o), sub1(s1), sub2(s2), sides(si)
  {
  }
};

bool cmp_ops(const ops *a, const ops *b) { return a->op_pos > b->op_pos; }

// op, sub and sides will be reset
void write_operation(vector<ops *> &operations, int &op_count, char &op,
                     int &sub1, int &sub2, int &sides)
{
  operations.push_back(new ops(op_count++, op, sub1, sub2, sides));
  op = '\0';
  sub1 = 1;
  sub2 = 1;
  sides = -1;
}

Status validate_cn_string(const string &cn_string, vector<ops *> &operations,
                          char &seed, int &seed_size, const string alpha_user)
{
  seed = '\0';
  seed_size = 0;
  int sub1 = 1;
  int sub2 = 1;
  int sides = -1;

  Status stat;
  bool end_of_command = false;

  int require_seed_sides = -1;
  int possible_sub = -1;
  int possible_sides = -1;
  int possible_repeats = -1;

  char super_symbol = '^';
  char sub_symbol = '_';
  char sides_symbol = ':';
  bool designated_sub = false;
  bool designated_sides = false;

  string digits = "0123456789";

  int num_val1 = -1;
  int num_val2 = -1;

  char pending_op = '\0';
  int pending_pos = 0;

  int op_count = 0;
  for (unsigned int i = 0; i < cn_string.length(); i++) {
    char current = cn_string[i];

    bool is_digit = (digits.find(current) != string::npos);

    bool is_delimiter = ((current == sub_symbol) || (current == sides_symbol));

    // seed with no N requirement has been found, no more characters expected
    if (end_of_command)
      return stat.set_error(
          msg_str("extra characters past seed \'%c\': \'%s\' at position %d",
                  seed, cn_string.substr(i).c_str(), i + 1));

    // seed with N requirement has been found, digits expected
    if ((require_seed_sides > -1) && !is_digit)
      return stat.set_error(msg_str(
          "seed \'%c\' at position %d requires a number of sides of %d or more",
          seed, i + 1, require_seed_sides));

    // a number found where it shouldn't be
    if (is_digit) {
      if (i == 0)
        return stat.set_error(msg_str(
            "digit \'%c\' unexpected in first position %d", current, i + 1));

      if (((possible_sub == -1) && (possible_sides == -1) &&
           (require_seed_sides == -1) && (possible_repeats == -1)))
        return stat.set_error(
            msg_str("operator \'%c\' at position %d unexpected number value "
                    "specified: %c",
                    pending_op, i + 1, current));
    }

    // a number must follow a superscript
    if ((possible_repeats == 1) && !is_digit)
      return stat.set_error(msg_str(
          "superscript \'%c\' at position %d requires a number of %d or more",
          super_symbol, i + 1, possible_repeats));

    // if no more digits or delimiter, handle delayed write with current values,
    // process next character
    if (((possible_sub > -1) || (possible_sides > -1)) && !is_digit &&
        !is_delimiter) {
      write_operation(operations, op_count, pending_op, sub1, sub2, sides);
      possible_sub = -1;
      possible_sides = -1;
      num_val1 = -1;
      num_val2 = -1;
    }

    // its a sides symbol
    if (current == sides_symbol) {
      // PATCH: L0 does not have sides value
      if (pending_op == 'L' && sub1 == 0)
        return stat.set_error(
            msg_str("sides symbol \'%c\' at position %d not allowed for L0",
                    sides_symbol, i + 1));

      // sides not already set, or set and encountered again
      if (possible_sides == -1 || sides != -1) {
        string sides_str;
        if (sides != -1)
          sides_str = msg_str("(sides already set: %d)", sides);
        return stat.set_error(msg_str("sides symbol \'%c\' at position %d for "
                                      "operator \'%c\' unexpected %s",
                                      sides_symbol, i + 1, pending_op,
                                      sides_str.c_str()));
      }
      designated_sides = true;
    }
    // its a subscript
    else if (current == sub_symbol) {
      if (possible_sub == -1)
        return stat.set_error(
            msg_str("subscript symbol \'%c\' at position %d unexpected",
                    sub_symbol, i + 1));
      if (sides != -1)
        return stat.set_error(
            msg_str("subscript symbol \'%c\' at position %d found after sides",
                    sub_symbol, i + 1));
      designated_sub = true;
    }
    // its a superscript
    else if (current == super_symbol) {
      pending_op = current;
      pending_pos = i + 1;

      if (i == 0)
        return stat.set_error(
            msg_str("superscript \'%c\' unexpected in first position %d",
                    pending_op, i + 1));
      if (designated_sub || designated_sides)
        return stat.set_error(msg_str(
            "superscript \'%c\' unexpected in position %d, digit expected",
            pending_op, i + 1));

      possible_repeats = 1;
    }
    // its a digit
    else if (is_digit) {
      // find full number string and advance counter
      int digits_start = i;
      int digits_end = i;
      while ((digits_end + 1 < (int)cn_string.length()) &&
             (digits.find(cn_string[digits_end + 1]) != string::npos))
        digits_end++;
      // set iterator to one before next character to test
      i = digits_end;

      int num_val = std::stoi(
          cn_string.substr(digits_start, (digits_end - digits_start + 1)));

      // if sides symbol was seen, we need to use num_val1
      if (num_val1 != -1 && designated_sides)
        num_val1 = -1;

      if (num_val1 == -1)
        num_val1 = num_val;
      // check for extra subscripts here
      else if (num_val2 == -1) {
        int min_allowed_sub2 = find_sub2_allowed(pending_op);
        if (min_allowed_sub2 == -1)
          return stat.set_error(msg_str("operator \'%c\' at position %d, "
                                        "second subscript of %d unexpected",
                                        pending_op, pending_pos, num_val));
        else {
          if (num_val < min_allowed_sub2)
            return stat.set_error(msg_str(
                "operator \'%c\' at position %d, subscript must be %d or more",
                pending_op, pending_pos, num_val));
          else
            num_val2 = num_val;
        }
      }
      else
        return stat.set_error(msg_str("too many subscripts for operator \'%c\' "
                                      "at position %d, value %d unexpected",
                                      pending_op, pending_pos, num_val));

      // its a pending superscript
      if (possible_repeats == 1) {
        if (num_val1 < 1)
          return stat.set_error(
              msg_str("superscript \'%c\' at position %d requires a number "
                      "of %d or more",
                      pending_op, pending_pos, possible_repeats));

        // repeat last operator
        char op_last = operations[operations.size() - 1]->op;
        int sub1_last = operations[operations.size() - 1]->sub1;
        int sub2_last = operations[operations.size() - 1]->sub2;
        int sides_last = operations[operations.size() - 1]->sides;
        for (int i = 0; i < num_val1 - 1; i++)
          operations.push_back(
              new ops(op_count++, op_last, sub1_last, sub2_last, sides_last));

        possible_repeats = -1;
      }
      // its a pending seed
      else if (require_seed_sides > -1) {
        seed_size = num_val1;
        if (seed_size < require_seed_sides)
          return stat.set_error(
              msg_str("seed \'%c\' at position %d requires a number of "
                      "sides of %d or more, got %d",
                      seed, pending_pos, require_seed_sides, num_val1));
        require_seed_sides = -1;
        end_of_command = true;
      }
      // its a designated subscript, check before passive
      else if (designated_sub) {
        sub1 = num_val1;
        if (num_val2 != -1)
          sub2 = num_val2;
        if (sub1 < possible_sub)
          return stat.set_error(msg_str(
              "operator \'%c\' at position %d, subscript must be %d or more",
              pending_op, pending_pos, possible_sub));
        designated_sub = false;
      }
      // its a designated sides, check before passive
      else if (designated_sides) {
        sides = num_val1;
        if (sides < possible_sides)
          return stat.set_error(
              msg_str("operator \'%c\' at position %d requires a number of "
                      "sides of %d or more, got %d",
                      pending_op, pending_pos, possible_sides, num_val1));
        designated_sides = false;
      }
      // its a pending sides, passive
      // checking sides first before subs gives sides precedence
      else if (possible_sides > -1) {
        sides = num_val1;
        if (sides < possible_sides)
          return stat.set_error(
              msg_str("operator \'%c\' at position %d requires a number of "
                      "sides of %d or more, got %d",
                      pending_op, pending_pos, possible_sides, num_val1));
      }
      // its a pending subscript, passive
      else if (possible_sub > -1) {
        sub1 = num_val1;
        if (num_val2 != -1)
          sub2 = num_val2;
        if (sub1 < possible_sub)
          return stat.set_error(msg_str("operator \'%c\' at position %d "
                                        "subscript must be %d or more, got %d",
                                        pending_op, pending_pos, possible_sub,
                                        num_val1));
      }
    }
    // its an operator
    else if (find_operator(current)) {
      pending_op = current;
      pending_pos = i + 1;

      possible_sub = find_sub1_allowed(current);
      possible_sides = find_sides_allowed(current);

      // if an operator needs no value number, write immediately
      if (possible_sub == -1 && possible_sides == -1)
        write_operation(operations, op_count, current, sub1, sub2, sides);
    }
    // its an alpha, needs no value number, write immediately
    else if (alpha_user.find(current) != string::npos) {
      pending_op = current; // for reporting
      pending_pos = i + 1;

      write_operation(operations, op_count, current, sub1, sub2, sides);
    }
    // its a seed
    else if (find_seed(current)) {
      seed = current;
      pending_pos = i + 1;
      require_seed_sides = find_seed_sides(seed);
      if (require_seed_sides == -1)
        end_of_command = true;
    }
    // fell through, unused character
    else
      return stat.set_error(msg_str(
          "unexpected character \'%c\' at position %d", current, i + 1));
  }

  // if superscript was last character
  if ((pending_op == super_symbol) && (possible_repeats == 1))
    return stat.set_error(msg_str(
        "superscript \'%c\' at position %d requires a number of %d or more",
        super_symbol, pending_pos, possible_repeats));

  // if this happened its probably a dangling seed that needs N
  if (require_seed_sides > -1) {
    return stat.set_error(msg_str(
        "seed \'%c\' at position %d requires a number of sides of %d or more",
        seed, pending_pos, require_seed_sides));
  }

  // possible pending passive operation if no seed was specified
  if ((possible_sub > -1) || (possible_sides > -1)) {
    write_operation(operations, op_count, pending_op, sub1, sub2, sides);
  }

  return Status::ok();
}

// G. Hart Commentary
// P4 --> C    (C is prism)
// A3 --> O    (O is antiprism)
// Y3 --> T    (T is pyramid)
// e  --> aa   (abbr. for explode)
// b  --> ta   (abbr. for bevel)
// o  --> jj   (abbr. for ortho)
// m  --> kj   (abbr. for meta)
// t(n) --> dk(n)d  (dual operations) (see special case)
// j  --> dad  (dual operations)
// s  --> dgd  (dual operations)
// dd --> null (order 2)
// ad --> a    (a_ = ad_)
// gd --> g    (g_ = gd_)
// aY --> A    (interesting fact)
// dT --> T    (self-dual)
// gT --> D    (symm change)
// aT --> O    (symm change)
// dC --> O    (dual pair)
// dO --> C    (dual pair)
// dI --> D    (dual pair)
// dD --> I    (dual pair)
// aO --> aC   (for uniqueness)
// aI --> aD   (for uniqueness)
// gO --> gC   (for uniqueness)
// gI --> gD   (for uniqueness)

struct ResolveItem {
  const char *target;
  const char *resolve;
};

// clang-format off
// order sensitive, do not sort
ResolveItem resolve_item_list[] = {
    {  "P4", "C",   },
    {  "A3", "O",   },
    {  "Y3", "T",   },
    {  "e",  "aa",  },
    {  "b",  "ta",  },
    {  "o",  "jj",  },
    {  "m",  "kj",  },
    {  "t",  "dk",  },
    {  "j",  "dad", },
    {  "s",  "dgd", },
    {  "dd", "",    },
    {  "ad", "a",   },
    {  "gd", "g",   },
    {  "aY", "A",   },
    {  "dT", "T",   },
    {  "gT", "D",   },
    {  "aT", "O",   },
    {  "dC", "O",   },
    {  "dO", "C",   },
    {  "dI", "D",   },
    {  "dD", "I",   },
    {  "aO", "aC",  },
    {  "aI", "aD",  },
    {  "gO", "gC",  },
    {  "gI", "gD",  },
};
// clang-format on

string resolved_cn_string(const string &cn_string)
{
  string resolve_string = cn_string;

  int num_subst = sizeof(resolve_item_list) / sizeof(resolve_item_list[0]);

  string target;
  string resolve;

  string digits = "0123456789";

  // first 3 targets are positional
  if (resolve_string.length() > 1) {
    string tmp = resolve_string.substr(resolve_string.length() - 2, 2);
    for (unsigned int i = 0; i < 3; i++) {
      target = resolve_item_list[i].target;
      resolve = resolve_item_list[i].resolve;
      if (tmp == target)
        resolve_string.replace(resolve_string.length() - 2, target.length(),
                               resolve);
    }
  }

  // 4 to 6
  for (unsigned int i = 3; i < 7; i++) {
    target = resolve_item_list[i].target;
    resolve = resolve_item_list[i].resolve;
    for (size_t pos = resolve_string.find_first_of(target, 0);
         pos != string::npos; pos = resolve_string.find_first_of(target, 0)) {
      // RK - have to check now since some operators can have a number
      if ((pos + 1 == resolve_string.length()) ||
          (digits.find(resolve_string[pos + 1]) == string::npos))
        resolve_string.replace(pos, target.length(), resolve);
      else {
        // RK - if so e, b, o and m are temporarily replaced with @, #, $, %
        // so that loop can continue
        if (resolve_string[pos] == 'e')
          resolve_string.replace(pos, 1, "@");
        else if (resolve_string[pos] == 'b')
          resolve_string.replace(pos, 1, "#");
        else if (resolve_string[pos] == 'o')
          resolve_string.replace(pos, 1, "$");
        else if (resolve_string[pos] == 'm')
          resolve_string.replace(pos, 1, "%");
      }
    }
  }

  // if temporary characters exists
  replace(resolve_string.begin(), resolve_string.end(), '@', 'e');
  replace(resolve_string.begin(), resolve_string.end(), '#', 'b');
  replace(resolve_string.begin(), resolve_string.end(), '$', 'o');
  replace(resolve_string.begin(), resolve_string.end(), '%', 'm');

  // 7 is special case
  target = resolve_item_list[7].target;
  resolve = resolve_item_list[7].resolve;
  for (size_t pos = resolve_string.find_first_of(target, 0);
       pos != string::npos; pos = resolve_string.find_first_of(target, 0)) {
    resolve_string.replace(pos, target.length(), resolve);
    size_t pos2 = resolve_string.find_first_not_of(digits, pos + 2);
    if (pos2 != string::npos)
      resolve_string.insert(pos2, "d");
    else
      resolve_string.append("d");
  }

  // 8 to end. Because of length greater than 1 of the target, string.find and
  // replace cannot be used
  for (int i = 8; i < num_subst; i++) {
    target = resolve_item_list[i].target;
    resolve = resolve_item_list[i].resolve;
    int j = 0;
    int stop = resolve_string.length() - target.length() + 1;
    while (j <= stop) {
      string tmp = resolve_string.substr(j, target.length());
      if (tmp == target) {
        resolve_string.replace(j, target.length(), resolve);
        j = -1; // we have to begin looking from the beginning once again, ok if
                // we modify stop
        stop = resolve_string.length() - target.length() + 1;
      }
      j++;
    }
  }

  if (resolve_string != cn_string) {
    fprintf(stderr, "%s resolved to: ", cn_string.c_str());
    if (resolve_string.length() == 0)
      fprintf(stderr, "%s", "NOTHING");
    else
      fprintf(stderr, "%s", resolve_string.c_str());
    fprintf(stderr, "\n");
  }

  return resolve_string;
}

class cn_opts : public ProgramOpts {
public:
  IterationControl it_ctrl;

  string ifile;
  string ofile;

  string cn_string;
  bool resolve_ops;
  bool hart_mode;
  bool tile_mode;
  bool reverse_ops;
  char seed;
  int seed_size;
  char planarize_method;
  bool unitize;
  bool verbosity;
  char face_coloring_method;
  int face_opacity;
  string face_pattern;
  int seed_coloring_method;

  double epsilon;

  Color vert_col;
  Color edge_col;

  TilingColoring col_type;

  ColorMapMulti map;

  vector<ops *> operations;

  // for on the fly user operators
  string alpha_user;
  std::map<char, vector<ops *>> operations_user;

  cn_opts()
      : ProgramOpts("conway"), cn_string(""), resolve_ops(false),
        hart_mode(false), tile_mode(false), reverse_ops(false), seed('\0'),
        seed_size(0), planarize_method('q'), unitize(false), verbosity(false),
        face_coloring_method('n'), face_opacity(-1), face_pattern("1"),
        seed_coloring_method(1), epsilon(0),
        vert_col(Color(255, 215, 0)),  // gold
        edge_col(Color(211, 211, 211)) // lightgray
  {
    it_ctrl.set_max_iters(1000);
    it_ctrl.set_status_checks("-1,1");
    it_ctrl.set_sig_digits(int(-log(::epsilon) / log(10) + 0.5));
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void extended_help()
{
  fprintf(stdout, R"(
Conway Notation was described by Mathematician John Conway to George Hart in
the late 1990's for a book they planned to coauthor. Due to an illness the book
never came to fruition and John Conway did not think there was enough a for a
separate publication. Conway gave George Hart permission to present it in
"Sculpture Based on Propellerized Polyhedra in the Proceedings of MOSAIC 2000"
The paper can be viewed here: http://www.georgehart.com/propello/propello.html
The project was expected to encourage more operations to be developed which has
happened in various places including here at Antiprism. (www.antiprism.com)

The following is a description of Conway Notation edited from the Conway
Notation web page by George W. Hart (http://www.georgehart.com)

More detailed information and examples can be found at
http://www.georgehart.com/virtual-polyhedra/conway_notation.html
and at
http://en.wikipedia.org/wiki/Conway_polyhedron_notation
and at
http://antitile.readthedocs.io/en/latest/conway.html

Basics: In this notation, one specifies a "seed" polyhedron with a capital
letter. Operations to perform on any polyhedron are specified with lower-case
letters preceding it. This program contains a small set of seeds and operators
from which an infinite number of derived polyhedra can be generated.

Note: This C++ port of Conway Notation can also operate on OFF files from
standard input if the seed polyhedron is not specified. (Antiprism Extension)

Seeds: The platonic solids are denoted T, O, C, I, and D, according to their
first letter. Other polyhedra which are implemented here include prisms, Pn,
antiprisms, An, and pyramids, Yn, where n is a number (3 or greater) which you
specify to indicate the size of the base you want, e.g., Y3=T, P4=C, and A3=O.

Operations: Currently, abdegjkmoprst are defined. They are motivated by the
operations needed to create the Archimedean solids and their duals from the
platonic solids.  Try each on a cube:

(Antiprism Extension: note that more operations have since been defined)

a = ambo   The ambo operation can be thought of as truncating to the edge
midpoints.  It produces a polyhedron, aX, with one vertex for each edge of X.
There is one face for each face of X and one face for each vertex of X.
Notice that for any X, the vertices of aX are all 4-fold, and that aX=adX.
If two mutually dual polyhedra are in "dual position", with all edges tangent
to a common sphere, the ambo of either is their intersection.  For example
aC=aO is the cuboctahedron.
Note: ambo is also known as "rectifying" the polyhedron, or rectification

b = bevel  The bevel operation can be defined by bX=taX.  bC is the truncated
cuboctahedron.  (Antiprism Extension: or "bn" where n is 0 or greater)
Note: bevel is also known as "omnitruncating" the polyhedron, or omnitruncation

d = dual   The dual of a polyhedron has a vertex for each face, and a face for
each vertex, of the original polyhedron, e.g., dC=O.  Duality is an operation
of order two, meaning for any polyhedron X, ddX=X, e.g., ddC=dO=C. 

e = expand This is Mrs. Stott's expansion operation.  Each face of X is
separated from all its neighbors and reconnected with a new 4-sided face,
corresponding to an edge of X.  An n-gon is then added to connect the 4-sided
faces at each n-fold vertex.  For example, eC is the rhombicuboctahedron.  It
turns out that eX=aaX and so eX=edX (Antiprism: "en" where n is 0 or greater)
Note: expand is also known as "cantellating" the polyhedron, or cantellation

g = gyro   The dual operation to s is g. gX=dsdX=dsX, with all 5-sided faces.
The gyrocube, gC=gO="pentagonal icositetrahedron", is dual to the snub cube.
g is like k but with the new edges connecting the face centers to the 1/3
points on the edges rather than the vertices. (Antiprism Extension: or "gn"
where n is 1 or greater)

j = join   The join operator is dual to ambo, so jX=dadX=daX.  jX is like kX
without the original edges of X.  It produces a polyhedron with one 4-sided
face for each edge of X.  For example, jC=jO is the rhombic dodecahedron.

k = kis    All faces are processed or kn = just n-sided faces are processed
The kis operation divides each n-sided face into n triangles.  A new vertex is
added in the center of each face, e.g., the kiscube, kC, has 24 triangular
faces.  The k operator is dual to t, meaning kX=dtdX.

m = meta   Dual to b, mX=dbX=kjX.  mC has 48 triangular faces.  m is like k
and o combined; new edges connect new vertices at the face centers to the old
vertices and new vertices at the edge midpoints.  mX=mdX.  mC is the
"hexakis octahedron".  (Antiprism Extension: or "mn" where n is 0 or greater)

o = ortho  Dual to e, oX=deX=jjX.  oC is the trapezoidal icositetrahedron, with
24 kite-shaped faces.  oX has the effect of putting new vertices in the middle
of each face of X and connecting them, with new edges, to the edge midpoints of
X.  (Antiprism Extension: or "on" where n is 0 or greater)

p = propeller    Makes each n-gon face into a "propeller" of an n-gon
surrounded by n quadrilaterals, e.g., pT is the tetrahedrally stellated
icosahedron. Try pkD and pt6kT. p is a self-dual operation, i.e., dpdX=pX and
dpX=pdX, and p also commutes with a and j, i.e., paX=apX. (This and the next
are extensions were added by George Hart and not specified by Conway)

r = reflect   Changes a left-handed solid to right handed, or vice versa, but
has no effect on a reflexible solid. So rC=C, but compare sC and rsC.

s = snub   The snub operation produces the snub cube, sC, from C.  It can be
thought of as eC followed by the operation of slicing each of the new 4-fold
faces along a diagonal into two triangles.  With a consistent handedness to
these cuts, all the vertices of sX are 5-fold.  Note that sX=sdX.
(Antiprism Extension: or "sn" where n is 1 or greater)

t = truncate  All faces are processed or tn = just n-sided faces are processed
Truncating a polyhedron cuts off each vertex, producing a new n-sided face for
each n-fold vertex.  The faces of the original polyhedron still appear, but
have twice as many sides, e.g., the tC has six octagonal sides corresponding to
the six squares of the C, and eight triangles corresponding to the cube's eight
vertices.


Antiprism Extension: Further operations added. Also see
http://en.wikipedia.org/wiki/Conway_polyhedron_notation

c = chamfer   New hexagonal faces are added in place of edges

B = bowtie    Bowtie like triangles divide pentagonal faces

E = ethyl     like expand but triangles are divided into 3 kites

G = opposite-lace Similar to lace, except triangles split opposite lace

J = joined-medial Like medial, but new rhombic faces in place of original edges

K = stake     Subdivide faces with central quads, and triangles
              All faces processed or can be "Km" where m is 3 or greater

L0 = joined-lace  Similar to lace, except with new quad faces across original
                  edges

L = lace      An augmentation of each face by an antiprism, adding a twist
              smaller copy of each face, and triangles between
              Subscript as "Ln" where n is 0 or greater
              All faces processed or can be "L:m" where m is 3 or greater
              Both may be specified as "L_n:m"

l = loft      An augmentation of each face by prism, adding a smaller copy of
              each face with trapezoids between the inner and outer ones
              Subscript as "ln" where n is 0 or greater
              All faces processed or can be "l:m" where m is 3 or greater
              Both may be specified as "l_n:m"

M = medial    Similar to meta except no diagonal edges added, creating quad
              faces. All faces processed or can be "Mm" where m is 0 or greater

n = needle    Dual of truncation, triangulate with 2 triangles across every
              edge. This bisect faces across all vertices and edges, while
              removing original edges

q = quinto    ortho followed by truncation of vertices centered on original
              faces. This create 2 new pentagons for every original edge

S = seed      Seed form

u = subdivide Ambo while retaining original vertices. Similar to Loop
              subdivision surface for triangle face
              One subscript as "un" where n is 1 or greater
              Two subscript as "un_m" or "u_n_m" where n and m are 1 or greater

W = waffle    Truncation on all vertices and then all faces split into sections

w = whirl     Gyro followed by truncation of vertices centered on original
              faces. This create 2 new hexagons for every original edge

X = cross     Combination of kis and subdivide operation. Original edges are
              divided in half, with triangle and quad faces

z = zip       Dual of kis or truncation of the dual. This create new edges
              perpendicular to original edges, a truncation beyond "ambo" with
              new edges "zipped" between original faces. It is also called
              bitruncation

Orientation of the input model will have an effect on chiral operations such as
snub or whirl. The orientation mode is set to positive by default. Operations
have been added to control orientation mode. The mode will remain until changed
+ (plus sign) = positive orientation  - (minus sign) = negative orientation
Changing orientation mode can be placed anywhere in the operation string

Summary of operators which can take n as subscript or r as face/vertex number

b  - n may be 0 or greater (default: 1)
e  - n may be 0 or greater (default: 1)
g  - n may be 1 or greater (default: 1)
K  - r may be 3 or greater representing face sides
k  - r may be 3 or greater representing face sides
L  - n,r n may be 0 or greater, r may be 3 or greater
l  - n,r n may be 0 or greater, r may be 3 or greater
       both n amd r may be used togeter as L_n:r or l_n:r (L_0: may not have r)
       without delimiters Lr and lr, r is face sides and subscript default to 1
M  - n may be 0 or greater (default: 1)
m  - n may be 0 or greater (default: 1)
o  - n may be 0 or greater (default: 1)
s  - n may be 1 or greater (default: 1)
t  - r may be 2 or greater representing vertex connections (2 in tiles)
u  - n,m n and m may be 1 or greater (default: _1_1);

Antiprism Extension: any operation can be repeated N time by following it with
the superscript symbol ^ and a number greater than 0. Examples: a^3C M0^2T

Seeds which require a number n, 3 or greater

P  - Prism
A  - Antiprism
Y  - Pyramid
Z  - Polygon (Antiprism Extension)
R  - Random Convex Polyhedron (Antiprism Extension)

Note: Antiprism Extensions will work on tilings. Hart algorithms (-d) will not
e.g.: unitile2d 3 | conway p -t | antiview -v 0.1 (-t for tile mode)

Regular 2D tilings can be constructed from base polygons. The basic tilings are

            One Layer  Two Layers  Three Layers...
Square:     oZ4        o2Z4        o3Z4
Hexagonal:  tkZ6       ctkZ6       cctkZ6
Triangular: ktkZ6      kctkZ6      kcctkZ6 (kis operation on Hexagonal)

Name                   Vertex Fig  Op     String Dual Name              String
Square                 4,4,4,4            oZ4    Square                 do2Z4
Truncated Square       4,8,8       trunc  toZ4   Tetrakis Square        dto2Z4
Snub Square            3,3,4,3,4   snub   soZ4   Cairo Pentagonal       dso2Z4
Triangular             3,3,3,3,3,3 kis    ktkZ6  Hexagonal              ddctkZ6
Hexagonal              6,6,6              tkZ6   Triangular             dctkZ6
Trihexagonal           3,6,3,6     ambo   atkZ6  Rhombille              dactkZ6
Snub Trihexagonal      3,3,3,3,6   snub   stkZ6  Floret Pentagonal      dsctkZ6
Truncated Hexagonal    3,12,12     trunc  ttkZ6  Triakis triangular     dtctkZ6
Rhombitrihexagonal     3,4,6,4     expand etkZ6  Deltoidal Trihexagonal dectkZ6
Truncated Trihexagonal 4,6,12      bevel  btkZ6  Kisrhombille           dbctkZ6
Elongated Triangular   3,3,3,4,4   NonWythoffian Prismatic Triangular   none


Substitutions used by George Hart algorithms

)");

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j <= 20; j += 5) {
      int k = i + j;
      fprintf(stdout, "%-2s -> %-8s", resolve_item_list[k].target,
              resolve_item_list[k].resolve);
    }
    fprintf(stdout, "\n");
  }

  fprintf(stdout, R"(
Equivalent Operations

b0 = z        e0 = d        o0 = S        m0 = k        M0 = o
b1 = b        e1 = e        o1 = o        m1 = m        M1 = M
)");
}

/*
Various equivalent forms

jT    = C
sT    = I
dA3   = C
k5A5  = I (A special gyroelongated dipyramid)
t5dA5 = D (A special truncated trapezohedron)
t4daC = cC
t4kC  = lC
daC   = jC
t5daaD or t5deD or t5oD = qD
dedD   = oD
*/

void cn_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [Conway Notation string] [input_file]

Conway Notation uses algorithms by George W. Hart (http://www.georgehart.com)
http://www.georgehart.com/virtual-polyhedra/conway_notation.html

Antiprism Extensions: Further operations added. See
http://en.wikipedia.org/wiki/Conway_polyhedron_notation
and
http://antitile.readthedocs.io/en/latest/conway.html

Read a polyhedron from a file in OFF format.
If input_file is not given and no seed polyhedron is given in the notation
string then the program reads from standard input.

Options
%s
  -H        Conway Notation detailed help. seeds and operator descriptions
  -s        apply Conway Notation string substitutions
  -g        use George Hart algorithms (sets -s)
  -c <op=s> user defined operation strings in the form of op,string
              op can be any operation letter not currently in use
              string can be any operations. More than one <op=s> can be used
              Examples: -c x=kt,y=tk,v=dwd or -c x=kt -c y=tk -c v=dwd
  -t        tile mode. when input is a 2D tiling. unsets -g
              set if seed of Z is detected
  -r        execute operations in reverse order (left to right)
  -u        make final product be averge unit edge length
  -v        verbose output
  -i <itrs> maximum planarize iterations. -1 for unlimited (default: %d)
            WARNING: unstable models may not finish unless -i is set
  -l <lim>  minimum distance change to terminate planarization, as negative
               exponent (default: %d giving %.0e)
            WARNING: high values can cause non-terminal behaviour. Use -i
  -z <nums> number of iterations between status reports (implies termination
            check) (0 for final report only, -1 for no report), optionally
            followed by a comma and the number of iterations between
            termination checks (0 for report checks only) (default: %d,%d)
  -o <file> write output to file (default: write to standard output)

Coloring Options (run 'off_util -H color' for help on color formats)
  -V <col>  vertex color (default: gold)
  -E <col>  edge color   (default: lightgray)
  -f <mthd> mthd is face coloring method using color in map (default: n)
               key word: none - sets no color
               n - by number of sides
               s - symmetric coloring
               u - unique coloring
               o - newly created faces by operation
               w - resolve color indexes (overrides -V)
%s
               (when -f w is set)
  -R <opt>  built in seed coloring: one=1, unique=2, symmetry=3 (default: 1)
  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)
  -O <strg> face transparency pattern string (-f n only). valid values
               0 - map color alpha value, 1 -T alpha applied (default: '1')
  -m <maps> color maps for faces to be tried in turn (default: m1, for -g, m2)
               keyword m1: red,darkorange1,yellow,darkgreen,cyan,blue,magenta,
                           white,gray50,black
               keyword m2: red,blue,green,yellow,brown,magenta,purple,grue,
                           gray,orange (from George Hart's original applet)

)",
          prog_name(), help_ver_text, it_ctrl.get_max_iters(),
          it_ctrl.get_sig_digits(), it_ctrl.get_test_val(),
          it_ctrl.get_status_check_and_report_iters(),
          it_ctrl.get_status_check_only_iters(),
          TilingColoring::get_option_help('C').c_str());
}

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;
  int op_term = 0;

  string alphas_in_use;

  string map_file;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hHsgtruvc:l:i:z:f:C:R:V:E:T:O:m:o:")) !=
         -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'H':
      extended_help();
      exit(0);

    case 's':
      resolve_ops = true;
      break;

    case 'g':
      hart_mode = true;
      resolve_ops = true;
      break;

    case 't':
      tile_mode = true;
      break;

    case 'r':
      reverse_ops = true;
      break;

    case 'u':
      unitize = true;
      break;

    case 'v':
      verbosity = true;
      break;

    case 'c': {
      char operation = '\0';
      string user_op;

      Split parts(optarg, ",");
      unsigned int parts_sz = parts.size();

      for (unsigned int i = 0; i < parts_sz; i++) {
        op_term++;

        Split parts2(parts[i], "=");
        unsigned int parts2_sz = parts2.size();

        for (unsigned int j = 0; j < parts2_sz; j++) {
          if (j == 0) {
            if (strlen(parts2[j]) != 1)
              error(msg_str("term %d: operation must be one character '%s'",
                            op_term, parts2[j]),
                    c);
            else if (!isalpha(parts2[j][0])) {
              error(msg_str("term %d: operation must be alphabetic '%s'",
                            op_term, parts2[j]),
                    c);
            }
            else if (find_operator(parts2[j][0]) || find_seed(parts2[j][0]) ||
                     (alphas_in_use.find(parts2[j][0]) != string::npos)) {
              error(msg_str("term %d: operation character already in use '%s'",
                            op_term, parts2[j]),
                    c);
            }
            else {
              operation = parts2[j][0];
              alphas_in_use += operation;
              alpha_user += operation;
            }
          }
          else if (j == 1) {
            user_op = parts2[j];
          }
          else {
            error(msg_str("term %d: unexpected parameter '%s'", op_term,
                          parts2[j]),
                  c);
          }
        }

        char seed_test = '\0';
        int dummy1;
        string dummy2;
        string error_msg =
            validate_cn_string(user_op, operations_user[operation], seed_test,
                               dummy1, dummy2)
                .msg();
        if (!error_msg.empty())
          error(msg_str("term %d: %s", op_term, error_msg.c_str()));
        if (seed_test)
          error(
              msg_str("term %d: cannot contain a seed: %c", op_term, seed_test),
              c);
      }

      break;
    }

    case 'l':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_sig_digits(num), c);
      break;

    case 'i':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_max_iters(num), c);
      break;

    case 'z':
      print_status_or_exit(it_ctrl.set_status_checks(optarg), c);
      break;

    case 'f':
      if (!strcasecmp(optarg, "none"))
        face_coloring_method = '\0';
      else if (strspn(optarg, "nsuow") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("invalid face Coloring method '%s'", optarg), c);
      else {
        face_coloring_method = *optarg;
      }
      break;

    case 'C': {
      print_status_or_exit(col_type.read(optarg), c);
      break;
    }

    case 'R': {
      string arg_id;
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "one=1|unique=2|symmetry=3",
                                      argmatch_add_id_maps),
                           c);
      seed_coloring_method = atoi(arg_id.c_str());
      break;
    }

    case 'V':
      print_status_or_exit(vert_col.read(optarg), c);
      break;

    case 'E':
      print_status_or_exit(edge_col.read(optarg), c);
      break;

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
      break;

    case 'O':
      if (strspn(optarg, "01") != strlen(optarg))
        error(msg_str("transparency string is '%s' must consist of "
                      "0 and 1's",
                      optarg),
              c);
      face_pattern = optarg;
      break;

    case 'm':
      map_file = optarg;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 2)
    error("too many arguments");

  if (argc - optind > 0)
    cn_string = argv[optind];
  else
    error("no Conway Notation string given");

  print_status_or_exit(
      validate_cn_string(cn_string, operations, seed, seed_size, alpha_user));

  if (resolve_ops) {
    cn_string = resolved_cn_string(cn_string);

    // operations is a vector of objects
    for (auto &operation : operations)
      delete operation;
    operations.clear();

    // revalidate (should be valid) to rebuild operations table
    print_status_or_exit(
        validate_cn_string(cn_string, operations, seed, seed_size, alpha_user));
  }

  if (argc - optind == 2) {
    ifile = argv[++optind];
    if (seed)
      error(msg_str("seed '%c' was specified so input file '%s' is unexpected",
                    seed, ifile.c_str()));
  }

  // operations can be done in reverse order
  if (!reverse_ops) {
    sort(operations.begin(), operations.end(), cmp_ops);
    for (unsigned int i = 0; i < alpha_user.length(); i++)
      sort(operations_user[alpha_user[i]].begin(),
           operations_user[alpha_user[i]].end(), cmp_ops);
  }

  // force tile mode if using polygon
  if (seed == 'Z')
    tile_mode = true;

  if (tile_mode) {
    warning("in tile mode", 't');
    if (hart_mode) {
      warning("polygons will not process correctly with George Hart "
              "algorithms. turned off",
              'g');
      hart_mode = false;
    }
    planarize_method = 'a';
  }

  if (hart_mode) {
    if ((face_coloring_method == 'o') || (face_coloring_method == 'w'))
      error("when -g set, face coloring methods o and w are invalid", 'f');
  }

  // when use George Hart algorithms, use map he used on line
  if (!map_file.size())
    map_file = (hart_mode) ? "m2" : "m1";

  if (map_file == "m1" || map_file == "m2") {
    auto *col_map = new ColorMapMap;
    if (map_file == "m1") {
      col_map->set_col(0, Color(255, 0, 0));     // 3-sided faces red
      col_map->set_col(1, Color(255, 127, 0));   // 4-sided faces darkoranage1
      col_map->set_col(2, Color(255, 255, 0));   // 5-sided faces yellow
      col_map->set_col(3, Color(0, 100, 0));     // 6-sided faces darkgreen
      col_map->set_col(4, Color(0, 255, 255));   // 7-sided faces cyan
      col_map->set_col(5, Color(0, 0, 255));     // 8-sided faces blue
      col_map->set_col(6, Color(255, 0, 255));   // 9-sided faces magenta
      col_map->set_col(7, Color(255, 255, 255)); // 10-sided faces white
      col_map->set_col(8, Color(127, 127, 127)); // 11-sided faces gray50
      col_map->set_col(9, Color(0, 0, 0));       // 12-sided faces black
    }
    else if (map_file == "m2") {
      auto *col_map0 = new ColorMapMap;
      col_map0->set_col(0, Color(0.9, 0.3, 0.3));   // 3-sided faces red
      col_map0->set_col(1, Color(0.4, 0.4, 1.0));   // 4-sided faces blue
      col_map0->set_col(2, Color(0.2, 0.9, 0.3));   // 5-sided faces green
      col_map0->set_col(3, Color(0.9, 0.9, 0.2));   // 6-sided faces yellow
      col_map0->set_col(4, Color(0.5, 0.25, 0.25)); // 7-sided faces brown
      col_map0->set_col(5, Color(0.8, 0.2, 0.8));   // 8-sided faces magenta
      col_map0->set_col(6, Color(0.5, 0.2, 0.8));   // 9-sided faces purple
      col_map0->set_col(7, Color(0.1, 0.9, 0.9));   // 10-sided faces grue
      col_map0->set_col(8, Color(0.5, 0.5, 0.5));   // 11-sided faces gray
      col_map0->set_col(9, Color(1.0, 0.6, 0.1));   // 12-sided faces orange
      map.add_cmap(col_map0);

      // George Hart had all higher faces at gray
      col_map->set_col(0, Color(0.5, 0.5, 0.5)); // 13-sided faces and higher
    }
    map.add_cmap(col_map);

    // append a spread map instead of wrapping
    ColorMap *spread_map = colormap_from_name("spread");
    map.add_cmap(spread_map);
  }
  else
    print_status_or_exit(map.init(map_file.c_str()), 'm');

  epsilon = it_ctrl.get_test_val();
}

void verbose(const char &operation, const cn_opts &opts, const int &sub1 = 1,
             const int &sub2 = 2, const int &sides = -1)
{
  if (opts.verbosity) {
    string operator_name;
    string operator_short;

    int last_op =
        sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
    for (int i = 0; i < last_op; i++) {
      if (operation == conway_operator_list[i].operator_short[0]) {
        operator_name = conway_operator_list[i].operator_name;
        operator_short = conway_operator_list[i].operator_short;
        break;
      }
    }

    string hart_operators = "agkp";
    string hart_string;
    if (opts.hart_mode && (hart_operators.find(operation) != string::npos))
      hart_string = "(hart)";

    // special cases
    string buf;
    bool user_op = false;
    if (opts.alpha_user.find(operation) != string::npos) {
      user_op = true;
      operator_name = "user operation";
    }
    else if (operation == '_')
      operator_name = "planarizing ...";
    else if (operation == '+')
      operator_name = "orient positive mode";
    else if (operation == '-')
      operator_name = "orient negative mode";
    else if (operation == '@')
      operator_name = "non-orientable geometry";
    else if (operation == '$')
      operator_name = "done.";
    // if L0 is specially named
    else if (operation == 'L' && sub1 == 0) {
      operator_name = "joined-lace";
    }
    // all other case show numbers when allowed
    else {
      int sub1a = find_sub1_allowed(operation);
      int sub2a = find_sub2_allowed(operation);
      int sidesa = find_sides_allowed(operation);
      if (sub1a != -1 || sub2a != -1) {
        buf = "(";
        if (sub1a != -1)
          buf += std::to_string(sub1);
        if (sub2a != -1)
          buf += "," + std::to_string(sub2);
        if (sub1a != -1 || sub2a != -1)
          buf += ") ";
      }
      if (sidesa != -1 && sides != -1)
        buf += "sides = " + std::to_string(sides);
    }

    string op_string;
    if (user_op) {
      op_string = "(";
      op_string += operation;
      op_string += ")";
    }
    else if (!operator_short.empty())
      op_string = "(" + operator_short + ")";
    if (!op_string.empty())
      op_string += " ";

    fprintf(stderr, "%s%s %s %s\n", op_string.c_str(), operator_name.c_str(),
            buf.c_str(), hart_string.c_str());
  }
}

void unitize_edges(Geometry &geom)
{
  GeometryInfo info(geom);
  if (info.num_iedges() > 0) {
    double val = info.iedge_length_lims().sum / info.num_iedges();
    geom.transform(Trans3d::scale(1 / val));
  }
}

void centroid_to_origin(Geometry &geom)
{
  geom.transform(Trans3d::translate(-centroid(geom.verts())));
}

/*
// RK - average radius rather than maximum has more reliability than max
void unitize_vertex_radius(Geometry &geom)
{
  GeometryInfo info(geom);
  info.set_center(geom.centroid());
  // geom.transform(Trans3d::scale(1 / info.vert_dist_lims().max));
  double avg = info.vert_dist_lims().sum / info.num_verts();
  geom.transform(Trans3d::scale(1 / avg));
}
*/

void cn_planarize(Geometry &geom, char planarize_method, const cn_opts &opts)
{
  // if the model becomes open mid-processing, sand and fill planar works
  GeometryInfo info(geom);
  if ((planarize_method != 'a') && !info.is_closed()) {
    planarize_method = 'a';
  }

  if (opts.it_ctrl.get_max_iters() != 0) {
    verbose('_', opts);
    if (planarize_method == 'q') {
      planarize_bd(geom, opts.it_ctrl);
    }
    else if (planarize_method == 'a') {
      planarize_unit(geom, opts.it_ctrl);
    }
  }
}

void get_seed(Geometry &geom, const cn_opts &opts)
{
  string uniforms = "TCOID";

  if (uniforms.find(opts.seed) != string::npos) {
    switch (opts.seed) {
    case 'T':
      geom.read_resource("std_tet");
      break;

    case 'C':
      geom.read_resource("std_cube");
      break;

    case 'O':
      geom.read_resource("std_oct");
      break;

    case 'I':
      geom.read_resource("std_ico");
      break;

    case 'D':
      geom.read_resource("std_dod");
      break;
    }
  }
  else {
    Polygon pgon(opts.seed_size, 1);

    switch (opts.seed) {
    case 'P':
      pgon.set_type(Polygon::prism);
      break;

    case 'A':
      pgon.set_type(Polygon::antiprism);
      break;

    case 'Y':
      pgon.set_type(Polygon::pyramid);
      break;

    case 'Z':
      pgon.set_type(Polygon::dihedron);
      pgon.set_subtype(Polygon::sub_dihedron_polygon);
      break;

    // a polyhedron of random points
    case 'R':
      Random rnd;
      rnd.time_seed();
      for (int i = 0; i < opts.seed_size; i++)
        geom.add_vert(Vec3d::random(rnd).unit());
      geom.set_hull();
      break;
    }

    pgon.set_edge(0, 1.0);

    if (opts.seed == 'Y' && opts.seed_size > 5)
      // Based on circumradius
      pgon.set_height(0, (1 / sin(M_PI / opts.seed_size)) / 2);
    // inradius
    // poly->set_height((1/tan(M_PI/seed_size))/2);
    else
      pgon.set_edge(1, 1.0);

    pgon.make_poly(geom);

    /* RK - if polygon size 2 was allowed, caused too much trouble
    if ((opts.seed_size == 2) && (opts.seed == 'P' || opts.seed == 'Y'))
      geom.transform(Trans3d::rotate(deg2rad(90), 0, 0));
    */
  }

  // by default seed will be all one color
  Coloring clrng(&geom);
  if (opts.seed_coloring_method == 1) {
    Color col = opts.map.get_col(0);
    clrng.vef_one_col(col, col, col);
  }
  else if (opts.seed_coloring_method == 2) {
    // unique
    clrng.v_unique(true);
    clrng.e_unique(true);
    clrng.f_unique(true);
  }
  else if (opts.seed_coloring_method == 3) {
    // color by symmetry
    Symmetry sym;
    vector<vector<set<int>>> sym_equivs;
    sym.init(geom, &sym_equivs);
    clrng.v_sets(sym_equivs[0], true);
    clrng.e_sets(sym_equivs[1], true);
    clrng.f_sets(sym_equivs[2], true);
  }
}

// RK - for hart code
void build_new_faces(map<string, map<string, string>> &faces_table,
                     map<string, int> verts_table,
                     vector<vector<int>> &faces_new)
{
  map<string, map<string, string>>::iterator ft;
  map<string, string>::iterator ftm;
  string face_name;
  for (ft = faces_table.begin(); ft != faces_table.end(); ft++) {
    for (ftm = ft->second.begin(); ftm != ft->second.end(); ftm++) {
      if (face_name != ft->first) {
        face_name = ft->first;
        string v0 = faces_table[face_name][ftm->first];
        string v = v0;
        vector<int> face;
        do {
          face.push_back(verts_table[v]);
          v = faces_table[face_name][v];
        } while (v != v0);
        if (face.size() > 2) // make sure face is valid
          faces_new.push_back(face);
        face.clear();
      }
    }
  }
}

// hart_ code ported from George Hart java
void hart_ambo(Geometry &geom)
{
  vector<vector<int>> &faces = geom.raw_faces();
  vector<Vec3d> &verts = geom.raw_verts();

  map<string, int> verts_table;
  map<string, map<string, string>> faces_table;
  vector<Vec3d> verts_new;

  string buf1;
  string buf2;
  string buf3;
  unsigned int vert_num = 0;
  for (unsigned int i = 0; i < faces.size(); i++) {
    int v1 = faces[i].at(faces[i].size() - 2);
    int v2 = faces[i].at(faces[i].size() - 1);
    for (unsigned int j = 0; j < faces[i].size(); j++) {
      int v3 = faces[i].at(j);
      if (v1 < v2) {
        buf1 = std::to_string((v1 < v2) ? v1 : v2) + "_" +
               std::to_string((v1 > v2) ? v1 : v2);
        verts_table[buf1] = vert_num++;
        verts_new.push_back((verts[v1] + verts[v2]) * 0.5);
      }
      buf1 = "f" + std::to_string(i);
      buf2 = std::to_string((v1 < v2) ? v1 : v2) + "_" +
             std::to_string((v1 > v2) ? v1 : v2);
      buf3 = std::to_string((v2 < v3) ? v2 : v3) + "_" +
             std::to_string((v2 > v3) ? v2 : v3);
      faces_table[buf1][buf2] = buf3;
      buf1 = "v" + std::to_string(v2);
      buf2 = std::to_string((v2 < v3) ? v2 : v3) + "_" +
             std::to_string((v2 > v3) ? v2 : v3);
      buf3 = std::to_string((v1 < v2) ? v1 : v2) + "_" +
             std::to_string((v1 > v2) ? v1 : v2);
      faces_table[buf1][buf2] = buf3;
      v1 = v2;
      v2 = v3;
    }
  }

  geom.clear_all();
  verts = verts_new;
  verts_new.clear();

  build_new_faces(faces_table, verts_table, faces);
  faces_table.clear();
  verts_table.clear();
}

void hart_gyro(Geometry &geom)
{
  vector<vector<int>> &faces = geom.raw_faces();
  vector<Vec3d> &verts = geom.raw_verts();

  map<string, int> verts_table;
  map<string, map<string, string>> faces_table;
  vector<Vec3d> verts_new;

  string buf1;
  string buf2;
  string buf3;
  unsigned int vert_num = 0;
  vector<Vec3d> centers;
  geom.face_cents(centers);
  for (unsigned int i = 0; i < faces.size(); i++) {
    buf1 = "f" + std::to_string(i);
    verts_table[buf1] = vert_num++;
    verts_new.push_back(centers[i].unit());
  }
  centers.clear();

  for (unsigned int i = 0; i < verts.size(); i++) {
    buf1 = "v" + std::to_string(i);
    verts_table[buf1] = vert_num++;
    verts_new.push_back(verts[i]);
  }

  for (unsigned int i = 0; i < faces.size(); i++) {
    int v1 = faces[i].at(faces[i].size() - 2);
    int v2 = faces[i].at(faces[i].size() - 1);
    for (unsigned int j = 0; j < faces[i].size(); j++) {
      int v3 = faces[i].at(j);
      buf1 = std::to_string(v1) + "~" + std::to_string(v2);
      verts_table[buf1] = vert_num++;
      // approx. (2/3)v1 + (1/3)v2
      verts_new.push_back(verts[v1] * 0.7 + verts[v2] * 0.3);

      buf1 = std::to_string(i) + "f" + std::to_string(v1);
      buf2 = "f" + std::to_string(i);
      buf3 = std::to_string(v1) + "~" + std::to_string(v2);
      faces_table[buf1][buf2] = buf3;
      buf2 = std::to_string(v1) + "~" + std::to_string(v2);
      buf3 = std::to_string(v2) + "~" + std::to_string(v1);
      faces_table[buf1][buf2] = buf3;
      buf2 = std::to_string(v2) + "~" + std::to_string(v1);
      buf3 = "v" + std::to_string(v2);
      faces_table[buf1][buf2] = buf3;
      buf2 = "v" + std::to_string(v2);
      buf3 = std::to_string(v2) + "~" + std::to_string(v3);
      faces_table[buf1][buf2] = buf3;
      buf2 = std::to_string(v2) + "~" + std::to_string(v3);
      buf3 = "f" + std::to_string(i);
      faces_table[buf1][buf2] = buf3;

      v1 = v2;
      v2 = v3;
    }
  }

  geom.clear_all();
  verts = verts_new;
  verts_new.clear();

  build_new_faces(faces_table, verts_table, faces);
  faces_table.clear();
  verts_table.clear();
}

void hart_kisN(Geometry &geom, int n)
{
  // default num_val was changed to 1
  if (n < 3)
    n = 0;

  vector<vector<int>> &faces = geom.raw_faces();
  vector<Vec3d> &verts = geom.raw_verts();

  vector<Vec3d> centers;
  geom.face_cents(centers);

  vector<vector<int>> faces_new;
  vector<int> face;
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (n > 0 && (int)faces[i].size() != n) {
      faces_new.push_back(faces[i]);
      continue;
    }

    verts.push_back(centers[i]);
    for (unsigned int j = 0; j < faces[i].size() - 1; j++) {
      face.push_back(verts.size() - 1);
      face.push_back(faces[i][j]);
      face.push_back(faces[i][j + 1]);
      faces_new.push_back(face);
      face.clear();
    }

    face.push_back(verts.size() - 1);
    face.push_back(faces[i][faces[i].size() - 1]);
    face.push_back(faces[i][0]);
    faces_new.push_back(face);
    face.clear();
  }
  centers.clear();

  if (faces_new.size() > 0) {
    faces.clear();
    faces = faces_new;
    faces_new.clear();
  }
}

void hart_propellor(Geometry &geom)
{
  vector<vector<int>> &faces = geom.raw_faces();
  vector<Vec3d> &verts = geom.raw_verts();

  map<string, int> verts_table;
  map<string, map<string, string>> faces_table;
  vector<Vec3d> verts_new;

  string buf1;
  string buf2;
  string buf3;
  unsigned int vert_num = 0;
  for (unsigned int i = 0; i < verts.size(); i++) {
    buf1 = "v" + std::to_string(i);
    verts_table[buf1] = vert_num++;
    verts_new.push_back(verts[i].unit());
  }

  for (unsigned int i = 0; i < faces.size(); i++) {
    int v1 = faces[i].at(faces[i].size() - 2);
    int v2 = faces[i].at(faces[i].size() - 1);
    for (unsigned int j = 0; j < faces[i].size(); j++) {
      int v3 = faces[i].at(j);
      buf1 = std::to_string(v1) + "~" + std::to_string(v2);
      verts_table[buf1] = vert_num++;
      // approx. (2/3)v1 + (1/3)v2
      verts_new.push_back(verts[v1] * 0.7 + verts[v2] * 0.3);

      buf1 = "v" + std::to_string(i);
      buf2 = std::to_string(v1) + "~" + std::to_string(v2);
      buf3 = std::to_string(v2) + "~" + std::to_string(v3);
      faces_table[buf1][buf2] = buf3;
      buf1 = std::to_string(i) + "f" + std::to_string(v2);
      buf2 = std::to_string(v1) + "~" + std::to_string(v2);
      buf3 = std::to_string(v2) + "~" + std::to_string(v1);
      faces_table[buf1][buf2] = buf3;
      buf2 = std::to_string(v2) + "~" + std::to_string(v1);
      buf3 = "v" + std::to_string(v2);
      faces_table[buf1][buf2] = buf3;
      buf2 = "v" + std::to_string(v2);
      buf3 = std::to_string(v2) + "~" + std::to_string(v3);
      faces_table[buf1][buf2] = buf3;
      buf2 = std::to_string(v2) + "~" + std::to_string(v3);
      buf3 = std::to_string(v1) + "~" + std::to_string(v2);
      faces_table[buf1][buf2] = buf3;

      v1 = v2;
      v2 = v3;
    }
  }

  geom.clear_all();
  verts = verts_new;
  verts_new.clear();

  build_new_faces(faces_table, verts_table, faces);
  faces_table.clear();
  verts_table.clear();
}

/*
// chamfer for hart code
void hart_chamfer(Geometry &geom, const cn_opts &opts)
{
  vector<Vec3d> &verts = geom.raw_verts();
  unsigned int sz = verts.size();

  // make all edges a face
  // this is a join operation but retains the original vertex indexes
  make_edges_to_faces(geom);

  project_onto_sphere(geom);

  // truncate only on the new vertices
  vector<int> v_idxs;
  for (unsigned int i = sz; i < verts.size(); i++)
    v_idxs.push_back(i);

  verbose('t', opts);
  truncate_verts(geom, v_idxs, CN_ONE_HALF);
}

// whirl for hart code
void hart_whirl(Geometry &geom, bool orientation_positive, const cn_opts &opts)
{
  unsigned int num_faces = geom.raw_faces().size();

  verbose('g', opts);
  hart_gyro(geom);

  // after intra-step operation
  GeometryInfo info(geom);
  if (!info.is_orientable())
    verbose('@', opts);
  else
    // orientation is reversed if reflected 1=positive 2=negative
    geom.orient((orientation_positive) ? 1 : 2);

  cn_planarize(geom, opts);

  // only truncate on original face centers
  vector<int> v_idxs;
  for (unsigned int i = 0; i < num_faces; i++)
    v_idxs.push_back(i);

  verbose('t', opts);
  truncate_verts(geom, v_idxs, CN_ONE_HALF, nullptr);
}
*/

// operations which can use Antiprism built in features
/*
void antiprism_dual(Geometry &geom)
{
  Geometry dual;
  centroid_to_origin(geom);
  get_dual(dual, geom, 1, Vec3d(0, 0, 0));
  geom = dual;
}
*/

void antiprism_reflect(Geometry &geom, const cn_opts &opts)
{
  if (opts.tile_mode)
    geom.transform(Trans3d::reflection(Vec3d(1, 0, 0)));
  else
    geom.transform(Trans3d::inversion());
}

// built in truncate from off_util
void antiprism_truncate(Geometry &geom, double ratio, int n)
{
  truncate_verts(geom, ratio, n);
}

void orient_planar(Geometry &geom, bool &is_orientable,
                   bool &orientation_positive, const cn_opts &opts)
{
  // local copy
  char planarize_method = opts.planarize_method;

  GeometryInfo info(geom);
  is_orientable = info.is_orientable();
  if (!is_orientable) {
    verbose('@', opts);
    // default planarization method for non-orientable geometry set to unit edge
    if (opts.planarize_method != 'a') {
      planarize_method = 'a';
    }
  }
  else
    // orientation is reversed if reflected 1=positive 2=negative
    geom.orient((orientation_positive) ? 1 : 2);

  // planarize after each step
  cn_planarize(geom, planarize_method, opts);
}

// is_orientable and orientation_positive can change
void wythoff(Geometry &geom, char operation, int sub1, int sub2, int sides,
             int &operation_number, bool &is_orientable,
             bool &orientation_positive, const cn_opts &opts)
{
  operation_number++;

  string non_color_ops = "r+-";

  // if coloring new faces, track color of current faces
  // skip for reflections, orientation (when colors are not altered)
  vector<pair<Vec3d, Color>> color_centroids;
  if ((opts.face_coloring_method == 'o') &&
      (non_color_ops.find(operation) == string::npos)) {
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      pair<Vec3d, Color> color_cent;
      color_cent.first = geom.face_cent(i);
      color_cent.second = geom.colors(FACES).get(i);
      color_centroids.push_back(color_cent);
    }
  }

  // truncate with N>1 uses Hart algorithm
  if (operation == 't' && sides > 1)
    antiprism_truncate(geom, CN_ONE_THIRD, sides);
  else if (operation == 'r') {
    antiprism_reflect(geom, opts);
    // decrimenting operation number gives consistent colors
    operation_number--;
  }
  else if (operation == '+') {
    orientation_positive = true;
    operation_number--;
  }
  else if (operation == '-') {
    orientation_positive = false;
    operation_number--;
  }
  else {
    Geometry geom_save;
    vector<int> dels;
    if (sides > 1) {
      geom_save = geom;
      // remove all faces of size sides from geom_save
      for (unsigned int i = 0; i < geom.faces().size(); i++) {
        if ((int)geom_save.faces(i).size() == sides)
          dels.push_back(i);
      }
      // if matching faces found
      if (dels.size()) {
        geom_save.del(FACES, dels);
        // place only those faces in geom
        geom = faces_to_geom(geom, dels);
      }
      // no faces to act on, loop
      else {
        return;
      }
    }

    // wythoff call requires a string
    string wythoff_op;
    wythoff_op.push_back(operation);

    if ((sub1 != 1) || (sub2 != 1))
      wythoff_op += "_" + std::to_string(sub1);
    if (sub2 != 1)
      wythoff_op += "_" + std::to_string(sub2);

    // for tile mode, use old wythoff truncate
    if (opts.tile_mode && wythoff_op == "t")
      wythoff_op = "[VE]0v0e,0V,0E";

    // fprintf(stderr, "wythoff_op = %s\n", wythoff_op.c_str());
    opts.print_status_or_exit(wythoff_make_tiling(
        geom, geom, wythoff_op, is_orientable, false, opts.col_type));

    // remove digons
    dels.clear();
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      if (geom.faces(i).size() < 3) {
        // if coloring model like wythoff, move digons to edges
        if (opts.face_coloring_method == 'w') {
          Color col = geom.colors(FACES).get(i);
          geom.add_edge(geom.faces(i), col);
        }
        dels.push_back(i);
      }
    }
    geom.del(FACES, dels);
    // remove any free vertices that were formed
    geom.del(VERTS, geom.get_info().get_free_verts());

    // check for 3 faces at an edge
    auto efpairs = geom.get_edge_face_pairs(false);
    map<vector<int>, vector<int>>::const_iterator ei;
    for (ei = efpairs.begin(); ei != efpairs.end(); ++ei) {
      if (ei->second.size() > 2) {
        opts.warning("3 or more faces to an edge");
        break;
      }
    }

    // if geom_save has geometry, part of geom was saved. remerge
    if (geom_save.verts().size())
      geom.append(geom_save);
    // set to merge in any case if duplicates are encountered
    merge_coincident_elements(geom, "vef", opts.epsilon);
  }

  // if coloring new faces, restore color of previous faces
  if ((opts.face_coloring_method == 'o') &&
      (non_color_ops.find(operation) == string::npos)) {
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      Vec3d face_centroid = geom.face_cent(i);
      bool found = false;
      for (unsigned int j = 0; j < color_centroids.size(); j++) {
        pair<Vec3d, Color> color_cent = color_centroids[j];
        if (!compare(face_centroid, color_cent.first, opts.epsilon)) {
          geom.colors(FACES).set(i, color_cent.second);
          found = true;
          break;
        }
      }
      if (!found)
        geom.colors(FACES).set(i, opts.map.get_col(operation_number));
    }
  }

  orient_planar(geom, is_orientable, orientation_positive, opts);
}

void do_operations(Geometry &geom, cn_opts &opts)
{
  bool is_orientable = true;
  bool orientation_positive = true;
  int operation_number = 0;

  // the program works better with oriented input, centroid at the origin
  verbose('+', opts);
  GeometryInfo info(geom);
  is_orientable = info.is_orientable();
  if (!is_orientable)
    opts.warning("input file contains a non-orientable geometry. output is "
                 "unpredictable");
  else
    geom.orient(1); // 1=positive

  centroid_to_origin(geom);

  for (auto operation : opts.operations) {
    verbose(operation->op, opts, operation->sub1, operation->sub2,
            operation->sides);

    bool hart_operation_done = false;

    // reflection is done in wythoff
    if (opts.hart_mode) {
      hart_operation_done = true;

      switch (operation->op) {
      // ambo
      case 'a':
        hart_ambo(geom);
        break;

      // gyro
      case 'g':
        hart_gyro(geom);
        break;

      // kis
      case 'k':
        hart_kisN(geom, operation->sides);
        break;

      // propellor
      case 'p':
        hart_propellor(geom);
        break;

      default:
        hart_operation_done = false;
      }
    }

    if (hart_operation_done) {
      // these steps are needed for hart_mode
      operation_number++;
      orient_planar(geom, is_orientable, orientation_positive, opts);
    }
    else {
      // wythoff mode
      if (opts.alpha_user.find(operation->op) == string::npos)
        wythoff(geom, operation->op, operation->sub1, operation->sub2,
                operation->sides, operation_number, is_orientable,
                orientation_positive, opts);
      else {
        for (auto operation_user : opts.operations_user[operation->op]) {
          verbose(operation_user->op, opts, operation_user->sub1,
                  operation->sub2, operation_user->sides);
          wythoff(geom, operation_user->op, operation_user->sub1,
                  operation_user->sub2, operation_user->sides, operation_number,
                  is_orientable, orientation_positive, opts);
        }
      }
    }
  }
}

void cn_coloring(Geometry &geom, const cn_opts &opts)
{
  // can't color an empty geom. on -f s it will cause segfault
  if (!geom.verts().size())
    return;

  if (opts.face_coloring_method == 'n') {
    bool trans_success = true;
    const vector<vector<int>> &faces = geom.faces();
    for (unsigned int i = 0; i < faces.size(); i++) {
      unsigned int fsz = faces[i].size();
      Color col = opts.map.get_col(fsz - 3);
      if (opts.face_opacity > -1) {
        int opq = (opts.face_pattern[fsz % opts.face_pattern.size()] == '1')
                      ? opts.face_opacity
                      : col[3];
        // map colors always are set but can be map index or invisible
        if (!col.set_alpha(opq))
          trans_success = false;
      }
      geom.colors(FACES).set(i, col);
    }
    if (!trans_success)
      opts.warning("some faces could not be made transparent", 'T');
  }
  else if (opts.face_coloring_method == 's') {
    Symmetry sym;
    vector<vector<set<int>>> sym_equivs;
    sym.init(geom, &sym_equivs);

    Coloring clrng;
    clrng.add_cmap(opts.map.clone());
    clrng.set_geom(&geom);
    clrng.f_sets(sym_equivs[2], true);

    if (opts.face_opacity > -1) {
      ColorValuesToRangeHsva valmap(msg_str("A%g", opts.face_opacity / 255.0));
      valmap.apply(geom, FACES);

      for (const auto &kp : geom.colors(FACES).get_properties()) {
        if (kp.second.is_index()) {
          opts.warning("map indexes cannot be made transparent", 'T');
          break;
        }
      }
    }
  }
  else if (opts.face_coloring_method == 'u') {
    Coloring clrng(&geom);
    clrng.add_cmap(opts.map.clone());
    clrng.f_unique(true);
    clrng.f_apply_cmap();
  }
  // set color values from map indexes from wythoff call
  else if (opts.face_coloring_method == 'w') {
    Coloring clrng(&geom);
    clrng.add_cmap(opts.map.clone());
    clrng.v_apply_cmap();
    clrng.e_apply_cmap();
    clrng.f_apply_cmap();
  }

  // check if some faces are not set for transparency warning
  if (opts.face_opacity > -1) {
    if (geom.colors(FACES).get_properties().size() < geom.faces().size())
      opts.warning("unset faces cannot be made transparent", 'T');
  }

  // color vertices and edges
  // wythoff coloring doesn't color edges
  geom.add_missing_impl_edges();
  Coloring(&geom).e_one_col(opts.edge_col);

  // wythoff coloring does color vertices
  if (opts.face_coloring_method != 'w')
    Coloring(&geom).v_one_col(opts.vert_col);
}

int main(int argc, char *argv[])
{
  cn_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  if (opts.seed)
    get_seed(geom, opts);
  else
    opts.read_or_error(geom, opts.ifile);

  // if input model is not closed, Base/Dual Planarization will not work. Switch
  // to unit edge
  GeometryInfo info(geom);
  if ((opts.planarize_method == 'q') && !info.is_closed()) {
    opts.planarize_method = 'a';
  }

  do_operations(geom, opts);

  if (opts.unitize)
    unitize_edges(geom);

  cn_coloring(geom, opts);

  opts.write_or_error(geom, opts.ofile);

  verbose('$', opts);

  return 0;
}
