/*
   Copyright (c) 2007-2023, Roger Kaufman
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
#include "canonical_common.h"
#include "color_common.h"

#include <cstdio>
#include <string>
#include <vector>

using std::pair;
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
  int sub1_default;
  int sub2_default;
  int sides;
  bool hart_operator;
};

// clang-format off
ConwayOperator conway_operator_list[]{
    {"a",  "ambo",            -1, -1, -1, -1, -1, true  },
    {"B",  "bowtie",          -1, -1, -1, -1, -1, false },
    {"b",  "bevel",            1, -1,  2, -1, -1, false }, // subscript >= 1
    {"c",  "chamfer",         -1, -1, -1, -1, -1, true  },
    {"d",  "dual",            -1, -1, -1, -1, -1, true  },
    {"E",  "ethyl",           -1, -1, -1, -1, -1, false },
    {"e",  "expand",           0,  0,  2,  0, -1, false }, // two subscripts n>=0 m>=0
    {"G",  "opposite-lace",   -1, -1, -1, -1, -1, false },
    {"g",  "gyro",             1, -1,  1, -1, -1, true  }, // subscript >= 1
    {"J",  "joined-medial",   -1, -1, -1, -1, -1, false }, // replaces wiki M0
    {"j",  "join",            -1, -1, -1, -1, -1, false },
    {"K",  "stake",           -1, -1, -1, -1,  3, false }, // face sides >= 3
    {"k",  "kis",             -1, -1, -1, -1,  3, true  }, // vertex sides >= 3
    {"L_0","joined-lace",     -1, -1, -1, -1, -1, false }, // placeholder
    {"L",  "lace",             0, -1,  1, -1,  3, false }, // subscript >= 0, face sides >= 3
    {"l",  "loft",             0, -1,  1, -1,  3, false }, // subscript >= 0, face sides >= 3
    {"M",  "medial",           1, -1,  2, -1, -1, false }, // subscript >= 1
    {"m",  "meta",             1, -1,  2, -1, -1, false }, // subscript >= 1
    {"n",  "needle",          -1, -1, -1, -1, -1, false },
    {"o",  "ortho",            0,  0,  2,  0, -1, false }, // two subscripts n>=0 m>=0
    {"p",  "propeller",       -1, -1, -1, -1, -1, true  },
    {"q",  "quinto",          -1, -1, -1, -1, -1, false },
    {"r",  "reflect",         -1, -1, -1, -1, -1, false }, // only in conway
    {"S",  "seed",            -1, -1, -1, -1, -1, false },
    {"s",  "snub",             1, -1,  2, -1, -1, false }, // subscript >= 1
    {"t",  "truncate",        -1, -1, -1, -1,  2, false }, // vertex sides >= 2
    {"u",  "subdivide",        0,  0,  2,  0, -1, false }, // two subscripts n>=0 m>=0
    {"W",  "waffle",          -1, -1, -1, -1, -1, false },
    {"w",  "whirl",           -1, -1, -1, -1, -1, true  },
    {"X",  "cross",            1, -1,  2, -1, -1, false }, // subscript >= 1
    {"z",  "zip",             -1, -1, -1, -1, -1, false },
    {"+",  "orient positive", -1, -1, -1, -1, -1, false },
    {"-",  "orient negative", -1, -1, -1, -1, -1, false },
};

// operators available: FfHhiNQUVvxy

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

bool find_operator(const string &op)
{
  bool found = false;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short) {
      found = true;
      break;
    }
  }
  return found;
}

int find_sub1_allowed(const string &op)
{
  int found = -1;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short) {
      if (conway_operator_list[i].sub1 > -1)
        found = conway_operator_list[i].sub1;
      break;
    }
  }
  return found;
}

int find_sub2_allowed(const string &op)
{
  int found = -1;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short) {
      if (conway_operator_list[i].sub2 > -1)
        found = conway_operator_list[i].sub2;
      break;
    }
  }
  return found;
}

int sub1_default(const string &op)
{
  int def = -1;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short) {
      def = conway_operator_list[i].sub1_default;
      break;
    }
  }
  return def;
}

int sub2_default(const string &op)
{
  int def = -1;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short) {
      def = conway_operator_list[i].sub2_default;
      break;
    }
  }
  return def;
}

int find_sides_allowed(const string &op)
{
  int found = -1;
  unsigned int last_op =
      sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (op == conway_operator_list[i].operator_short) {
      if (conway_operator_list[i].sides > -1)
        found = conway_operator_list[i].sides;
      break;
    }
  }
  return found;
}

bool find_seed(const string &seed)
{
  bool found = false;
  unsigned int last_op = sizeof(conway_seed_list) / sizeof(conway_seed_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (seed == conway_seed_list[i].seed) {
      found = true;
      break;
    }
  }
  return found;
}

int find_seed_sides(const string &seed)
{
  int found = -1;
  unsigned int last_op = sizeof(conway_seed_list) / sizeof(conway_seed_list[0]);
  for (unsigned int i = 0; i < last_op; i++) {
    if (seed == conway_seed_list[i].seed) {
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
  string op;
  int sub1;
  int sub2;
  int sides;
  ops(int n, string o, int s1, int s2, int si)
      : op_pos(n), op(o), sub1(s1), sub2(s2), sides(si)
  {
  }
};

bool cmp_ops(const ops *a, const ops *b) { return a->op_pos > b->op_pos; }

// iterate map for user operations, to find operation because of const
bool find_user_operation(const std::map<string, vector<ops *>> &operations_user,
                         const string &op)
{
  bool found = false;
  for (std::pair<string, vector<ops *>> element : operations_user) {
    string key = element.first;
    if (key == op) {
      found = true;
      break;
    }
  }
  return found;
}

void reinitialize_vars(int &sub1, int &sub2, int &sides)
{
  sub1 = -2; // -2 means lookup not done
  sub2 = -2; // -2 means lookup not done
  sides = -1;
}

// op, sub and sides will be reset
Status write_operation(vector<ops *> &operations, int &op_count, string &op,
                       int &sub1, int &sub2, int &sides)
{
  // PATCH:: check two subscript operators for _0 not allowed (sub2 not set)
  if (find_sub2_allowed(op) != -1) {
    if (sub1 == 0 && sub2 == -2)
      return Status::error(msg_str("subscript of %d alone not allowed", sub1));
    // set sub2 default here
    if (sub2 == -2)
      sub2 = sub2_default(op);
  }

  // PATCH: special condition, subscripts of _0_0 are not allowed
  if (sub1 == 0 && sub2 == 0)
    return Status::error(
        msg_str("two subscripts of %d and %d not allowed", sub1, sub2));

  operations.push_back(new ops(op_count++, op, sub1, sub2, sides));

  op = "";
  reinitialize_vars(sub1, sub2, sides);

  return Status::ok();
}

Status
validate_cn_string(const string &cn_string, vector<ops *> &operations,
                   string &seed, int &seed_size,
                   const std::map<string, vector<ops *>> &operations_user,
                   const string &second_char)
{
  Status stat;
  if (cn_string.empty())
    return stat.set_error("no Conway Notation string given");

  int sub1;
  int sub2;
  int sides;
  reinitialize_vars(sub1, sub2, sides);

  seed = "";
  seed_size = 0;

  bool end_of_command = false;

  int require_seed_sides = -1;
  int possible_sub = -1;
  int possible_sides = -1;
  int possible_repeats = -1;

  string super_symbol = "^";
  string sub_symbol = "_";
  string sides_symbol = ":";
  bool designated_sub = false;
  bool designated_sides = false;

  string digits = "0123456789";

  int num_val1 = -1;
  int num_val2 = -1;

  string pending_op = "";
  int pending_pos = 0;

  int op_count = 0;
  for (unsigned int i = 0; i < cn_string.length(); i++) {
    // check if it is a 2 character operator
    int j = 1;
    if ((i + 1 < cn_string.length()) &&
        (isalpha(cn_string[i]) && cn_string.substr(i + 1, 1) == second_char))
      j++;

    string current = cn_string.substr(i, j);

    // if operator was 2 characters, advance i
    if (j == 2)
      i++;

    bool is_digit = (digits.find(current) != string::npos);

    // seed with no N requirement has been found, no more characters expected
    if (end_of_command)
      return stat.set_error(
          msg_str("extra characters past seed \'%s\': \'%s\' at position %d",
                  seed.c_str(), cn_string.substr(i).c_str(), i + 1));

    // a number found where it shouldn't be
    if (is_digit) {
      if (i == 0)
        return stat.set_error(
            msg_str("digit \'%s\' unexpected in first position %d",
                    current.c_str(), i + 1));

      // there is no operator or seed pending for the digit
      if (((possible_sub == -1) && (possible_sides == -1) &&
           (require_seed_sides == -1) && (possible_repeats == -1)))
        return stat.set_error(msg_str("unexpected number value "
                                      "specified: \'%s\' at position %d",
                                      current.c_str(), i + 1));
    }
    // not a number where it should be
    else {
      // a number must follow sides
      if (designated_sides)
        return stat.set_error(msg_str(
            "sides \'%s\' at position %d requires a number of %d or more",
            sides_symbol.c_str(), i, possible_sides));

      // a number must follow a subscript
      if (designated_sub)
        return stat.set_error(msg_str(
            "subscript \'%s\' at position %d requires a number of %d or more",
            sub_symbol.c_str(), i, possible_sub));

      // a number must follow a seed with N requirement
      if (require_seed_sides > -1)
        return stat.set_error(msg_str("seed \'%s\' at position %d requires a "
                                      "number of sides of %d or more",
                                      seed.c_str(), i, require_seed_sides));

      // a number must follow a superscript
      if (possible_repeats == 1)
        return stat.set_error(msg_str(
            "superscript \'%s\' at position %d requires a number of %d or more",
            super_symbol.c_str(), pending_pos, possible_repeats));
    }

    // if no more digits or delimiter, handle delayed write with current values,
    // process next character
    if (((possible_sub > -1) || (possible_sides > -1)) && !is_digit &&
        !((current == sub_symbol) || (current == sides_symbol))) {
      string error_msg =
          write_operation(operations, op_count, pending_op, sub1, sub2, sides)
              .msg();
      if (!error_msg.empty())
        return stat.set_error(msg_str("operator \'%s\' at position %d, %s",
                                      pending_op.c_str(), pending_pos,
                                      error_msg.c_str()));
      possible_sub = -1;
      possible_sides = -1;
      num_val1 = -1;
      num_val2 = -1;
    }

    // set sub1 subscript defaults here, sub2 later
    if (sub1 == -2)
      sub1 = sub1_default(current);

    // it is a sides symbol
    if (current == sides_symbol) {
      // PATCH: L_0 does not have sides value
      if (pending_op == "L" && sub1 == 0)
        return stat.set_error(
            msg_str("sides symbol \'%s\' at position %d not allowed for L_0",
                    sides_symbol.c_str(), i + 1));

      if (i == 0)
        return stat.set_error(
            msg_str("sides \'%s\' unexpected in first position %d",
                    sides_symbol.c_str(), i + 1));
      if (possible_sides == -1 || sides != 1) {
        if (possible_sides == -1)
          return stat.set_error(
              msg_str("sides symbol \'%s\' at position %d unexpected",
                      sides_symbol.c_str(), i + 1));
        // sides set and encountered again
        else if (sides != -1) {
          return stat.set_error(
              msg_str("sides symbol \'%s\' at position %d for "
                      "operator \'%s\' unexpected (sides already set: %d)",
                      sides_symbol.c_str(), i + 1, pending_op.c_str(), sides));
        }
      }
      designated_sides = true;
    }
    // it is a subscript
    else if (current == sub_symbol) {
      if (i == 0)
        return stat.set_error(
            msg_str("subscript \'%s\' unexpected in first position %d",
                    sub_symbol.c_str(), i + 1));
      if (possible_sub == -1)
        return stat.set_error(
            msg_str("subscript symbol \'%s\' at position %d unexpected",
                    sub_symbol.c_str(), i + 1));
      if (sides != -1)
        return stat.set_error(
            msg_str("subscript symbol \'%s\' at position %d found after sides",
                    sub_symbol.c_str(), i + 1));
      designated_sub = true;
    }
    // it is a superscript
    else if (current == super_symbol) {
      pending_op = current;
      pending_pos = i + 1;

      if (i == 0)
        return stat.set_error(
            msg_str("superscript \'%s\' unexpected in first position %d",
                    pending_op.c_str(), i + 1));
      if (designated_sub || designated_sides)
        return stat.set_error(msg_str(
            "superscript \'%s\' unexpected in position %d, digit expected",
            pending_op.c_str(), i + 1));

      possible_repeats = 1;
    }
    // it is a digit
    else if (is_digit) {
      // find full number string and advance counter
      int digits_start = i;
      int digits_end = i;
      while ((digits_end + 1 < (int)cn_string.length()) &&
             (digits.find(cn_string[digits_end + 1]) != string::npos))
        digits_end++;
      // set iterator to one before next character in main loop
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
          return stat.set_error(msg_str("operator \'%s\' at position %d, "
                                        "second subscript of %d unexpected",
                                        pending_op.c_str(), pending_pos,
                                        num_val));
        else {
          if (num_val < min_allowed_sub2)
            return stat.set_error(msg_str(
                "operator \'%s\' at position %d, subscript must be %d or more",
                pending_op.c_str(), pending_pos, num_val));
          else
            num_val2 = num_val;
        }
      }
      else {
        return stat.set_error(msg_str("too many subscripts for operator \'%s\' "
                                      "at position %d, value %d unexpected",
                                      pending_op.c_str(), pending_pos,
                                      num_val));
      }

      // it is a pending superscript
      if (possible_repeats == 1) {
        if (num_val1 < 1)
          return stat.set_error(
              msg_str("superscript \'%s\' at position %d requires a number "
                      "of %d or more",
                      pending_op.c_str(), pending_pos, possible_repeats));

        // repeat last operator
        string op_last = operations[operations.size() - 1]->op;
        int sub1_last = operations[operations.size() - 1]->sub1;
        int sub2_last = operations[operations.size() - 1]->sub2;
        int sides_last = operations[operations.size() - 1]->sides;
        for (int i = 0; i < num_val1 - 1; i++)
          operations.push_back(
              new ops(op_count++, op_last, sub1_last, sub2_last, sides_last));

        // need to reset as though write_operation occured
        reinitialize_vars(sub1, sub2, sides);

        num_val1 = -1;
        possible_repeats = -1;
      }
      // it is a pending seed
      else if (require_seed_sides > -1) {
        seed_size = num_val1;
        if (seed_size < require_seed_sides)
          return stat.set_error(
              msg_str("seed \'%s\' at position %d requires a number of "
                      "sides of %d or more, got %d",
                      seed.c_str(), pending_pos, require_seed_sides, num_val1));
        require_seed_sides = -1;
        end_of_command = true;
      }
      // it is a designated subscript, check before passive
      else if (designated_sub) {
        sub1 = num_val1;
        if (num_val2 != -1)
          sub2 = num_val2;
        if (sub1 < possible_sub)
          return stat.set_error(msg_str(
              "operator \'%s\' at position %d, subscript must be %d or more",
              pending_op.c_str(), pending_pos, possible_sub));
        designated_sub = false;
      }
      // it is a designated sides, check before passive
      else if (designated_sides) {
        sides = num_val1;
        if (sides < possible_sides)
          return stat.set_error(msg_str(
              "operator \'%s\' at position %d requires a number of "
              "sides of %d or more, got %d",
              pending_op.c_str(), pending_pos, possible_sides, num_val1));
        designated_sides = false;
      }
      // it is a pending sides, passive
      // checking sides first before subs gives sides precedence
      else if (possible_sides > -1) {
        sides = num_val1;
        if (sides < possible_sides)
          return stat.set_error(msg_str(
              "operator \'%s\' at position %d requires a number of "
              "sides of %d or more, got %d",
              pending_op.c_str(), pending_pos, possible_sides, num_val1));
      }
      // it is a pending subscript, passive
      else if (possible_sub > -1) {
        sub1 = num_val1;
        if (num_val2 != -1)
          sub2 = num_val2;
        if (sub1 < possible_sub)
          return stat.set_error(msg_str("operator \'%s\' at position %d "
                                        "subscript must be %d or more, got %d",
                                        pending_op.c_str(), pending_pos,
                                        possible_sub, num_val1));
      }
    }
    // it is an operator
    else if (find_operator(current)) {
      pending_op = current;
      pending_pos = i + 1;

      possible_sub = find_sub1_allowed(current);
      possible_sides = find_sides_allowed(current);

      // if an operator needs no value number, write immediately
      if (possible_sub == -1 && possible_sides == -1) {
        string error_msg =
            write_operation(operations, op_count, pending_op, sub1, sub2, sides)
                .msg();
        if (!error_msg.empty())
          return stat.set_error(msg_str("operator \'%s\' at position %d, %s",
                                        pending_op.c_str(), pending_pos,
                                        error_msg.c_str()));
      }
    }
    // it is an alpha, needs no value number, write immediately
    else if (find_user_operation(operations_user, current)) {
      pending_op = current;
      pending_pos = i + 1;

      string error_msg =
          write_operation(operations, op_count, pending_op, sub1, sub2, sides)
              .msg();
      if (!error_msg.empty())
        return stat.set_error(msg_str("operator \'%s\' at position %d, %s",
                                      pending_op.c_str(), pending_pos,
                                      error_msg.c_str()));
    }
    // it is a seed
    else if (find_seed(current)) {
      seed = current;
      pending_pos = i + 1;
      require_seed_sides = find_seed_sides(seed);
      if (require_seed_sides == -1)
        end_of_command = true;
    }
    // fell through, unused character
    else
      return stat.set_error(
          msg_str("unexpected character \'%s\' at position %d", current.c_str(),
                  i + 1));
  }

  // end of string occured

  // if sides was last character
  if (designated_sides)
    return stat.set_error(
        msg_str("sides \'%s\' at position %d requires a number of %d or more",
                sides_symbol.c_str(), pending_pos + 1, possible_sides));

  // if subscript was last character
  if (designated_sub)
    return stat.set_error(msg_str(
        "subscript \'%s\' at position %d requires a number of %d or more",
        sub_symbol.c_str(), pending_pos + 1, possible_sub));

  // it is a dangling seed that needs N
  if (require_seed_sides > -1)
    return stat.set_error(msg_str(
        "seed \'%s\' at position %d requires a number of sides of %d or more",
        seed.c_str(), pending_pos, require_seed_sides));

  // if superscript was last character
  if ((pending_op == super_symbol) && (possible_repeats == 1))
    return stat.set_error(msg_str(
        "superscript \'%s\' at position %d requires a number of %d or more",
        super_symbol.c_str(), pending_pos, possible_repeats));

  // possible pending passive operation if no seed was specified
  if ((possible_sub > -1) || (possible_sides > -1)) {
    string error_msg =
        write_operation(operations, op_count, pending_op, sub1, sub2, sides)
            .msg();
    if (!error_msg.empty())
      return stat.set_error(msg_str("operator \'%s\' at position %d, %s",
                                    pending_op.c_str(), pending_pos,
                                    error_msg.c_str()));
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
  const string target;
  const string resolve;
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

// cn_string is validated input
string resolved_cn_string(const string &cn_string)
{
  string resolve_string = cn_string;

  char super_symbol = '^';
  string digits = "0123456789";

  // string multipliers (^) must be expanded before substitution
  if (resolve_string.find_first_of(super_symbol, 0) != string::npos) {
    string expanded_string;
    char current_op;
    bool superscript_found = false;
    for (unsigned int i = 0; i < resolve_string.size(); i++) {
      if (resolve_string[i] == super_symbol) {
        superscript_found = true;
        continue;
      }

      // if superscript is found, the next characters are digits
      if (superscript_found) {
        int digits_start = i;
        int digits_end = i;
        while ((digits_end + 1 < (int)resolve_string.length()) &&
               (digits.find(resolve_string[digits_end + 1]) != string::npos))
          digits_end++;
        // set iterator to one before next character in main loop
        i = digits_end;

        int num_val = std::stoi(resolve_string.substr(
            digits_start, (digits_end - digits_start + 1)));

        for (int j = 0; j < (num_val - 1); j++)
          expanded_string += current_op;

        superscript_found = false;
      }
      else {
        // real operator so append it to expanded
        current_op = resolve_string[i];
        expanded_string += current_op;
      }
    }
    // final product
    resolve_string = expanded_string;
  }

  int num_subst = sizeof(resolve_item_list) / sizeof(resolve_item_list[0]);

  string target;
  string resolve;

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
        // RK - if so e, b, o and m are temporarily replaced with @, &, $, %
        // so that loop can continue
        if (resolve_string[pos] == 'e')
          resolve_string.replace(pos, 1, "@");
        else if (resolve_string[pos] == 'b')
          resolve_string.replace(pos, 1, "&");
        else if (resolve_string[pos] == 'o')
          resolve_string.replace(pos, 1, "$");
        else if (resolve_string[pos] == 'm')
          resolve_string.replace(pos, 1, "%");
      }
    }
  }

  // if temporary characters exists
  replace(resolve_string.begin(), resolve_string.end(), '@', 'e');
  replace(resolve_string.begin(), resolve_string.end(), '&', 'b');
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

  return resolve_string;
}

class cn_opts : public ProgramOpts {
public:
  IterationControl it_ctrl;

  string ifile;
  string ofile;

  string cn_string;            // the conway notation string
  bool resolve_ops = false;    // resolve operations in hart mode
  bool hart_mode = false;      // George Hart legacy mode
  bool tile_mode = false;      // tile mode for flat surfaces
  bool reverse_ops = false;    // reverse the notation string
  string seed;                 // the initial seed for operations
  int seed_size = 0;           // for seeds that can have a size
  char planarize_method = 'c'; // internal default planarization method
  bool unitize = false;        // sets the edge lengths to average of 1
  bool verbosity = false;      // output on screen
  string second_char = "#";    // for 2 character operations

  double eps = anti::epsilon;

  vector<ops *> operations; // operations list

  // for on the fly user operators
  std::map<string, vector<ops *>> operations_user;

  // maps are managed, set no map
  OffColor off_color = OffColor("");

  int opacity[3] = {-1, -1, -1}; // transparency from 0 to 255, for v,e,f

  TilingColoring col_type; // for wythoff

  cn_opts() : ProgramOpts("conway")
  {
    it_ctrl.set_max_iters(1000);
    it_ctrl.set_status_checks("-1,1");
    it_ctrl.set_sig_digits(int(-log(anti::epsilon) / log(10) + 0.5));
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

(Antiprism Extension: note that more seeds have since been defined)

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
cuboctahedron.  (Antiprism Extension: or "bn" where n is 1 or greater)
Note: bevel is also known as "omnitruncating" the polyhedron, or omnitruncation

d = dual   The dual of a polyhedron has a vertex for each face, and a face for
each vertex, of the original polyhedron, e.g., dC=O.  Duality is an operation
of order two, meaning for any polyhedron X, ddX=X, e.g., ddC=dO=C. 

e = expand This is Mrs. Stott's expansion operation.  Each face of X is
separated from all its neighbors and reconnected with a new 4-sided face,
corresponding to an edge of X.  An n-gon is then added to connect the 4-sided
faces at each n-fold vertex.  For example, eC is the rhombicuboctahedron.  It
turns out that eX=aaX and so eX=edX (Antiprism Extension: One subscript as "en"
or "e_n" where n is 0 or greater. Two subscripts as "en_m" or "e_n_m" where
n and m are 0 or greater)
Note: expand is also known as "cantellating" the polyhedron, or cantellation

g = gyro   The dual operation to s is g. gX=dsdX=dsX, with all 5-sided faces.
The gyrocube, gC=gO="pentagonal icositetrahedron", is dual to the snub cube.
g is like k but with the new edges connecting the face centers to the 1/3
points on the edges rather than the vertices. (Antiprism Extension: or "gn"
where n is 1 or greater)

j = join   The join operator is dual to ambo, so jX=dadX=daX.  jX is like kX
without the original edges of X.  It produces a polyhedron with one 4-sided
face for each edge of X.  For example, jC=jO is the rhombic dodecahedron.

k = kis    All faces are processed or kr = only r-sided faces are processed
The kis operation divides each n-sided face into n triangles.  A new vertex is
added in the center of each face, e.g., the kiscube, kC, has 24 triangular
faces.  The k operator is dual to t, meaning kX=dtdX.

m = meta   Dual to b, mX=dbX=kjX.  mC has 48 triangular faces.  m is like k
and o combined; new edges connect new vertices at the face centers to the old
vertices and new vertices at the edge midpoints.  mX=mdX.  mC is the
"hexakis octahedron".  (Antiprism Extension: or "mn" where n is 1 or greater)

o = ortho  Dual to e, oX=deX=jjX.  oC is the trapezoidal icositetrahedron, with
24 kite-shaped faces.  oX has the effect of putting new vertices in the middle
of each face of X and connecting them, with new edges, to the edge midpoints of
X. (Antiprism Extension: One subscript as "on" or "o_n" where n is 0 or greater
Two subscripts as "on_m" or "o_n_m" where n and m are 0 or greater)

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

t = truncate  All faces are processed or tr = only vertices of order r are
processed. Truncating a polyhedron cuts off each vertex, producing a new
n-sided face for each n-fold vertex.  The faces of the original polyhedron
still appear, but have twice as many sides, e.g., the tC has six octagonal
sides corresponding to the six squares of the C, and eight triangles
corresponding to the cube's eight vertices.


Antiprism Extension: Further operations added. Also see
http://en.wikipedia.org/wiki/Conway_polyhedron_notation

c = chamfer   New hexagonal faces are added in place of edges

B = bowtie    Bowtie like triangles divide pentagonal faces

E = ethyl     like expand but triangles are divided into 3 kites

G = opposite-lace  Similar to lace, except with new quad faces split opposite
                   L_1. (has also been referred to as L_-1, not supported)

J = joined-medial  Like medial but new rhombic faces in place of original edges

K = stake     Subdivide faces with central quads, and triangles
              All faces processed or can be "Kr" where r is 3 or greater

L_0 = joined-lace  Similar to lace, except with new quad faces produced in L_1
                   are not split

L = lace      An augmentation of each face by an antiprism, adding a twist
              smaller copy of each face, and triangles between
              Subscript as "L_n" where n is 0 or greater
              All faces processed or can be "L:r" where r is 3 or greater
              Both may be specified as "L_n:r"

l = loft      An augmentation of each face by prism, adding a smaller copy of
              each face with trapezoids between the inner and outer ones
              Subscript as "l_n" where n is 0 or greater
              All faces processed or can be "l:r" where r is 3 or greater
              Both may be specified as "l_n:r"

M = medial    Similar to meta except no diagonal edges added, creating quad
              faces. Can be "Mn" where n is 1 or greater

n = needle    Dual of truncation, triangulate with 2 triangles across every
              edge. This bisect faces across all vertices and edges, while
              removing original edges

q = quinto    ortho followed by truncation of vertices centered on original
              faces. This create 2 new pentagons for every original edge
              it effectively lines the original faces with pentagons
              Lei Williams called this operation "Pental"

S = seed      Seed form

u = subdivide Ambo while retaining original vertices. Similar to Loop
              subdivision surface for triangle face
              One subscript as "un" or "u_n" where n is 1 or greater
              Two subscript as "un_m" or "u_n_m" where n and m are 1 or greater

W = waffle    Truncation on all vertices and then all faces split into sections

w = whirl     Gyro followed by truncation of vertices centered on original
              faces. This create 2 new hexagons for every original edge
              it effectively lines the original faces with hexagons

X = cross     Combination of kis and subdivide operation. Original edges are
              divided in half, with triangle and quad faces
              Can be "Xn" where n is 1 or greater

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

b  - n may be 1 or greater (default: 2)
e  - n,m n and m may be 0 or greater except for _0 and _0_0 (default: _2_0)
g  - n may be 1 or greater (default: 1)
K  - r may be 3 or greater representing face sides
k  - r may be 3 or greater representing face sides
L  - n,r n may be 0 or greater, r may be 3 or greater
l  - n,r n may be 0 or greater, r may be 3 or greater
       both n and r may be used together as L_n:r or l_n:r (L_0 may not have r)
       without delimiters Lr and lr, r is face sides and subscript default to 1
M  - n may be 1 or greater (default: 2)
m  - n may be 1 or greater (default: 2)
o  - n,m n and m may be 0 or greater except for _0 and _0_0 (default: _2_0)
s  - n may be 1 or greater (default: 2)
t  - r may be 2 or greater representing vertex connections (2 in tiles)
u  - n,m n and m may be 0 or greater except for _0 and _0_0 (default: _2_0)
X  - n may be 1 or greater (default: 2)

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
      fprintf(stdout, "%-2s -> %-8s", resolve_item_list[k].target.c_str(),
              resolve_item_list[k].resolve.c_str());
    }
    fprintf(stdout, "\n");
  }

  fprintf(stdout, R"(
Equivalent Operations (Antiprism)

b1 = z        e1 = d        M1 = o        m1 = k        o1 = S
b2 = b        e2 = e        M2 = M        m2 = m        o2 = o
s1 = d        u1 = S        X1 = k
s2 = s        u2 = u        X2 = X

Equal but opposite handed: e_n_m = e_m_n, o_n_m = o_m_n, u_n_m = u_m_n

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
  -v        verbose output
  -z <nums> number of iterations between status reports (implies termination
            check) (0 for final report only, -1 for no report), optionally
            followed by a comma and the number of iterations between
            termination checks (0 for report checks only) (default: %d,%d)
  -l <lim>  minimum distance change to terminate planarization, as negative
              exponent (default: %d giving %.0e)
            WARNING: high values can cause non-terminal behaviour. Use -i
  -o <file> write output to file (default: write to standard output)

Conway Notation Options
  -s        apply Conway Notation string substitutions
  -g        use George Hart algorithms (sets -s)
  -c <op=s> user defined operation strings in the form of op,string
              op can be any operation letter not currently in use
              two character operations can be as alpha# e.g z#
              string can be any operations. More than one <op=s> can be used
              Examples: -c x=kt,y=tk,a#=dwd or -c x=kt -c y=tk -c a#=dwd
  -t        tile mode. when input is a 2D tiling. unsets -g
              set if seed of Z is detected
  -r        execute operations in reverse order (left to right)
  -u        make final product be of averge unit edge length
  -p <opt>  planarize method (default: c, forced to x for open geometry)
              b - base/dual (reciprocate on face centroids magnitude squared)
              c - canonicalize
              x - none
  -i <itrs> maximum planarize iterations. -1 for unlimited (default: %d)
            WARNING: unstable models may not finish unless -i is set

Coloring Options (run 'off_util -H color' for help on color formats)
keyword: none - sets no color
  -F <col>  color the faces according to: (default: n)
              a color value - apply to all faces
              n - color by number of sides
              u - unique color
              v - color with average adjacent vertex color
              s - symmetric coloring [,sub_group,conj_type]
              b - color faces by convexity using all dihedral angles
              o - newly created faces by operation
              w - use wythoff colors (overrides -V and -E)
  -E <col>  color the edges according to: (default: lightgray)
              a color value - apply to all edges
              s - symmetric coloring [,sub_group,conj_type]
              b - color edges by convexity of dihedral angle
  -V <col>  color the vertices according to: (default: gold)
              a color value - apply to all vertices
              n - color by order of vertex
              u - unique color
              s - symmetric coloring [,sub_group,conj_type]
%s
  -T <t,e>  transparency. from 0 (invisible) to 255 (opaque). element is any
            or all of, v - vertices, e - edges, f - faces, a - all (default: f)
  -m <maps> a comma separated list of color maps used to transform color
            indexes (default: colorful), a part consisting of letters from
            v, e, f, selects the element types to apply the map list to
            (default 'vef'). use map name of 'index' to output index numbers
              colorful:   red,darkorange1,yellow,darkgreen,cyan,blue,magenta,
                          white,gray50,black
              ghart:      red,blue,green,yellow,brown,magenta,purple,grue,
                          gray,orange (from George Hart's original applet)
              convexity:  white,gray50,gray25 (for -F b, -E b)
              (maps are shifted '+-3' when -F n)
              (no effect when using -F w which uses internal wythoff maps)

)",
          prog_name(), help_ver_text,
          it_ctrl.get_status_check_and_report_iters(),
          it_ctrl.get_status_check_only_iters(), it_ctrl.get_sig_digits(),
          it_ctrl.get_test_val(), it_ctrl.get_max_iters(),
          TilingColoring::get_option_help('C').c_str());
}

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;
  int op_term = 0;

  Split parts;
  Color col;
  vector<string> map_files;

  off_color.set_f_col_op('n');
  off_color.set_e_col(Color(211, 211, 211)); // lightgray
  off_color.set_v_col(Color(255, 215, 0));   // gold

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hHsgtruvc:p:l:i:z:C:V:E:F:T:m:o:")) != -1) {
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
      string operation = "";
      string user_op;

      Split parts(optarg, ",");
      unsigned int parts_sz = parts.size();

      for (unsigned int i = 0; i < parts_sz; i++) {
        op_term++;

        Split parts2(parts[i], "=");
        unsigned int parts2_sz = parts2.size();

        for (unsigned int j = 0; j < parts2_sz; j++) {
          if (j == 0) {
            if (strlen(parts2[j]) > 2)
              error(msg_str(
                        "term %d: operation must be one or two characters '%s'",
                        op_term, parts2[j]),
                    c);
            else if (!isalpha(parts2[j][0])) {
              error(msg_str("term %d: first character must be alphabetic '%s'",
                            op_term, parts2[j]),
                    c);
            }
            else if ((strlen(parts2[j]) == 2) &&
                     (string(1, parts2[j][1]) != second_char)) {
              error(msg_str("term %d: second character must be '%s': '%s'",
                            op_term, second_char.c_str(), parts2[j]),
                    c);
            }
            else if (find_operator(parts2[j]) || find_seed(parts2[j]) ||
                     (find_user_operation(operations_user, parts2[j]))) {
              error(msg_str("term %d: operation string already in use '%s'",
                            op_term, parts2[j]),
                    c);
            }
            else
              operation = parts2[j];
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

        string seed_test;
        int dummy1;
        std::map<string, vector<ops *>> dummy2;
        string error_msg =
            validate_cn_string(user_op, operations_user[operation], seed_test,
                               dummy1, dummy2, second_char)
                .msg();
        if (!error_msg.empty())
          error(msg_str("term %d: %s", op_term, error_msg.c_str()), c);
        if (!seed_test.empty())
          error(msg_str("term %d: cannot contain a seed: %s", op_term,
                        seed_test.c_str()),
                c);
      }

      break;
    }

    case 'p':
      if (strlen(optarg) == 1 && strchr("bcx", int(*optarg)))
        planarize_method = *optarg;
      else
        error("planarize method type must be b, c or x", c);
      break;

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

    case 'C': {
      print_status_or_exit(col_type.read_coloring(optarg), c);
      break;
    }

    case 'V':
      if (col.read(optarg)) {
        off_color.set_v_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.v_op_check((char *)parts[0], "nus"))
        off_color.set_v_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_v_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_v_col_op()))
        off_color.set_v_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'E':
      if (col.read(optarg)) {
        off_color.set_e_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.e_op_check((char *)parts[0], "asb"))
        off_color.set_e_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_e_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_e_col_op()))
        off_color.set_e_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'F':
      if (col.read(optarg)) {
        off_color.set_f_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.f_op_check((char *)parts[0], "nuvsbow"))
        off_color.set_f_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_f_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_f_col_op()))
        off_color.set_f_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'T': {
      int parts_sz = parts.init(optarg, ",");
      if (parts_sz > 2)
        error("the argument has more than 2 parts", c);

      print_status_or_exit(read_int(parts[0], &num), c);
      if (num < 0 || num > 255)
        error("face transparency must be between 0 and 255", c);

      // if only one part, apply to faces as default
      if (parts_sz == 1) {
        opacity[FACES] = num;
      }
      else if (parts_sz > 1) {
        if (strspn(parts[1], "vefa") != strlen(parts[1]))
          error(msg_str("transparency elements are '%s' must be any or all "
                        "from  v, e, f, a",
                        parts[1]),
                c);

        string str = parts[1];
        if (str.find_first_of("va") != string::npos)
          opacity[VERTS] = num;
        if (str.find_first_of("ea") != string::npos)
          opacity[EDGES] = num;
        if (str.find_first_of("fa") != string::npos)
          opacity[FACES] = num;
      }
      break;
    }

    case 'm':
      map_files.push_back(optarg);
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

  print_status_or_exit(validate_cn_string(
      cn_string, operations, seed, seed_size, operations_user, second_char));

  if (resolve_ops) {
    string resolve_string = resolved_cn_string(cn_string);

    if (resolve_string != cn_string) {
      warning(msg_str("operator substitutions: '%s' resolved to '%s'",
                      cn_string.c_str(),
                      (resolve_string.empty() ? "NOTHING"
                                              : resolve_string.c_str())),
              's');

      cn_string = resolve_string;

      // if the string has not been resolved to empty
      if (!cn_string.empty()) {
        // operations is a vector of objects
        for (auto &operation : operations)
          delete operation;
        operations.clear();

        // revalidate (should be valid) to rebuild operations table
        print_status_or_exit(validate_cn_string(cn_string, operations, seed,
                                                seed_size, operations_user,
                                                second_char));
      }
    }
  }

  if (argc - optind == 2) {
    ifile = argv[++optind];
    if (!seed.empty())
      error(msg_str("seed '%s' was specified so input file '%s' is unexpected",
                    seed.c_str(), ifile.c_str()));
  }

  // operations can be done in reverse order
  if (!reverse_ops) {
    sort(operations.begin(), operations.end(), cmp_ops);
    // iterate map for user operations
    for (std::pair<string, vector<ops *>> element : operations_user) {
      string key = element.first;
      vector<ops *> ops = element.second;
      sort(ops.begin(), ops.end(), cmp_ops);
      operations_user[key] = ops;
    }
  }

  // force tile mode if using polygon
  if (seed == "Z") {
    if (!tile_mode) {
      warning("using Z polygons seed changes to tile mode", 't');
      tile_mode = true;
    }
  }

  if (tile_mode) {
    if (hart_mode) {
      warning("polygons will not process correctly with George Hart "
              "algorithms. turned off",
              'g');
      hart_mode = false;
    }
    if (planarize_method != 'x')
      warning("in tile mode, planarization changed to none", 'p');
    planarize_method = 'x';
  }

  if (hart_mode) {
    if ((strchr("oO", off_color.get_f_col_op())) ||
        (strchr("wW", off_color.get_f_col_op()))) {
      warning("when -g is set, face coloring methods o/O and W are invalid",
              'f');
      off_color.set_f_col_op('n');
    }
  }

  if (strchr("wW", off_color.get_f_col_op())) {
    warning("edge and vertex coloring is overridden by wythoff colors", 'f');
    off_color.set_e_col_op('\0');
    off_color.set_v_col_op('\0');
    map_files.clear();
  }

  // set all maps in list
  for (unsigned int i = 0; i < map_files.size(); i++)
    print_status_or_exit(read_colorings(off_color.clrngs, map_files[i].c_str()),
                         'm');

  // fill in missing maps
  string default_map_name = "colorful";
  for (unsigned int i = 0; i < 3; i++) {
    string map_name = default_map_name;
    // if map is already set, skip
    if (off_color.clrngs[i].get_cmaps().size())
      continue;
    if (i == EDGES) {
      char op = off_color.get_e_col_op();
      if (op && strchr("bB", op))
        map_name = "convexity";
    }
    // faces
    else if (i == FACES) {
      char op = off_color.get_f_col_op();
      if (op && strchr("bB", op))
        map_name = "convexity";
      else if (op && strchr("nN", op)) {
        if (hart_mode)
          // when use George Hart algorithms, use map he used on line, advance
          map_name = "ghart+-3";
        else
          // advance the default map
          map_name = default_map_name + "+-3";
      }
    }
    off_color.clrngs[i].add_cmap(colormap_from_name(map_name.c_str()));
  }

  eps = it_ctrl.get_test_val();
}

void verbose(const string &operation, const cn_opts &opts, const int &sub1 = -1,
             const int &sub2 = -1, const int &sides = -1)
{
  string op = operation;

  // PATCH: L_0 is a special case
  if (operation == "L" && sub1 == 0)
    op = "L_0";

  if (opts.verbosity) {
    string operator_name;
    string operator_short;

    int last_op =
        sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
    for (int i = 0; i < last_op; i++) {
      if (op == conway_operator_list[i].operator_short) {
        operator_name = conway_operator_list[i].operator_name;
        operator_short = conway_operator_list[i].operator_short;
        break;
      }
    }

    string hart_operators = "agkp";
    string hart_string;
    if (opts.hart_mode && (hart_operators.find(op) != string::npos))
      hart_string = "(hart)";

    // iterate map for user operations, use find because of const
    bool user_op = find_user_operation(opts.operations_user, op);

    // special cases
    string buf;
    if (user_op)
      operator_name = "user operation";
    else if (op == "&")
      operator_name = "end user operation";
    else if (op == "_")
      operator_name = "planarizing " + std::to_string(sub1) + " faces ...";
    else if (op == "+")
      operator_name = "orient positive mode";
    else if (op == "-")
      operator_name = "orient negative mode";
    else if (op == "@")
      operator_name = "non-orientable geometry";
    else if (op == "$")
      operator_name = "done.";
    // all other case show numbers when allowed
    else {
      int sub1a = find_sub1_allowed(op);
      int sub2a = find_sub2_allowed(op);
      int sidesa = find_sides_allowed(op);
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
      op_string += op;
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
  // if the model becomes open mid-processing, turn off planarization
  GeometryInfo info(geom);
  if ((planarize_method != 'x') && !info.is_closed()) {
    planarize_method = 'x';
  }

  Symmetry sym;
  // opts.print_status_or_exit(sym.init(geom), 'y');
  if ((opts.it_ctrl.get_max_iters() != 0) && (planarize_method != 'x')) {
    verbose("_", opts, (int)geom.faces().size());
    if (planarize_method == 'b')
      planarize_bd(geom, opts.it_ctrl);
    /* can cause coplanarity and nonconvexity
    else if (planarize_method == 'p') {
      double plane_factor = 20.0;
      make_planar(geom, opts.it_ctrl, plane_factor / 100, sym);
    }
    */
    else if (planarize_method == 'c') {
      char initial_point_type = 'c';
      double factor = 1.0;
      double factor_max = 50.0;
      make_canonical(geom, opts.it_ctrl, factor / 100, factor_max / 100,
                     initial_point_type, sym);
    }
  }
}

void get_seed(Geometry &geom, cn_opts &opts)
{
  string uniforms = "TCOID";

  if (uniforms.find(opts.seed) != string::npos) {
    if (opts.seed == "T")
      geom.read_resource("std_tet");
    else if (opts.seed == "C")
      geom.read_resource("std_cube");
    else if (opts.seed == "O")
      geom.read_resource("std_oct");
    else if (opts.seed == "I")
      geom.read_resource("std_ico");
    else if (opts.seed == "D")
      geom.read_resource("std_dod");
  }
  else {
    Polygon pgon(opts.seed_size, 1);

    if (opts.seed == "P")
      pgon.set_type(Polygon::prism);
    else if (opts.seed == "A")
      pgon.set_type(Polygon::antiprism);
    else if (opts.seed == "Y")
      pgon.set_type(Polygon::pyramid);
    else if (opts.seed == "Z") {
      pgon.set_type(Polygon::dihedron);
      pgon.set_subtype(Polygon::sub_dihedron_polygon);
    }
    // a polyhedron of random points
    else if (opts.seed == "R") {
      Random rnd;
      rnd.time_seed();
      for (int i = 0; i < opts.seed_size; i++)
        geom.add_vert(Vec3d::random(rnd).unit());
      geom.set_hull();
    }

    pgon.set_edge(0, 1.0);

    if (opts.seed == "Y" && opts.seed_size > 5)
      // Based on circumradius
      pgon.set_height(0, (1 / sin(M_PI / opts.seed_size)) / 2);
    // inradius
    // poly->set_height((1/tan(M_PI/seed_size))/2);
    else
      pgon.set_edge(1, 1.0);

    pgon.make_poly(geom);

    /* RK - if polygon size 2 was allowed, caused too much trouble
    if ((opts.seed_size == 2) && (opts.seed == "P" || opts.seed == "Y"))
      geom.transform(Trans3d::rotate(deg2rad(90), 0, 0));
    */
  }

  // by default seed face colors will first map color
  int n = (strchr("nN", opts.off_color.get_f_col_op())) ? 3 : 0;
  Coloring(&geom).f_one_col(opts.off_color.clrngs[FACES].get_col(n));
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

  verbose("t", opts);
  truncate_verts(geom, v_idxs, CN_ONE_HALF);
}

// whirl for hart code
void hart_whirl(Geometry &geom, bool orientation_positive, const cn_opts &opts)
{
  unsigned int num_faces = geom.raw_faces().size();

  verbose("g", opts);
  hart_gyro(geom);

  // after intra-step operation
  GeometryInfo info(geom);
  if (!info.is_orientable())
    verbose("@", opts);
  else
    // orientation is reversed if reflected 1=positive 2=negative
    geom.orient((orientation_positive) ? 1 : 2);

  cn_planarize(geom, opts);

  // only truncate on original face centers
  vector<int> v_idxs;
  for (unsigned int i = 0; i < num_faces; i++)
    v_idxs.push_back(i);

  verbose("t", opts);
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
    verbose("@", opts);
    // no planarization for non-orientable geometry
    if (opts.planarize_method != 'x') {
      planarize_method = 'x';
    }
  }
  else
    // orientation is reversed if reflected 1=positive 2=negative
    geom.orient((orientation_positive) ? 1 : 2);

  // planarize after each step
  cn_planarize(geom, planarize_method, opts);
}

// is_orientable and orientation_positive can change
void wythoff(Geometry &geom, string operation, int sub1, int sub2, int sides,
             int &operation_number, bool &is_orientable,
             bool &orientation_positive, cn_opts &opts)
{
  // for color by operation number
  operation_number++;

  string non_color_ops = "r+-";

  // if coloring new faces, track color of current faces
  // skip for reflections, orientation (when colors are not altered)
  vector<pair<Vec3d, Color>> color_centroids;
  if ((strchr("oO", opts.off_color.get_f_col_op())) &&
      (non_color_ops.find(operation) == string::npos)) {
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      pair<Vec3d, Color> color_cent;
      color_cent.first = geom.face_cent(i);
      color_cent.second = geom.colors(FACES).get(i);
      color_centroids.push_back(color_cent);
    }
  }

  // truncate with N>1 uses Hart algorithm
  if (operation == "t" && sides > 1)
    antiprism_truncate(geom, CN_ONE_THIRD, sides);
  else if (operation == "r") {
    antiprism_reflect(geom, opts);
    // decrimenting operation number gives consistent colors
    operation_number--;
  }
  else if (operation == "+") {
    orientation_positive = true;
    operation_number--;
  }
  else if (operation == "-") {
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
    string wythoff_op = operation;

    if (find_sub1_allowed(operation) != -1)
      wythoff_op += "_" + std::to_string(sub1);
    if (find_sub2_allowed(operation) != -1)
      wythoff_op += "_" + std::to_string(sub2);

    // for tile mode, use old wythoff truncate
    if (opts.tile_mode && wythoff_op == "t")
      wythoff_op = "[VE]0v0e,0V,0E";

    // fprintf(stderr, "wythoff_op = %s\n", wythoff_op.c_str());
    opts.print_status_or_exit(wythoff_make_tiling(
        geom, geom, wythoff_op, is_orientable, false, opts.col_type));

    // remove digons since they can cause processing errors
    dels.clear();
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      if (geom.faces(i).size() < 3) {
        // but if coloring model like wythoff, move digons to edges
        if (strchr("wW", opts.off_color.get_f_col_op())) {
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
    merge_coincident_elements(geom, "vef", opts.eps);
  }

  // if coloring new faces, restore color of previous faces
  if ((strchr("oO", opts.off_color.get_f_col_op())) &&
      (non_color_ops.find(operation) == string::npos)) {
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      Vec3d face_centroid = geom.face_cent(i);
      bool found = false;
      for (unsigned int j = 0; j < color_centroids.size(); j++) {
        pair<Vec3d, Color> color_cent = color_centroids[j];
        if (!compare(face_centroid, color_cent.first, opts.eps)) {
          geom.colors(FACES).set(i, color_cent.second);
          found = true;
          break;
        }
      }
      if (!found) {
        Color col =
            (opts.off_color.get_f_col_op() == 'o')
                ? operation_number
                : opts.off_color.clrngs[FACES].get_col(operation_number);
        geom.colors(FACES).set(i, col);
      }
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
  verbose("+", opts);
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

    // reflection is done in wythoff()
    if (opts.hart_mode) {
      hart_operation_done = true;

      // ambo
      if (operation->op == "a")
        hart_ambo(geom);
      // gyro
      else if (operation->op == "g")
        hart_gyro(geom);
      // kis
      else if (operation->op == "k")
        hart_kisN(geom, operation->sides);
      // propellor
      else if (operation->op == "p")
        hart_propellor(geom);
      else
        hart_operation_done = false;
    }

    if (hart_operation_done) {
      // these steps are needed for hart_mode
      operation_number++;
      orient_planar(geom, is_orientable, orientation_positive, opts);
    }
    else {
      // wythoff mode
      // check size, or an empty map ends up with a map entry. memory bug?
      if (find_user_operation(opts.operations_user, operation->op)) {
        for (auto operation_user : opts.operations_user[operation->op]) {
          verbose(operation_user->op, opts, operation_user->sub1,
                  operation_user->sub2, operation_user->sides);
          wythoff(geom, operation_user->op, operation_user->sub1,
                  operation_user->sub2, operation_user->sides, operation_number,
                  is_orientable, orientation_positive, opts);
        }
        verbose("&", opts); // mention end user operation
      }
      else {
        wythoff(geom, operation->op, operation->sub1, operation->sub2,
                operation->sides, operation_number, is_orientable,
                orientation_positive, opts);
      }
    }
  }
}

void cn_coloring(Geometry &geom, cn_opts &opts)
{
  // can't color an empty geom. on -F s it will cause segfault
  if (!geom.verts().size())
    return;

  // need to create edges so all calls work
  geom.add_missing_impl_edges();

  // set color values from map indexes from wythoff call
  // overrides all other coloring
  if (strchr("wW", opts.off_color.get_f_col_op())) {
    Tiling tiling;
    Status stat = tiling.set_geom(geom, false, 0.0);
    if (stat.is_error())
      opts.print_status_or_exit(stat, 'C');
    tiling.set_coloring(opts.col_type);

    Coloring clrngs[3];
    for (int i = 0; i < 3; i++) {
      auto &clrng = clrngs[i];
      clrng.set_geom(&geom);
      if (i == VERTS)
        clrng.add_cmap(tiling.get_default_point_colormap());
      else // FACES and EDGES
        clrng.add_cmap(tiling.get_default_tile_colormap());
    }
    clrngs[FACES].f_apply_cmap();
    clrngs[EDGES].e_apply_cmap();
    clrngs[VERTS].v_apply_cmap();

    // for conway
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      Color c = geom.colors(EDGES).get(i);
      // all edges need color for later logic
      if (!c.is_set())
        geom.colors(EDGES).set(i, opts.off_color.get_e_col());
    }
  }
  else {
    char op = opts.off_color.get_f_col_op();
    if (op && strchr("bB", op))
      color_faces_by_convexity(geom, opts.off_color.clrngs[FACES], (op == 'B'),
                               opts.eps);

    // color edges 3rd
    op = opts.off_color.get_e_col_op();
    if (op && strchr("bB", op))
      color_edges_by_dihedral(geom, opts.off_color.clrngs[EDGES], (op == 'B'),
                              opts.eps);
  }

  // any other color options done by class
  Status stat;
  if (!(stat = opts.off_color.off_color_main(geom)))
    opts.error(stat.msg());

  // apply all element transparencies
  apply_transparencies(geom, opts.opacity);
}

int main(int argc, char *argv[])
{
  cn_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  if (!opts.seed.empty())
    get_seed(geom, opts);
  else
    opts.read_or_error(geom, opts.ifile);

  // if input model is not closed, Base/Dual Planarization will not work.
  // Switch to unit edge
  GeometryInfo info(geom);
  if ((opts.planarize_method != 'x') && !info.is_closed()) {
    if (opts.planarize_method != 'x')
      opts.warning("open model detected. planarization changed to none", 'p');
    opts.planarize_method = 'x';
    // turn on tile mode so truncate works correctly
    opts.tile_mode = true;
  }

  do_operations(geom, opts);

  if (opts.unitize)
    unitize_edges(geom);

  cn_coloring(geom, opts);

  opts.write_or_error(geom, opts.ofile);

  verbose("$", opts);

  return 0;
}
