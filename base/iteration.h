/*
   Copyright (c) 2020, Adrian Rossiter

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

/*!\file iteration.h
   \brief iteration control and support
*/

#ifndef ITERATION_H
#define ITERATION_H

#include "../base/const.h"
#include "../base/status.h"

#include <math.h>

namespace anti {

/// Iteration control for iterative algorithms
class IterationControl {
public:
  /// Set maximum number of iterations
  /**\param iters the maximum number of iterations
   * \return status, evaluates to \c true if a valid integer
   *  was read, otherwise \c false.*/
  Status set_max_iters(int iters);

  /// Set maximum number of iterations to be unlimited
  /**\c is_end_iter() and \c is_last_iter() will always return \c false */
  void set_max_iters_unlimited() { max_iters = unlimited; }

  /// Get maximum number of iterations
  /**\return maximum number of iterations, or -1 if unlimited. */
  int get_max_iters() const
  {
    return (max_iters == unlimited) ? -1 : max_iters;
  }

  /// Set number of iterations between status checks with reports
  /**\param iters the number of iterations between status checks with reports
   * \return status, evaluates to \c true if a valid integer
   *  was read, otherwise \c false.*/
  Status set_status_check_and_report_iters(int iters);

  /// Get number of iterations between status checks with reports
  /**\return number of iterations between status checks with reports */
  int get_status_check_and_report_iters() const
  {
    return status_check_and_report_iters;
  }

  /// Set number of iterations between status checks only
  /**\param iters the number of iterations between status checks
   * \return status, evaluates to \c true if a valid integer
   *  was read, otherwise \c false.*/
  Status set_status_check_only_iters(int iters);

  /// Get number of iterations between status checks only
  /**\return number of iterations between status reports */
  int get_status_check_only_iters() const { return status_check_only_iters; }

  /// Set status check periods for reporting and check-only
  /**\param periods the number of iterations for status check with report,
   *  optionally followed by a comma and the number of iterations for
   *  status only
   * \return status, evaluates to \c true if a valid string
   *  was read, otherwise \c false.*/
  Status set_status_checks(std::string iters_str);

  /// Set number of significant digits for testing completion
  /**\param digits number of significant digits
   * \return status, evaluates to \c true if a valid integer
   *  was read, otherwise \c false.*/
  Status set_sig_digits(int digits);

  /// Get number of significant digits for testing completion
  /**\return number of significant digits. */
  int get_sig_digits() const { return sig_digits; }

  /// Set output stream for reporting
  /**\param strm output stream (or nullptr for no reporting)*/
  void set_stream(FILE *strm) { stream = strm; }

  /// Get output stream for reporting
  /**\return output stream (or nullptr for no reporting).*/
  FILE *get_stream() const { return stream; }

  /// Set iteration counter to start, with initial loop to setup variables
  /**Test for the initial setup loop with \c is_iterating() */
  void start_iter_with_setup() { set_current_iter(0); }

  /// Set iteration counter to start
  void start_iter() { set_current_iter(1); }

  /// Increment iteration counter
  void next_iter() { current_iter++; }

  /// Indicate iteration loop should finish
  bool is_done() { return is_end_iter() || is_finished(); }

  /// Get current iteration counter
  unsigned int get_current_iter() const { return current_iter; }

  /// Is the current iteration larger than the maximum iteration
  bool is_end_iter() const { return current_iter > max_iters; }

  /// Is the current iteration the maximum iteration
  bool is_last_iter() const
  {
    return current_iter == max_iters && max_iters != unlimited;
  }

  /// Is finished
  /**\return \c true if finished, otherwise \c false.*/
  bool is_finished() const { return finished; }

  /// Set finished state
  /**\param finishd \c true to set finished, \c false to set not finished */
  void set_finished(bool finishd = true) { finished = finishd; }

  /// Is the current iteration just setting up variable
  /**\return \c true a variable setup iteration, \c false a normal
   *         model-modifying iteration. */
  bool is_setup_iter() const { return current_iter == 0; }

  /// Get a value to test for completion (1e^-sig_digits)
  /**\return a test value.*/
  double get_test_val() const { return pow(10, -sig_digits); }

  /// Should the status be reported on the current iteration
  /**\return \c true the status should be checked, \c false the
   *         status does not need to be checked.*/
  bool is_status_report_iter() const
  {
    return
        // true if periodic report (>0), and
        // iters is 1 or is at the report period
        (status_check_and_report_iters > 0 &&
         ((current_iter == 1 ||
           current_iter % status_check_and_report_iters == 0)))

        || // or

        // true if periodic or final report (>=0), and
        // is the last iteration or has finished
        (status_check_and_report_iters >= 0 &&
         (is_last_iter() || is_finished()));
  }

  /// Should the status be checked on the current iteration
  /**\return \c true the status should be checked, \c false the
   *         status does not need to be checked.*/
  bool is_status_check_iter() const
  {
    return is_status_report_iter()         // status check from reporting
           ||                              // or
           (status_check_only_iters > 0 && // periodic checking and
            current_iter % status_check_only_iters == 0); // check-only period
  }

  /// Print a message to the report stream
  /**\param fmt a printf-style format string
   * \param ... the values for the format
   * \return number of characters printed */
  int print(const char *fmt, ...) const;

private:
  void set_current_iter(unsigned int iter)
  {
    current_iter = iter;
    finished = false;
  }
  Status check_reporting();

  static const unsigned unlimited;
  unsigned int current_iter = 0;
  unsigned int max_iters = 10000;
  int status_check_and_report_iters = 1000;
  int status_check_only_iters = 0;
  int sig_digits = 13;
  FILE *stream = stderr;
  bool finished = false;
};

}; // namespace anti

#endif // ITERATION_H
