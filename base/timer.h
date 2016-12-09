/*
   Copyright (c) 2003-2016, Adrian Rossiter

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

/**\file timer.h
   \brief Timing utilities
*/

#ifndef TIMER_H
#define TIMER_H

namespace anti {

struct time_val {
  long tv_sec;
  long tv_usec;
};

time_val &tv_normalise(time_val &tv);
bool operator>(const time_val &t0, const time_val &t1);
time_val operator+(const time_val &t0, const time_val &t1);
time_val operator-(const time_val &t0, const time_val &t1);
time_val to_time_val(double tm);
long to_long_usecs(time_val tv);

/// A subsecond %Timer
class Timer {
private:
  time_val end;

  int get_time(struct time_val *tp); // wrapper for subsecond current time
  void u_sleep(unsigned long usecs); // wrapper for subsecond sleep

public:
  /// Constructor
  Timer() { set_timer(0.0); }

  /// Set the %Timer.
  /**\param interval length of time the %Timer should run. */
  void set_timer(time_val interval)
  {
    time_val tv;
    get_time(&tv);
    end = tv + interval;
  }

  /// Set the %Timer.
  /**\param interval length of time in microseconds that the %Timer
   * should run. */
  void set_timer(double interval) { set_timer(to_time_val(interval)); }

  /// Increment the %Timer.
  /**\param inc length of time in microseconds that the %Timer
   * should be extended. */
  void inc_timer(double inc) { end = end + to_time_val(inc); }

  /// Check whether the %timer has finished.
  /**\return \c true if the %timer has finished, otherwise \c false. */
  bool finished()
  {
    time_val tv;
    get_time(&tv);
    return tv > end;
  }

  /// Sleep until finished.
  /**Pause program execution for the amount of time remaining
   * on the %timer. */
  void sleep_until_finished()
  {
    time_val tv;
    get_time(&tv);
    if (end > tv)
      u_sleep(to_long_usecs(end - tv));
  }
};

} // namespace anti

#endif // TIMER_H
