/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ctcdf.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 24-Apr-2018 18:23:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "ctcdf.h"
#include "gammaln.h"
#include "ctcdf_emxutil.h"

namespace matlabautogen
{

/* Function Declarations */
static double b_eml_gammainc(double x, double a, double la, double lgap1);
static double betainc_cf(double x, double a, double b);
static double eml_gammainc(double x, double la, double lgap1);
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double x
 *                double a
 *                double la
 *                double lgap1
 * Return Type  : double
 */
static double b_eml_gammainc(double x, double a, double la, double lgap1)
{
  double rval;
  double asq;
  double stirlerr;
  static const double dv0[31] = { 0.0, 0.15342640972002736, 0.081061466795327261,
    0.054814121051917651, 0.0413406959554093, 0.033162873519936291,
    0.027677925684998338, 0.023746163656297496, 0.020790672103765093,
    0.018488450532673187, 0.016644691189821193, 0.015134973221917378,
    0.013876128823070748, 0.012810465242920227, 0.01189670994589177,
    0.011104559758206917, 0.010411265261972096, 0.0097994161261588039,
    0.0092554621827127329, 0.0087687001341393862, 0.00833056343336287,
    0.00793411456431402, 0.0075736754879518406, 0.007244554301320383,
    0.00694284010720953, 0.0066652470327076821, 0.0064089941880042071,
    0.0061717122630394576, 0.0059513701127588475, 0.0057462165130101155,
    0.0055547335519628011 };

  double t;
  double xD0;
  double vsq;
  double logpax;
  double dj;
  double n;
  int exitg1;
  boolean_T exitg2;
  double a1;
  int i2;
  int i;
  double b_logpax;
  double gold;
  if (!(a > 0.0)) {
    if (x == x) {
      rval = 1.0;
    } else {
      rval = rtNaN;
    }
  } else if (!(x > 0.0)) {
    if (x == 0.0) {
      rval = 0.0;
    } else {
      rval = rtNaN;
    }
  } else if (rtIsInf(a)) {
    if (rtIsInf(x)) {
      rval = rtNaN;
    } else {
      rval = 0.0;
    }
  } else if (rtIsInf(x)) {
    rval = 1.0;
  } else {
    if (a <= 15.0) {
      asq = 2.0 * a;
      if (asq == floor(asq)) {
        stirlerr = dv0[(int)(asq + 1.0) - 1];
      } else {
        stirlerr = ((lgap1 - (a + 0.5) * la) + a) - 0.91893853320467267;
      }
    } else {
      asq = a * a;
      stirlerr = (0.083333333333333329 + (-0.0027777777777777779 +
        (0.00079365079365079365 + (-0.00059523809523809529 +
        0.00084175084175084171 / asq) / asq) / asq) / asq) / a;
    }

    if (fabs(a - x) > 0.1 * (a + x)) {
      if (a < 2.2250738585072014E-308 * x) {
        xD0 = x;
      } else if ((x < 1.0) && (a > 1.7976931348623157E+308 * x)) {
        xD0 = (a * la - a * log(x)) - a;
      } else {
        xD0 = (a * log(a / x) + x) - a;
      }
    } else {
      t = x / a;
      asq = (1.0 - t) / (1.0 + t);
      vsq = asq * asq;
      xD0 = (a - x) * asq;
      t = xD0;
      asq = 2.0 * (a * asq);
      dj = 3.0;
      do {
        exitg1 = 0;
        asq *= vsq;
        xD0 += asq / dj;
        if (xD0 == t) {
          exitg1 = 1;
        } else {
          t = xD0;
          dj += 2.0;
        }
      } while (exitg1 == 0);
    }

    logpax = (-0.5 * (1.8378770664093453 + la) - stirlerr) - xD0;
    if (x > 1.0E+6) {
      stirlerr = sqrt(x);
      t = fabs((a - x) - 1.0) / stirlerr;
      xD0 = t * t;
      if (t < 15.0) {
        asq = 0.70710678118654746 * t;
        dj = 3.97886080735226 / (asq + 3.97886080735226);
        a1 = 0.5 * ((((((((((((((((((((((0.0012710976495261409 * (dj - 0.5) +
          0.00011931402283834095) * (dj - 0.5) + -0.0039638509736051354) * (dj -
          0.5) + -0.00087077963531729586) * (dj - 0.5) + 0.0077367252831352668) *
          (dj - 0.5) + 0.0038333512626488732) * (dj - 0.5) +
          -0.012722381378212275) * (dj - 0.5) + -0.013382364453346007) * (dj -
          0.5) + 0.016131532973325226) * (dj - 0.5) + 0.039097684558848406) *
          (dj - 0.5) + 0.0024936720005350331) * (dj - 0.5) + -0.0838864557023002)
                              * (dj - 0.5) + -0.11946395996432542) * (dj - 0.5)
                             + 0.016620792496936737) * (dj - 0.5) +
                            0.35752427444953105) * (dj - 0.5) +
                           0.80527640875291062) * (dj - 0.5) +
                          1.1890298290927332) * (dj - 0.5) + 1.3704021768233816)
                        * (dj - 0.5) + 1.313146538310231) * (dj - 0.5) +
                       1.0792551515585667) * (dj - 0.5) + 0.77436819911953858) *
                     (dj - 0.5) + 0.49016508058531844) * (dj - 0.5) +
                    0.27537474159737679) * dj * exp(-asq * asq) *
          2.5066282746310002 * exp(0.5 * xD0);
        dj = (a1 * ((xD0 - 3.0) * t) - (xD0 - 4.0)) / 6.0;
        vsq = (a1 * ((xD0 * xD0 - 9.0) * xD0 + 6.0) - ((xD0 - 1.0) * xD0 - 6.0) *
               t) / 72.0;
        asq = (a1 * (((((5.0 * xD0 + 45.0) * xD0 - 81.0) * xD0 - 315.0) * xD0 +
                      270.0) * t) - ((((5.0 * xD0 + 40.0) * xD0 - 111.0) * xD0 -
                 174.0) * xD0 + 192.0)) / 6480.0;
      } else {
        a1 = (1.0 + (-1.0 + (3.0 - 15.0 / xD0) / xD0) / xD0) / t;
        dj = (1.0 + (-4.0 + (25.0 - 210.0 / xD0) / xD0) / xD0) / xD0;
        vsq = (1.0 + (-11.0 + (130.0 - 1750.0 / xD0) / xD0) / xD0) / (xD0 * t);
        asq = (1.0 + (-26.0 + (546.0 - 11368.0 / xD0) / xD0) / xD0) / (xD0 * xD0);
      }

      if (x < a - 1.0) {
        asq = a * (((a1 / stirlerr - dj / x) + vsq / (x * stirlerr)) - asq / (x *
                    x));
        if (logpax < 709.782712893384) {
          rval = exp(logpax) * asq;
        } else {
          rval = exp(logpax + log(asq));
        }
      } else {
        asq = a * (((a1 / stirlerr + dj / x) + vsq / (x * stirlerr)) + asq / (x *
                    x));
        if (logpax < 709.782712893384) {
          b_logpax = exp(logpax) * asq;
        } else {
          b_logpax = exp(logpax + log(asq));
        }

        rval = 1.0 - b_logpax;
      }
    } else if ((x < a) || (x < 1.0)) {
      n = 1.0;
      if ((!(x < a)) && (a < 2.2250738585072014E-308) && (x >
           1.7976931348623157E+308 * a)) {
        rval = 1.0;
      } else {
        if (x > 2.2204460492503131E-16 * a) {
          dj = x / a;
          n = 2.0;
          do {
            exitg1 = 0;
            dj = x * dj / (a + (n - 1.0));
            if (dj < 2.2204460492503131E-16) {
              exitg1 = 1;
            } else {
              n++;
            }
          } while (exitg1 == 0);
        }

        asq = 0.0;
        i2 = (int)((1.0 + (-1.0 - (n - 1.0))) / -1.0);
        for (i = 0; i < i2; i++) {
          asq = x * (1.0 + asq) / (a + ((n - 1.0) + -(double)i));
        }

        asq++;
        if (logpax < 709.782712893384) {
          rval = exp(logpax) * asq;
        } else {
          rval = exp(logpax + log(asq));
        }

        if (rval > 1.0) {
          rval = 1.0;
        }
      }
    } else {
      dj = 1.0;
      n = 1.0;
      exitg2 = false;
      while ((!exitg2) && (n <= floor(a + x))) {
        dj = (a - n) * dj / x;
        if (fabs(dj) < 2.2204460492503131E-16) {
          exitg2 = true;
        } else {
          n++;
        }
      }

      if (n <= floor(a + x)) {
        asq = 1.0;
      } else {
        vsq = a - floor(a);
        if (vsq == 0.0) {
          asq = 1.0;
          n = floor(a);
        } else if (vsq == 0.5) {
          asq = 0.70710678118654746 * sqrt(2.0 * x);
          t = 3.97886080735226 / (asq + 3.97886080735226);
          asq = sqrt(3.1415926535897931 * x) * exp(x) * (2.0 * (0.5 *
            ((((((((((((((((((((((0.0012710976495261409 * (t - 0.5) +
            0.00011931402283834095) * (t - 0.5) + -0.0039638509736051354) * (t -
            0.5) + -0.00087077963531729586) * (t - 0.5) + 0.0077367252831352668)
            * (t - 0.5) + 0.0038333512626488732) * (t - 0.5) +
                             -0.012722381378212275) * (t - 0.5) +
                            -0.013382364453346007) * (t - 0.5) +
                           0.016131532973325226) * (t - 0.5) +
                          0.039097684558848406) * (t - 0.5) +
                         0.0024936720005350331) * (t - 0.5) +
                        -0.0838864557023002) * (t - 0.5) + -0.11946395996432542)
                      * (t - 0.5) + 0.016620792496936737) * (t - 0.5) +
                     0.35752427444953105) * (t - 0.5) + 0.80527640875291062) *
                   (t - 0.5) + 1.1890298290927332) * (t - 0.5) +
                  1.3704021768233816) * (t - 0.5) + 1.313146538310231) * (t -
            0.5) + 1.0792551515585667) * (t - 0.5) + 0.77436819911953858) * (t -
            0.5) + 0.49016508058531844) * (t - 0.5) + 0.27537474159737679) * t *
            exp(-asq * asq)));
          n = floor(a) + 1.0;
        } else {
          xD0 = 1.0;
          a1 = x;
          stirlerr = 0.0;
          t = 1.0;
          dj = 1.0 / x;
          n = 1.0;
          asq = dj;
          gold = 0.0;
          while (fabs(asq - gold) > 2.2204460492503131E-16 * asq) {
            gold = asq;
            asq = n - vsq;
            xD0 = (a1 + xD0 * asq) * dj;
            stirlerr = (t + stirlerr * asq) * dj;
            asq = n * dj;
            a1 = x * xD0 + asq * a1;
            t = x * stirlerr + asq * t;
            dj = 1.0 / a1;
            asq = t * dj;
            n++;
          }

          asq *= x;
          n = floor(a) + 1.0;
        }
      }

      i2 = (int)((1.0 + (-1.0 - (n - 1.0))) / -1.0);
      for (i = 0; i < i2; i++) {
        asq = 1.0 + (a - ((n - 1.0) + -(double)i)) * asq / x;
      }

      asq = asq * a / x;
      if (logpax < 709.782712893384) {
        rval = exp(logpax) * asq;
      } else {
        rval = exp(logpax + log(asq));
      }

      if (rval > 1.0) {
        rval = 1.0;
      }

      rval = 1.0 - rval;
    }
  }

  return rval;
}

/*
 * Arguments    : double x
 *                double a
 *                double b
 * Return Type  : double
 */
static double betainc_cf(double x, double a, double b)
{
  double y;
  double aplusb;
  double C;
  double Dinv;
  int m;
  int exitg1;
  double yold;
  int twom;
  double d;
  double b_y;
  aplusb = a + b;
  C = 1.0;
  Dinv = 1.0 - aplusb * x / (a + 1.0);
  y = 1.0 / Dinv;
  m = 0;
  do {
    exitg1 = 0;
    if (m < 1000) {
      yold = y;
      twom = (1 + m) << 1;
      d = (1.0 + (double)m) * (b - (1.0 + (double)m)) * x / (((a - 1.0) +
        (double)twom) * (a + (double)twom));
      b_y = d / C;
      C = 1.0 + d / C;
      Dinv = 1.0 + d / Dinv;
      y *= (1.0 + b_y) / Dinv;
      d = -(a + (1.0 + (double)m)) * (aplusb + (1.0 + (double)m)) * x / ((a +
        (double)twom) * ((a + 1.0) + (double)twom));
      C = 1.0 + d / C;
      Dinv = 1.0 + d / Dinv;
      y *= C / Dinv;
      if (fabs(y - yold) < 2.2204460492503131E-15) {
        exitg1 = 1;
      } else {
        m++;
      }
    } else {
      y = rtNaN;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return y;
}

/*
 * Arguments    : double x
 *                double la
 *                double lgap1
 * Return Type  : double
 */
static double eml_gammainc(double x, double la, double lgap1)
{
  double rval;
  double t;
  double v;
  double xD0;
  double vsq;
  double logpax;
  double term;
  double dj;
  double sqrtx;
  int exitg1;
  boolean_T exitg2;
  int i1;
  int i;
  if (!(x > 0.0)) {
    if (x == 0.0) {
      rval = 1.0;
    } else {
      rval = rtNaN;
    }
  } else if (rtIsInf(x)) {
    rval = 0.0;
  } else {
    if (fabs(0.5 - x) > 0.1 * (0.5 + x)) {
      if (0.5 < 2.2250738585072014E-308 * x) {
        xD0 = x;
      } else if ((x < 1.0) && (0.5 > 1.7976931348623157E+308 * x)) {
        xD0 = (0.5 * la - 0.5 * log(x)) - 0.5;
      } else {
        xD0 = (0.5 * log(0.5 / x) + x) - 0.5;
      }
    } else {
      t = x / 0.5;
      v = (1.0 - t) / (1.0 + t);
      vsq = v * v;
      xD0 = (0.5 - x) * v;
      t = xD0;
      term = 2.0 * (0.5 * v);
      dj = 3.0;
      do {
        exitg1 = 0;
        term *= vsq;
        xD0 += term / dj;
        if (xD0 == t) {
          exitg1 = 1;
        } else {
          t = xD0;
          dj += 2.0;
        }
      } while (exitg1 == 0);
    }

    logpax = (-0.5 * (1.8378770664093453 + la) - 0.15342640972002736) - xD0;
    if (x > 1.0E+6) {
      sqrtx = sqrt(x);
      t = fabs((0.5 - x) - 1.0) / sqrtx;
      term = t * t;
      if (t < 15.0) {
        v = 0.70710678118654746 * t;
        dj = 3.97886080735226 / (v + 3.97886080735226);
        dj = 0.5 * ((((((((((((((((((((((0.0012710976495261409 * (dj - 0.5) +
          0.00011931402283834095) * (dj - 0.5) + -0.0039638509736051354) * (dj -
          0.5) + -0.00087077963531729586) * (dj - 0.5) + 0.0077367252831352668) *
          (dj - 0.5) + 0.0038333512626488732) * (dj - 0.5) +
          -0.012722381378212275) * (dj - 0.5) + -0.013382364453346007) * (dj -
          0.5) + 0.016131532973325226) * (dj - 0.5) + 0.039097684558848406) *
          (dj - 0.5) + 0.0024936720005350331) * (dj - 0.5) + -0.0838864557023002)
                              * (dj - 0.5) + -0.11946395996432542) * (dj - 0.5)
                             + 0.016620792496936737) * (dj - 0.5) +
                            0.35752427444953105) * (dj - 0.5) +
                           0.80527640875291062) * (dj - 0.5) +
                          1.1890298290927332) * (dj - 0.5) + 1.3704021768233816)
                        * (dj - 0.5) + 1.313146538310231) * (dj - 0.5) +
                       1.0792551515585667) * (dj - 0.5) + 0.77436819911953858) *
                     (dj - 0.5) + 0.49016508058531844) * (dj - 0.5) +
                    0.27537474159737679) * dj * exp(-v * v) * 2.5066282746310002
          * exp(0.5 * term);
        vsq = (dj * ((term - 3.0) * t) - (term - 4.0)) / 6.0;
        xD0 = (dj * ((term * term - 9.0) * term + 6.0) - ((term - 1.0) * term -
                6.0) * t) / 72.0;
        v = (dj * (((((5.0 * term + 45.0) * term - 81.0) * term - 315.0) * term
                    + 270.0) * t) - ((((5.0 * term + 40.0) * term - 111.0) *
               term - 174.0) * term + 192.0)) / 6480.0;
      } else {
        dj = (1.0 + (-1.0 + (3.0 - 15.0 / term) / term) / term) / t;
        vsq = (1.0 + (-4.0 + (25.0 - 210.0 / term) / term) / term) / term;
        xD0 = (1.0 + (-11.0 + (130.0 - 1750.0 / term) / term) / term) / (term *
          t);
        v = (1.0 + (-26.0 + (546.0 - 11368.0 / term) / term) / term) / (term *
          term);
      }

      v = 0.5 * (((dj / sqrtx + vsq / x) + xD0 / (x * sqrtx)) + v / (x * x));
      if (logpax < 709.782712893384) {
        rval = exp(logpax) * v;
      } else {
        rval = exp(logpax + log(v));
      }
    } else if (x < 1.0) {
      v = 0.5 * -x;
      vsq = v / 1.5;
      dj = 2.0;
      do {
        exitg1 = 0;
        v = -x * v / dj;
        term = v / (0.5 + dj);
        vsq += term;
        if (fabs(term) <= fabs(vsq) * 2.2204460492503131E-16) {
          exitg1 = 1;
        } else {
          dj++;
        }
      } while (exitg1 == 0);

      v = 0.5 * log(x) - lgap1;
      dj = exp(v);
      if (!(dj == 1.0)) {
        if (dj - 1.0 == -1.0) {
          v = -1.0;
        } else {
          v = (dj - 1.0) * v / log(dj);
        }
      }

      rval = -((vsq + v) + vsq * v);
      if (rval < 0.0) {
        rval = 0.0;
      }
    } else {
      v = 1.0;
      dj = 1.0;
      exitg2 = false;
      while ((!exitg2) && (dj <= floor(0.5 + x))) {
        v = (0.5 - dj) * v / x;
        if (fabs(v) < 2.2204460492503131E-16) {
          exitg2 = true;
        } else {
          dj++;
        }
      }

      if (dj <= floor(0.5 + x)) {
        v = 1.0;
      } else {
        v = 0.70710678118654746 * sqrt(2.0 * x);
        t = 3.97886080735226 / (v + 3.97886080735226);
        v = sqrt(3.1415926535897931 * x) * exp(x) * (2.0 * (0.5 *
          ((((((((((((((((((((((0.0012710976495261409 * (t - 0.5) +
          0.00011931402283834095) * (t - 0.5) + -0.0039638509736051354) * (t -
          0.5) + -0.00087077963531729586) * (t - 0.5) + 0.0077367252831352668) *
          (t - 0.5) + 0.0038333512626488732) * (t - 0.5) + -0.012722381378212275)
                          * (t - 0.5) + -0.013382364453346007) * (t - 0.5) +
                         0.016131532973325226) * (t - 0.5) +
                        0.039097684558848406) * (t - 0.5) +
                       0.0024936720005350331) * (t - 0.5) + -0.0838864557023002)
                     * (t - 0.5) + -0.11946395996432542) * (t - 0.5) +
                    0.016620792496936737) * (t - 0.5) + 0.35752427444953105) *
                  (t - 0.5) + 0.80527640875291062) * (t - 0.5) +
                 1.1890298290927332) * (t - 0.5) + 1.3704021768233816) * (t -
          0.5) + 1.313146538310231) * (t - 0.5) + 1.0792551515585667) * (t - 0.5)
             + 0.77436819911953858) * (t - 0.5) + 0.49016508058531844) * (t -
          0.5) + 0.27537474159737679) * t * exp(-v * v)));
        dj = 1.0;
      }

      i1 = (int)((1.0 + (-1.0 - (dj - 1.0))) / -1.0);
      for (i = 0; i < i1; i++) {
        v = 1.0 + (0.5 - ((dj - 1.0) + -(double)i)) * v / x;
      }

      v = v * 0.5 / x;
      if (logpax < 709.782712893384) {
        rval = exp(logpax) * v;
      } else {
        rval = exp(logpax + log(v));
      }

      if (rval > 1.0) {
        rval = 1.0;
      }
    }
  }

  return rval;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = 1.0;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *v
 *                emxArray_real_T *res
 * Return Type  : void
 */
void ctcdf(const emxArray_real_T *x, const emxArray_real_T *v, emxArray_real_T
           *res)
{
  int c;
  int b_c;
  emxArray_real_T *ztemp;
  int i0;
  double y;
  double xsq;
  double absx;
  double glnb1;
  double log1mx;
  boolean_T b0;
  boolean_T guard1 = false;
  int eint;
  if (x->size[0] <= v->size[0]) {
    c = x->size[0];
  } else {
    c = v->size[0];
  }

  if (x->size[1] <= v->size[1]) {
    b_c = x->size[1];
  } else {
    b_c = v->size[1];
  }

  emxInit_real_T(&ztemp, 2);
  i0 = ztemp->size[0] * ztemp->size[1];
  ztemp->size[0] = c;
  ztemp->size[1] = b_c;
  emxEnsureCapacity((emxArray__common *)ztemp, i0, sizeof(double));
  i0 = res->size[0] * res->size[1];
  res->size[0] = c;
  res->size[1] = b_c;
  emxEnsureCapacity((emxArray__common *)res, i0, sizeof(double));
  i0 = ztemp->size[0] * ztemp->size[1];
  c = 0;
  emxFree_real_T(&ztemp);
  while (c <= i0 - 1) {
    if ((v->data[c] > 0.0) && (!rtIsNaN(x->data[c]))) {
      if (x->data[c] == 0.0) {
        res->data[c] = 0.5;
      } else if (v->data[c] > 1.0E+7) {
        y = -x->data[c] / 1.4142135623730951;

        /* ========================== COPYRIGHT NOTICE ============================ */
        /*  The algorithms for calculating ERF(X) and ERFC(X) are derived           */
        /*  from FDLIBM, which has the following notice:                            */
        /*                                                                          */
        /*  Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.       */
        /*                                                                          */
        /*  Developed at SunSoft, a Sun Microsystems, Inc. business.                */
        /*  Permission to use, copy, modify, and distribute this                    */
        /*  software is freely granted, provided that this notice                   */
        /*  is preserved.                                                           */
        /* =============================    END    ================================ */
        absx = fabs(y);
        if (rtIsNaN(y)) {
          absx = y;
        } else if (rtIsInf(y)) {
          if (y < 0.0) {
            absx = 2.0;
          } else {
            absx = 0.0;
          }
        } else if (absx < 0.84375) {
          if (absx < 1.3877787807814457E-17) {
            absx = 1.0 - y;
          } else {
            xsq = y * y;
            absx = (0.12837916709551256 + xsq * (-0.3250421072470015 + xsq *
                     (-0.02848174957559851 + xsq * (-0.0057702702964894416 + xsq
                       * -2.3763016656650163E-5)))) / (1.0 + xsq *
              (0.39791722395915535 + xsq * (0.0650222499887673 + xsq *
              (0.0050813062818757656 + xsq * (0.00013249473800432164 + xsq *
              -3.9602282787753681E-6)))));
            if (y < 0.25) {
              absx = 1.0 - (y + y * absx);
            } else {
              absx = 0.5 - (y * absx + (y - 0.5));
            }
          }
        } else if (absx < 1.25) {
          glnb1 = -0.0023621185607526594 + (absx - 1.0) * (0.41485611868374833 +
            (absx - 1.0) * (-0.37220787603570132 + (absx - 1.0) *
                            (0.31834661990116175 + (absx - 1.0) *
                             (-0.11089469428239668 + (absx - 1.0) *
                              (0.035478304325618236 + (absx - 1.0) *
                               -0.0021663755948687908)))));
          xsq = 1.0 + (absx - 1.0) * (0.10642088040084423 + (absx - 1.0) *
            (0.540397917702171 + (absx - 1.0) * (0.071828654414196266 + (absx -
            1.0) * (0.12617121980876164 + (absx - 1.0) * (0.013637083912029051 +
            (absx - 1.0) * 0.011984499846799107)))));
          if (y >= 0.0) {
            absx = 0.15493708848953247 - glnb1 / xsq;
          } else {
            absx = 1.0 + (0.84506291151046753 + glnb1 / xsq);
          }
        } else if (y < -6.0) {
          absx = 2.0;
        } else if (y >= 28.0) {
          absx = 0.0;
        } else {
          xsq = 1.0 / (absx * absx);
          if (absx < 2.8571414947509766) {
            log1mx = -0.0098649440348471482 + xsq * (-0.69385857270718176 + xsq *
              (-10.558626225323291 + xsq * (-62.375332450326006 + xsq *
              (-162.39666946257347 + xsq * (-184.60509290671104 + xsq *
              (-81.2874355063066 + xsq * -9.8143293441691455))))));
            glnb1 = 1.0 + xsq * (19.651271667439257 + xsq * (137.65775414351904
              + xsq * (434.56587747522923 + xsq * (645.38727173326788 + xsq *
              (429.00814002756783 + xsq * (108.63500554177944 + xsq *
              (6.5702497703192817 + xsq * -0.0604244152148581)))))));
          } else {
            log1mx = -0.0098649429247001 + xsq * (-0.799283237680523 + xsq *
              (-17.757954917754752 + xsq * (-160.63638485582192 + xsq *
              (-637.56644336838963 + xsq * (-1025.0951316110772 + xsq *
              -483.5191916086514)))));
            glnb1 = 1.0 + xsq * (30.338060743482458 + xsq * (325.79251299657392
              + xsq * (1536.729586084437 + xsq * (3199.8582195085955 + xsq *
              (2553.0504064331644 + xsq * (474.52854120695537 + xsq *
              -22.440952446585818))))));
          }

          if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
            xsq = frexp(absx, &eint);
            b_c = eint;
          } else {
            xsq = absx;
            b_c = 0;
          }

          xsq = floor(xsq * 2.097152E+6) / 2.097152E+6 * rt_powd_snf(2.0, b_c);
          absx = exp(-xsq * xsq - 0.5625) * exp((xsq - absx) * (xsq + absx) +
            log1mx / glnb1) / absx;
          if (y < 0.0) {
            absx = 2.0 - absx;
          }
        }

        res->data[c] = 0.5 * absx;
      } else if (v->data[c] == 1.0) {
        res->data[c] = atan(1.0 / -x->data[c]) / 3.1415926535897931;
        if (x->data[c] > 0.0) {
          res->data[c]++;
        }
      } else {
        xsq = x->data[c] * x->data[c];
        if (v->data[c] < xsq) {
          y = v->data[c] / (v->data[c] + xsq);
          absx = v->data[c] / 2.0;
          xsq = absx;
          gammaln(&xsq);
          glnb1 = 0.5;
          gammaln(&glnb1);
          log1mx = absx + 0.5;
          gammaln(&log1mx);
          xsq = (xsq + glnb1) - log1mx;
          if (!(absx >= 0.0)) {
            xsq = rtNaN;
          } else {
            if ((0.0 < y) && (y < 1.0)) {
              b0 = true;
            } else {
              b0 = false;
            }

            if (!b0) {
              if (y == 0.0) {
                xsq = 0.0;
              } else if (y == 1.0) {
                xsq = 1.0;
              } else {
                xsq = rtNaN;
              }
            } else if (absx == 0.0) {
              xsq = !(y == 0.0);
            } else if (rtIsInf(absx)) {
              xsq = !(y < 1.0);
            } else {
              guard1 = false;
              if (absx + 0.5 < 1.0E+7) {
                glnb1 = log(y);
                if (1.0 - y != 1.0) {
                  log1mx = log(1.0 - y) * (-y / ((1.0 - y) - 1.0));
                } else {
                  log1mx = -y;
                }

                if (y < (absx + 1.0) / ((absx + 0.5) + 2.0)) {
                  xsq = exp(((absx * glnb1 + 0.5 * log1mx) - log(absx)) - xsq) *
                    betainc_cf(y, absx, 0.5);
                } else {
                  xsq = 1.0 - exp(((absx * glnb1 + 0.5 * log1mx) -
                                   -0.69314718055994529) - xsq) * betainc_cf(1.0
                    - y, 0.5, absx);
                }

                if (xsq == xsq) {
                } else {
                  guard1 = true;
                }
              } else {
                guard1 = true;
              }

              if (guard1) {
                xsq = rt_powd_snf(0.5 * y, 0.33333333333333331);
                glnb1 = rt_powd_snf(absx * (1.0 - y), 0.33333333333333331);
                if (((absx + 0.5) - 1.0) * (1.0 - y) > 0.8) {
                  glnb1 = 0.70710678118654746 * -(3.0 * (0.77777777777777779 *
                    xsq - (1.0 - 1.0 / (9.0 * absx)) * glnb1) / sqrt(xsq * xsq /
                    0.5 + glnb1 * glnb1 / absx));
                  xsq = 3.97886080735226 / (fabs(glnb1) + 3.97886080735226);
                  xsq = 0.5 * ((((((((((((((((((((((0.0012710976495261409 * (xsq
                    - 0.5) + 0.00011931402283834095) * (xsq - 0.5) +
                    -0.0039638509736051354) * (xsq - 0.5) +
                    -0.00087077963531729586) * (xsq - 0.5) +
                    0.0077367252831352668) * (xsq - 0.5) + 0.0038333512626488732)
                    * (xsq - 0.5) + -0.012722381378212275) * (xsq - 0.5) +
                    -0.013382364453346007) * (xsq - 0.5) + 0.016131532973325226)
                    * (xsq - 0.5) + 0.039097684558848406) * (xsq - 0.5) +
                    0.0024936720005350331) * (xsq - 0.5) + -0.0838864557023002) *
                    (xsq - 0.5) + -0.11946395996432542) * (xsq - 0.5) +
                                        0.016620792496936737) * (xsq - 0.5) +
                                       0.35752427444953105) * (xsq - 0.5) +
                                      0.80527640875291062) * (xsq - 0.5) +
                                     1.1890298290927332) * (xsq - 0.5) +
                                    1.3704021768233816) * (xsq - 0.5) +
                                   1.313146538310231) * (xsq - 0.5) +
                                  1.0792551515585667) * (xsq - 0.5) +
                                 0.77436819911953858) * (xsq - 0.5) +
                                0.49016508058531844) * (xsq - 0.5) +
                               0.27537474159737679) * xsq * exp(-glnb1 * glnb1);
                  if (glnb1 < 0.0) {
                    xsq = 1.0 - xsq;
                  }
                } else {
                  xsq = 0.5;
                  gammaln(&xsq);
                  xsq = eml_gammainc(0.5 * (((absx + 0.5) - 1.0) * (3.0 - y) -
                    -0.5) * (1.0 - y), -0.69314718055994529, xsq);
                }
              }
            }
          }

          res->data[c] = xsq / 2.0;
        } else {
          y = xsq / (v->data[c] + xsq);
          absx = v->data[c] / 2.0;
          xsq = 0.5;
          gammaln(&xsq);
          glnb1 = absx;
          gammaln(&glnb1);
          log1mx = 0.5 + absx;
          gammaln(&log1mx);
          xsq = (xsq + glnb1) - log1mx;
          if (!(absx >= 0.0)) {
            xsq = rtNaN;
          } else {
            if ((0.0 < y) && (y < 1.0)) {
              b0 = true;
            } else {
              b0 = false;
            }

            if (!b0) {
              if (y == 0.0) {
                xsq = 1.0;
              } else if (y == 1.0) {
                xsq = 0.0;
              } else {
                xsq = rtNaN;
              }
            } else if (absx == 0.0) {
              xsq = (y < 1.0);
            } else if (rtIsInf(absx)) {
              xsq = (y == 0.0);
            } else {
              guard1 = false;
              if (0.5 + absx < 1.0E+7) {
                glnb1 = log(y);
                if (1.0 - y != 1.0) {
                  log1mx = log(1.0 - y) * (-y / ((1.0 - y) - 1.0));
                } else {
                  log1mx = -y;
                }

                if (y < 1.5 / ((0.5 + absx) + 2.0)) {
                  xsq = 1.0 - exp(((0.5 * glnb1 + absx * log1mx) -
                                   -0.69314718055994529) - xsq) * betainc_cf(y,
                    0.5, absx);
                } else {
                  xsq = exp(((0.5 * glnb1 + absx * log1mx) - log(absx)) - xsq) *
                    betainc_cf(1.0 - y, absx, 0.5);
                }

                if (xsq == xsq) {
                } else {
                  guard1 = true;
                }
              } else {
                guard1 = true;
              }

              if (guard1) {
                xsq = rt_powd_snf(absx * y, 0.33333333333333331);
                glnb1 = rt_powd_snf(0.5 * (1.0 - y), 0.33333333333333331);
                if (((0.5 + absx) - 1.0) * (1.0 - y) > 0.8) {
                  glnb1 = 0.70710678118654746 * (3.0 * ((1.0 - 1.0 / (9.0 * absx))
                    * xsq - 0.77777777777777779 * glnb1) / sqrt(xsq * xsq / absx
                    + glnb1 * glnb1 / 0.5));
                  xsq = 3.97886080735226 / (fabs(glnb1) + 3.97886080735226);
                  xsq = 0.5 * ((((((((((((((((((((((0.0012710976495261409 * (xsq
                    - 0.5) + 0.00011931402283834095) * (xsq - 0.5) +
                    -0.0039638509736051354) * (xsq - 0.5) +
                    -0.00087077963531729586) * (xsq - 0.5) +
                    0.0077367252831352668) * (xsq - 0.5) + 0.0038333512626488732)
                    * (xsq - 0.5) + -0.012722381378212275) * (xsq - 0.5) +
                    -0.013382364453346007) * (xsq - 0.5) + 0.016131532973325226)
                    * (xsq - 0.5) + 0.039097684558848406) * (xsq - 0.5) +
                    0.0024936720005350331) * (xsq - 0.5) + -0.0838864557023002) *
                    (xsq - 0.5) + -0.11946395996432542) * (xsq - 0.5) +
                                        0.016620792496936737) * (xsq - 0.5) +
                                       0.35752427444953105) * (xsq - 0.5) +
                                      0.80527640875291062) * (xsq - 0.5) +
                                     1.1890298290927332) * (xsq - 0.5) +
                                    1.3704021768233816) * (xsq - 0.5) +
                                   1.313146538310231) * (xsq - 0.5) +
                                  1.0792551515585667) * (xsq - 0.5) +
                                 0.77436819911953858) * (xsq - 0.5) +
                                0.49016508058531844) * (xsq - 0.5) +
                               0.27537474159737679) * xsq * exp(-glnb1 * glnb1);
                  if (glnb1 < 0.0) {
                    xsq = 1.0 - xsq;
                  }
                } else {
                  xsq = absx;
                  gammaln(&xsq);
                  xsq = b_eml_gammainc(0.5 * (((0.5 + absx) - 1.0) * (3.0 - y) -
                    (absx - 1.0)) * (1.0 - y), absx, log(absx), xsq);
                }
              }
            }
          }

          res->data[c] = xsq / 2.0;
        }

        if (x->data[c] > 0.0) {
          res->data[c] = 1.0 - res->data[c];
        }
      }
    } else {
      res->data[c] = rtNaN;
    }

    c++;
  }
}
}
/*
 * File trailer for ctcdf.c
 *
 * [EOF]
 */
