/**
 * @file ColorConverter.hxx
 * @author Thomas Lindemeier
 * @brief Collection of color conversion functions.
 * @date 2019-04-14
 *
 */

#ifndef COLOR_CONVERTER_HXX
#define COLOR_CONVERTER_HXX

#include <array>
#include <cmath>
#include <limits>

namespace color
{

/**
 * @brief Class that holds color conversion functions with reference to a given
 * white point (illuminant).
 *
 * Only data types of the template form vec<Type, size_t> are supported. Only
 * floating points and data types with 3 components are supported.
 *
 * Most of the conversion and constants were taken from
 * http://brucelindbloom.com/index.html (last accessed 2019 April 14th)
 * Credits to that guy!
 *
 * @tparam Scalar Data type of each component. Only floating point values
 * supported.
 *
 * @tparam N Number of components. Has to be equal to 3.
 *
 * @tparam Vec The data structure template class of the form Vec<Scalar, 3UL>
 */
template <
  class Scalar, class Vec3,
  typename std::enable_if_t<std::is_floating_point<Scalar>::value, int> = 0>
class ColorConverter
{
  static constexpr Scalar Pi =
    static_cast<Scalar>(3.1415926535897932384626433832795);

public:
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_A   = {1.09850, 1.00000,
                                                             0.35585};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_B   = {0.99072, 1.00000,
                                                             0.85223};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_C   = {0.98074, 1.00000,
                                                             1.18232};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_D50 = {0.96422, 1.00000,
                                                               0.82521};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_D55 = {0.95682, 1.00000,
                                                               0.92149};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_D65 = {0.95047, 1.00000,
                                                               1.08883};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_D75 = {0.94972, 1.00000,
                                                               1.22638};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_E   = {1.00000, 1.00000,
                                                             1.00000};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_2   = {0.99186, 1.00000,
                                                             0.67393};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_7   = {0.95041, 1.00000,
                                                             1.08747};
  static constexpr std::array<Scalar, 3U> IM_ILLUMINANT_11  = {1.00962, 1.00000,
                                                              0.64350};

private: // private constants
  // sRGB D65, http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  static constexpr double XYZ2RGB_MATRIX[3U][3U] = {
    {3.2404542, -1.5371385, -0.4985314},
    {-0.9692660, 1.8760108, 0.0415560},
    {0.0556434, -0.2040259, 1.0572252}};

  // sRGB D65, http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  static constexpr double RGB2XYZ_MATRIX[3U][3U] = {
    {0.4124564, 0.3575761, 0.1804375},
    {0.2126729, 0.7151522, 0.0721750},
    {0.0193339, 0.1191920, 0.9503041}};

public:
  ColorConverter(const std::array<Scalar, 3U>& whitepoint = IM_ILLUMINANT_D65)
    : illuminant(whitepoint)
  {
  }

  void lab2xyz(const Vec3& Lab, Vec3& XYZ) const
  {
    // chromatic adaption, reference white
    XYZ[1] = illuminant[1] * fi((1. / 116.) * (Lab[0] + 16.)); // Y
    XYZ[0] = illuminant[0] *
             fi((1. / 116.) * (Lab[0] + 16.) + (1. / 500.) * Lab[1]); // X
    XYZ[2] = illuminant[2] *
             fi((1. / 116.) * (Lab[0] + 16.) - (1. / 200.) * Lab[2]); // Z
  }

  void xyz2lab(const Vec3& XYZ, Vec3& Lab) const
  {
    Lab[0] = 116. * f(XYZ[1] / illuminant[1]) - 16.;
    Lab[1] = 500. * (f(XYZ[0] / illuminant[0]) - f(XYZ[1] / illuminant[1]));
    Lab[2] = 200. * (f(XYZ[1] / illuminant[1]) - f(XYZ[2] / illuminant[2]));
  }

  void lab2LCHab(const Vec3& Lab, Vec3& LCHab) const
  {
    LCHab[0] = Lab[0];                                       // [0,100]
    LCHab[1] = std::sqrt(Lab[1] * Lab[1] + Lab[2] * Lab[2]); // [0,100]

    LCHab[2] = std::atan2(Lab[2], Lab[1]);
    if (LCHab[2] < 0)
      {
        LCHab[2] += Pi * 2.; // [0, 2pi]
      }
  }

  void LCHab2lab(const Vec3& LCHab, Vec3& Lab) const
  {
    Lab[0]   = LCHab[0];
    Scalar h = LCHab[2];
    if (h > Pi)
      {
        h -= Pi * 2.; // [0, 2pi]
      }
    Lab[1] = LCHab[1] * std::cos(h);
    Lab[2] = LCHab[1] * std::sin(h);
  }

  // http://en.wikipedia.org/wiki/RYB_color_model
  // http://threekings.tk/mirror/ryb_TR.pdf
  void ryb2rgb(const Vec3& ryb, Vec3& rgb) const
  {
    auto cubicInt = [](Scalar t, Scalar A, Scalar B) {
      Scalar weight = t * t * (3 - 2 * t);
      return A + weight * (B - A);
    };
    Scalar x0, x1, x2, x3, y0, y1;
    // red
    x0     = cubicInt(ryb[2], 1., 0.163);
    x1     = cubicInt(ryb[2], 1., 0.);
    x2     = cubicInt(ryb[2], 1., 0.5);
    x3     = cubicInt(ryb[2], 1., 0.2);
    y0     = cubicInt(ryb[1], x0, x1);
    y1     = cubicInt(ryb[1], x2, x3);
    rgb[0] = cubicInt(ryb[0], y0, y1);
    // green
    x0     = cubicInt(ryb[2], 1., 0.373);
    x1     = cubicInt(ryb[2], 1., 0.66);
    x2     = cubicInt(ryb[2], 0., 0.);
    x3     = cubicInt(ryb[2], 0.5, 0.094);
    y0     = cubicInt(ryb[1], x0, x1);
    y1     = cubicInt(ryb[1], x2, x3);
    rgb[1] = cubicInt(ryb[0], y0, y1);
    // blue
    x0     = cubicInt(ryb[2], 1., 0.6);
    x1     = cubicInt(ryb[2], 0., 0.2);
    x2     = cubicInt(ryb[2], 0., 0.5);
    x3     = cubicInt(ryb[2], 0., 0.);
    y0     = cubicInt(ryb[1], x0, x1);
    y1     = cubicInt(ryb[1], x2, x3);
    rgb[2] = cubicInt(ryb[0], y0, y1);
  }

  // rgb to cmy
  void rgb2cmy(const Vec3& rgb, Vec3& cmy) const
  {
    cmy[0] = 1. - rgb[0];
    cmy[1] = 1. - rgb[1];
    cmy[2] = 1. - rgb[2];
  }

  // cmy to rgb
  void cmy2rgb(const Vec3& cmy, Vec3& rgb) const
  {
    rgb[0] = 1. - cmy[0];
    rgb[1] = 1. - cmy[1];
    rgb[2] = 1. - cmy[2];
  }

  // uses sRGB chromatic adapted matrix
  void rgb2xyz(const Vec3& rgb, Vec3& XYZ) const
  {
    XYZ[0] = RGB2XYZ_MATRIX[0][0] * rgb[0] + RGB2XYZ_MATRIX[0][1] * rgb[1] +
             RGB2XYZ_MATRIX[0][2] * rgb[2];
    XYZ[1] = RGB2XYZ_MATRIX[1][0] * rgb[0] + RGB2XYZ_MATRIX[1][1] * rgb[1] +
             RGB2XYZ_MATRIX[1][2] * rgb[2];
    XYZ[2] = RGB2XYZ_MATRIX[2][0] * rgb[0] + RGB2XYZ_MATRIX[2][1] * rgb[1] +
             RGB2XYZ_MATRIX[2][2] * rgb[2];
  }

  // uses sRGB chromatic adapted matrix
  void xyz2rgb(const Vec3& XYZ, Vec3& rgb) const
  {
    rgb[0] = XYZ2RGB_MATRIX[0][0] * XYZ[0] + XYZ2RGB_MATRIX[0][1] * XYZ[1] +
             XYZ2RGB_MATRIX[0][2] * XYZ[2];
    rgb[1] = XYZ2RGB_MATRIX[1][0] * XYZ[0] + XYZ2RGB_MATRIX[1][1] * XYZ[1] +
             XYZ2RGB_MATRIX[1][2] * XYZ[2];
    rgb[2] = XYZ2RGB_MATRIX[2][0] * XYZ[0] + XYZ2RGB_MATRIX[2][1] * XYZ[1] +
             XYZ2RGB_MATRIX[2][2] * XYZ[2];
  }

  // make linear rgb, no chromatic adaption
  void srgb2rgb(const Vec3& srgb, Vec3& rgb) const
  {
    for (auto i = 0U; i < 3U; i++)
      {
        if (srgb[i] <= 0.04045)
          rgb[i] = srgb[i] / 12.92;
        else
          rgb[i] = std::pow(((srgb[i] + 0.055) / 1.055), 2.4);
      }
  }

  // make sRGB, with gamma
  void rgb2srgb(const Vec3& rgb, Vec3& srgb) const
  {
    for (auto i = 0U; i < 3U; ++i)
      {
        if (rgb[i] <= 0.0031308)
          srgb[i] = rgb[i] * 12.92;
        else
          srgb[i] = 1.055 * std::pow(rgb[i], 1. / 2.4) - 0.055;
      }
  }

  // Lab (D50) -> XYZ -> rgb (D65) -> sRGB (D65)
  void lab2srgb(const Vec3& Lab, Vec3& srgb) const
  {
    Vec3 XYZ;
    lab2xyz(Lab, XYZ);
    Vec3 rgb;
    xyz2rgb(XYZ, rgb);
    rgb2srgb(rgb, srgb);
  }

  // sRGB (D65) -> rgb (D65) -> XYZ -> Lab
  void srgb2lab(const Vec3& srgb, Vec3& Lab) const
  {
    Vec3 rgb;
    srgb2rgb(srgb, rgb);
    Vec3 XYZ;
    rgb2xyz(rgb, XYZ);
    xyz2lab(XYZ, Lab);
  }

  // Lab  -> XYZ -> rgb (D65)
  void lab2rgb(const Vec3& Lab, Vec3& rgb) const
  {
    Vec3 XYZ;
    lab2xyz(Lab, XYZ);
    xyz2rgb(XYZ, rgb);
  }

  // rgb (D65) -> XYZ -> Lab
  void rgb2lab(const Vec3& rgb, Vec3& Lab) const
  {
    Vec3 XYZ;
    rgb2xyz(rgb, XYZ);
    xyz2lab(XYZ, Lab);
  }

  // XYZ -> rgb (D65) -> sRGB (D65)
  void xyz2srgb(const Vec3& XYZ, Vec3& srgb) const
  {
    Vec3 rgb;
    xyz2rgb(XYZ, rgb);
    rgb2srgb(rgb, srgb);
  }

  // [0..1] -> [0..1]
  void hsv2srgb(const Vec3& hsv, Vec3& srgb) const
  {
    const Scalar  h  = (360. * hsv[0]) / 60.;
    const int32_t hi = static_cast<int32_t>(std::floor(h));
    Scalar        f  = (h - hi);

    Scalar p = hsv[2] * (1 - hsv[1]);
    Scalar q = hsv[2] * (1 - hsv[1] * f);
    Scalar t = hsv[2] * (1 - hsv[1] * (1 - f));

    if (hi == 1)
      {
        srgb[0] = q;
        srgb[1] = hsv[2];
        srgb[2] = p;
      }
    else if (hi == 2)
      {
        srgb[0] = p;
        srgb[1] = hsv[2];
        srgb[2] = t;
      }
    else if (hi == 3)
      {
        srgb[0] = p;
        srgb[1] = q;
        srgb[2] = hsv[2];
      }
    else if (hi == 4)
      {
        srgb[0] = t;
        srgb[1] = p;
        srgb[2] = hsv[2];
      }
    else if (hi == 5)
      {
        srgb[0] = hsv[2];
        srgb[1] = p;
        srgb[2] = q;
      }
    else
      {
        srgb[0] = hsv[2];
        srgb[1] = t;
        srgb[2] = p;
      }
  }

  // [0..1] -> [0..1]
  void srgb2hsv(const Vec3& srgb, Vec3& hsv) const
  {
    Scalar min;
    Scalar max;
    Scalar delMax;

    min    = std::min<Scalar>(std::min<Scalar>(srgb[0], srgb[1]), srgb[2]);
    max    = std::max<Scalar>(std::max<Scalar>(srgb[0], srgb[1]), srgb[2]);
    delMax = 1. / (max - min);

    const Scalar fa = 1. / 360.0;

    if (fuzzyCompare(max, min))
      hsv[0] = 0;
    else if (fuzzyCompare(max, srgb[0]))
      hsv[0] = 60.0 * (0 + (srgb[1] - srgb[2]) * delMax);
    else if (fuzzyCompare(max, srgb[1]))
      hsv[0] = 60.0 * (2 + (srgb[2] - srgb[0]) * delMax);
    else if (fuzzyCompare(max, srgb[2]))
      hsv[0] = 60.0 * (4 + (srgb[0] - srgb[1]) * delMax);

    if (hsv[0] < 0.0)
      hsv[0] += 360.0;

    if (fuzzyCompare(max, 0.0))
      {
        hsv[1] = 0.0;
      }
    else
      {
        hsv[1] = (max - min) / max;
      }
    hsv[2] = max;

    hsv[0] *= fa;
  }

  // sRGB (D65) -> XYZ
  void srgb2xyz(const Vec3& srgb, Vec3& XYZ) const
  {
    Vec3 rgb;
    srgb2rgb(srgb, rgb);
    rgb2xyz(rgb, XYZ);
  }

  void xyz2xyY(const Vec3& XYZ, Vec3& xyY) const
  {
    xyY[0] = XYZ[0] / (XYZ[0] + XYZ[1] + XYZ[2]);
    xyY[1] = XYZ[1] / (XYZ[0] + XYZ[1] + XYZ[2]);
    xyY[2] = XYZ[1];
  }

  void xyY2xyz(const Vec3& xyY, Vec3& XYZ) const
  {
    XYZ[1] = xyY[2];
    XYZ[0] = (xyY[2] / xyY[1]) * xyY[0];
    XYZ[2] = (xyY[2] / xyY[1]) * (1 - xyY[0] - xyY[1]);
  }

  void Luv2XYZ(const Vec3& Luv, Vec3& XYZ) const
  {
    const Scalar eps  = 216. / 24389.;
    const Scalar k    = 24389. / 27.;
    const Scalar keps = k * eps;

    XYZ[1] =
      (Luv[0] > keps) ? (std::pow((Luv[0] + 16.) / 116., 3.)) : (Luv[1] / k);

    Scalar Xr, Yr, Zr;
    Xr = illuminant[0];
    Yr = illuminant[1];
    Zr = illuminant[2];

    Scalar u0, v0;
    u0 = (4. * Xr) / (Xr + 15. * Yr + 3. * Zr);
    v0 = (9. * Yr) / (Xr + 15. * Yr + 3. * Zr);

    Scalar a, b, c, d;
    a = (1. / 3.) * (((52. * Luv[0]) / (Luv[1] + 13. * Luv[0] * u0)) - 1.);
    b = -5. * XYZ[1];
    c = -(1. / 3.);
    d = XYZ[1] * (((39. * Luv[0]) / (Luv[2] + 13. * Luv[0] * v0)) - 5.);

    XYZ[0] = (d - b) / (a - c);
    XYZ[2] = XYZ[0] * a + b;
  }

  void Yuv2rgb(const Vec3& Yuv, Vec3& rgb) const
  {
    rgb[2] = 1.164 * (Yuv[0] - 16) + 2.018 * (Yuv[1] - 128.);
    rgb[1] =
      1.164 * (Yuv[0] - 16) - 0.813 * (Yuv[2] - 128.) - 0.391 * (Yuv[1] - 128.);
    rgb[0] = 1.164 * (Yuv[0] - 16) + 1.596 * (Yuv[2] - 128.);

    const Scalar s = (1. / 255.);
    rgb[0] *= s;
    rgb[1] *= s;
    rgb[2] *= s;
  }

  void rgb2Yuv(const Vec3& rgb, Vec3& Yuv) const
  {
    Vec3 rgb_scaled;
    rgb_scaled[0] = rgb[0] * 255.;
    rgb_scaled[1] = rgb[1] * 255.;
    rgb_scaled[2] = rgb[2] * 255.;
    Yuv[0]        = (0.257 * rgb_scaled[0]) + (0.504 * rgb_scaled[1]) +
             (0.098 * rgb_scaled[2]) + 16;
    Yuv[2] = (0.439 * rgb_scaled[0]) - (0.368 * rgb_scaled[1]) -
             (0.071 * rgb_scaled[2]) + 128;
    Yuv[1] = -(0.148 * rgb_scaled[0]) - (0.291 * rgb_scaled[1]) +
             (0.439 * rgb_scaled[2]) + 128;
  }

  void XYZ2Luv(const Vec3& XYZ, Vec3& Luv) const
  {
    const Scalar eps = 216. / 24389.;
    const Scalar k   = 24389. / 27.;

    // chromatic adaption, reference white
    Scalar Xr = illuminant[0];
    Scalar Yr = illuminant[1];
    Scalar Zr = illuminant[2];

    Scalar yr = XYZ[1] / Yr;

    Luv[0] = (yr > eps) ? (116. * std::pow(yr, 1. / 3.) - 16.) : k * yr;

    Scalar nen = XYZ[0] + 15. * XYZ[1] + 3 * XYZ[2];
    Scalar u_  = (4 * XYZ[0]) / (nen);
    Scalar v_  = (9 * XYZ[1]) / (nen);
    nen        = Xr + 15. * Yr + 3 * Zr;
    Scalar ur_ = (4 * Xr) / (nen);
    Scalar vr_ = (9 * Yr) / (nen);

    Luv[1] = 13. * Luv[0] * (u_ - ur_);
    Luv[2] = 13. * Luv[0] * (v_ - vr_);
  }

  void Luv2LCHuv(const Vec3& Luv, Vec3& LCHuv) const
  {
    LCHuv[0] = Luv[0];
    LCHuv[1] = std::sqrt((Luv[1] * Luv[1]) + (Luv[2] * Luv[2]));
    LCHuv[2] = atan2(Luv[2], Luv[1]);
  }

  void LCHuv2Luv(const Vec3& LCHuv, Vec3& Luv) const
  {
    Luv[0] = LCHuv[0];
    Luv[1] = LCHuv[1] * cos(LCHuv[2]);
    Luv[2] = LCHuv[1] * sin(LCHuv[2]);
  }

  void srgb2CIELCHab(const Vec3& srgb, Vec3& CIELCHab) const
  {
    Vec3 v;
    srgb2lab(srgb, v);
    lab2LCHab(v, CIELCHab);
  }
  void rgb2CIELCHab(const Vec3& rgb, Vec3& CIELCHab) const
  {
    Vec3 v;
    rgb2lab(rgb, v);
    lab2LCHab(v, CIELCHab);
  }
  void CIELCHab2srgb(const Vec3& LCHab, Vec3& srgb) const
  {
    Vec3 v;
    LCHab2lab(LCHab, v);
    lab2srgb(v, srgb);
  }
  void CIELCHab2rgb(const Vec3& CIELCHab, Vec3& rgb) const
  {
    Vec3 v;
    LCHab2lab(CIELCHab, v);
    lab2rgb(v, rgb);
  }

  /**
   * @brief Compute difference of two given colors.
   *
   * Ported to C++ from matlab script
   * (https://github.com/scienstanford/iqmetrics/blob/master/SCIELAB_FR/deltaE2000.m)
   *
   * Remarks:
   * Refer to https://en.wikipedia.org/wiki/Color_difference for details.
   *
   * @param lab1
   * @param lab2
   * @return Scalar
   */
  static Scalar ColorDifferenceCIEDE2000(const Vec3& lab1, const Vec3& lab2)
  {
    Scalar Lstd = lab1[0];
    Scalar astd = lab1[1];
    Scalar bstd = lab1[2];

    Scalar Lsample = lab2[0];
    Scalar asample = lab2[1];
    Scalar bsample = lab2[2];

    Scalar Cabstd    = std::sqrt(astd * astd + bstd * bstd);
    Scalar Cabsample = std::sqrt(asample * asample + bsample * bsample);

    Scalar Cabarithmean = (Cabstd + Cabsample) / 2.;

    Scalar G =
      0.5f * (1. - std::sqrt(std::pow(Cabarithmean, 7.) /
                             (std::pow(Cabarithmean, 7.) + std::pow(25., 7.))));

    Scalar apstd    = (1. + G) * astd;    // aprime in paper
    Scalar apsample = (1. + G) * asample; // aprime in paper
    Scalar Cpsample = std::sqrt(apsample * apsample + bsample * bsample);

    Scalar Cpstd = std::sqrt(apstd * apstd + bstd * bstd);
    // Compute product of chromas
    Scalar Cpprod = (Cpsample * Cpstd);

    // Ensure hue is between 0 and 2pi
    Scalar hpstd = std::atan2(bstd, apstd);
    if (hpstd < 0)
      hpstd += 2. * Pi; // rollover ones that come -ve

    Scalar hpsample = std::atan2(bsample, apsample);
    if (hpsample < 0)
      hpsample += 2. * Pi;
    if (fuzzyCompare((fabs(apsample) + fabs(bsample)), 0.))
      hpsample = 0.;

    Scalar dL = (Lsample - Lstd);
    Scalar dC = (Cpsample - Cpstd);

    // Computation of hue difference
    Scalar dhp = (hpsample - hpstd);
    if (dhp > Pi)
      dhp -= 2. * Pi;
    if (dhp < -Pi)
      dhp += 2. * Pi;
    // set chroma difference to zero if the product of chromas is zero
    if (fuzzyCompare(Cpprod, 0.))
      dhp = 0.;

    // Note that the defining equations actually need
    // signed Hue and chroma differences which is different
    // from prior color difference formulae

    Scalar dH = 2. * std::sqrt(Cpprod) * sin(dhp / 2.);
    //%dH2 = 4*Cpprod.*(sin(dhp/2)).^2;

    // weighting functions
    Scalar Lp = (Lsample + Lstd) / 2.;
    Scalar Cp = (Cpstd + Cpsample) / 2.;

    // Average Hue Computation
    // This is equivalent to that in the paper but simpler programmatically.
    // Note average hue is computed in radians and converted to degrees only
    // where needed
    Scalar hp = (hpstd + hpsample) / 2.;
    // Identify positions for which abs hue diff exceeds 180 degrees
    if (fabs(hpstd - hpsample) > Pi)
      hp -= Pi;
    // rollover ones that come -ve
    if (hp < 0)
      hp += 2. * Pi;

    // Check if one of the chroma values is zero, in which case set
    // mean hue to the sum which is equivalent to other value
    if (fuzzyCompare(Cpprod, 0.))
      hp = hpsample + hpstd;

    Scalar Lpm502 = (Lp - 50.) * (Lp - 50.);
    Scalar Sl     = 1. + 0.015f * Lpm502 / std::sqrt(20.0f + Lpm502);
    Scalar Sc     = 1. + 0.045f * Cp;
    Scalar Ta     = 1. - 0.17f * std::cos(hp - Pi / 6.) +
                0.24f * std::cos(2. * hp) +
                0.32f * std::cos(3. * hp + Pi / 30.) -
                0.20f * std::cos(4. * hp - 63. * Pi / 180.);
    Scalar Sh = 1. + 0.015f * Cp * Ta;
    Scalar delthetarad =
      (30. * Pi / 180.) *
      std::exp(-std::pow(((180. / Pi * hp - 275.) / 25.), 2.));
    Scalar Rc =
      2. * std::sqrt(std::pow(Cp, 7.) / (std::pow(Cp, 7.) + std::pow(25., 7.)));
    Scalar RT = -std::sin(2.0f * delthetarad) * Rc;

    // The CIE 00 color difference
    return std::sqrt(std::pow((dL / Sl), 2.) + std::pow((dC / Sc), 2.) +
                     std::pow((dH / Sh), 2.) + RT * (dC / Sc) * (dH / Sh));
  }

private:
  const std::array<Scalar, 3U> illuminant;

  static Scalar f(Scalar t)
  {
    return (t > std::pow<Scalar>(6. / 29., 3.))
             ? std::pow<Scalar>(t, 1. / 3.)
             : (1. / 3.) * std::pow<Scalar>(29. / 6., 2.) * t + (4. / 29.);
  }

  static Scalar fi(Scalar t)
  {
    return (t > 6. / 29.)
             ? std::pow<Scalar>(t, 3.)
             : 3. * std::pow<Scalar>(6. / 29., 2.) * (t - (4. / 29.));
  }

  /**
   * @brief Simple fuzzyCompare comparison of floating point values.
   *
   * @return true
   * @return false
   */
  static bool fuzzyCompare(const Scalar a, const Scalar b)
  {
    constexpr auto FuzzynessFactor = 10.0;
    return std::fabs(a - b) <
           (FuzzynessFactor * std::numeric_limits<Scalar>::epsilon());
  }
};

} // namespace color

#endif // COLOR_CONVERTER_HXX
