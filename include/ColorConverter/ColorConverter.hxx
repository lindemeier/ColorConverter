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

namespace color
{

/**
 * @brief Type trait that can be used for defining new types that should be
 * supported by the ColorConverter.
 *
 * @tparam T the Scalar
 */
template <typename T>
class VecType
{
public:
  using Scalar = T;
  static_assert(std::is_floating_point<T>::value,
                "only floating point scalar types supported.");

  static std::string getName() { return typeid(T).name(); }
};

/**
 * @brief Class that holds color cnoversion functions with reference to a given
 * white point (illuminant).
 *
 * Most of the conversion and constants were taken from
 * http://brucelindbloom.com/index.html (last accessed 2019 April 14th)
 * Credits to that guy!
 *
 * @tparam T precision
 */
template <class Vec3>
class ColorConverter
{
public:
  using Scalar = typename VecType<Vec3>::Scalar;

  /**
   * @brief Represents conversion source2target.
   */
  enum class Conversion
  {
    CIELab_2_XYZ,
    XYZ_2_CIELab,
    CIELab_2_LCHab, // CIE LCHab
    LCH_ab_2_CIELab,
    ryb_2_rgb, // http://en.wikipedia.org/wiki/RYB_color_model,
               // http://threekings.tk/mirror/ryb_TR.pdf
    rgb_2_cmy,
    cmy_2_rgb,
    rgb_2_XYZ,       // uses sRGB chromatic adapted matrix
    XYZ_2_rgb,       // uses sRGB chromatic adapted matrix
    srgb_2_rgb,      // make linear rgb, no chromatic adaption
    rgb_2_srgb,      // make sRGB, with gamma
    CIELab_2_srgb,   // Lab (whitepoint) -> XYZ -> rgb (D65) -> sRGB (D65)
    srgb_2_CIELab,   // sRGB (D65) -> rgb (D65) -> XYZ -> Lab (whitepoint)
    srgb_2_CIELCHab, // sRGB (D65) -> rgb (D65) -> XYZ -> Lab (whitepoint)
    CIELCHab_2_srgb, // sRGB (D65) -> rgb (D65) -> XYZ -> Lab (whitepoint)
    rgb_2_CIELCHab,
    CIELCHab_2_rgb,
    CIELab_2_rgb, // Lab (whitepoint) -> XYZ -> rgb (D65)
    rgb_2_CIELab, // rgb (D65) -> XYZ -> Lab (D50)
    XYZ_2_srgb,   // XYZ -> rgb (D65) -> sRGB (D65)
    hsv_2_srgb,   // [0..1] -> [0..1]
    srgb_2_hsv,   // [0..1] -> [0..1]
    srgb_2_XYZ,   // sRGB (D65) -> XYZ
    XYZ_2_xyY,
    xyY_2_XYZ,
    Luv_2_XYZ,
    XYZ_2_Luv,
    Luv_2_LCH_uv, // CIE LCHuv
    LCH_uv_2_Luv,
    Yuv_2_rgb, // Y Cb Cr
    rgb_2_Yuv
  };

  static constexpr Vec3 IM_ILLUMINANT_A   = {1.09850, 1.00000, 0.35585};
  static constexpr Vec3 IM_ILLUMINANT_B   = {0.99072, 1.00000, 0.85223};
  static constexpr Vec3 IM_ILLUMINANT_C   = {0.98074, 1.00000, 1.18232};
  static constexpr Vec3 IM_ILLUMINANT_D50 = {0.96422, 1.00000, 0.82521};
  static constexpr Vec3 IM_ILLUMINANT_D55 = {0.95682, 1.00000, 0.92149};
  static constexpr Vec3 IM_ILLUMINANT_D65 = {0.95047, 1.00000, 1.08883};
  static constexpr Vec3 IM_ILLUMINANT_D75 = {0.94972, 1.00000, 1.22638};
  static constexpr Vec3 IM_ILLUMINANT_E   = {1.00000, 1.00000, 1.00000};
  static constexpr Vec3 IM_ILLUMINANT_2   = {0.99186, 1.00000, 0.67393};
  static constexpr Vec3 IM_ILLUMINANT_7   = {0.95041, 1.00000, 1.08747};
  static constexpr Vec3 IM_ILLUMINANT_11  = {1.00962, 1.00000, 0.64350};

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
  ColorConverter(const Vec3& whitepoint = IM_ILLUMINANT_D65);

  inline void convert(const Vec3& input, Vec3& output,
                      Conversion conversion) const;

private:
  const Vec3 illuminant;

  inline void lab2xyz(const Vec3& Lab, Vec3& XYZ) const;
  inline void xyz2lab(const Vec3& XYZ, Vec3& Lab) const;
  inline void lab2LCHab(const Vec3& Lab, Vec3& LCHab) const;
  inline void LCHab2lab(const Vec3& LCHab, Vec3& Lab) const;
  inline void ryb2rgb(const Vec3& ryb, Vec3& rgb) const;
  inline void rgb2cmy(const Vec3& rgb, Vec3& cmy) const;
  inline void cmy2rgb(const Vec3& cmy, Vec3& rgb) const;
  inline void rgb2xyz(const Vec3& rgb, Vec3& XYZ) const;
  inline void xyz2rgb(const Vec3& XYZ, Vec3& rgb) const;
  inline void srgb2rgb(const Vec3& srgb, Vec3& rgb) const;
  inline void rgb2srgb(const Vec3& rgb, Vec3& srgb) const;
  inline void lab2srgb(const Vec3& Lab, Vec3& srgb) const;
  inline void srgb2lab(const Vec3& srgb, Vec3& Lab) const;
  inline void lab2rgb(const Vec3& Lab, Vec3& rgb) const;
  inline void rgb2lab(const Vec3& rgb, Vec3& Lab) const;
  inline void xyz2srgb(const Vec3& XYZ, Vec3& srgb) const;
  inline void hsv2srgb(const Vec3& hsv, Vec3& srgb) const;
  inline void srgb2hsv(const Vec3& srgb, Vec3& hsv) const;
  inline void srgb2xyz(const Vec3& srgb, Vec3& XYZ) const;
  inline void xyz2xyY(const Vec3& XYZ, Vec3& xyY) const;
  inline void xyY2xyz(const Vec3& xyY, Vec3& XYZ) const;
  inline void Luv2XYZ(const Vec3& Luv, Vec3& XYZ) const;
  inline void XYZ2Luv(const Vec3& XYZ, Vec3& Luv) const;
  inline void Luv2LCHuv(const Vec3& Luv, Vec3& LCHuv) const;
  inline void LCHuv2Luv(const Vec3& LCHuv, Vec3& Luv) const;
  inline void Yuv2rgb(const Vec3& Yuv, Vec3& rgb) const;
  inline void rgb2Yuv(const Vec3& rgb, Vec3& Yuv) const;

  static bool inline fuzzy(const Scalar a, const Scalar b);
};

template <class Vec3>
ColorConverter<Vec3>::ColorConverter(const Vec3& whitepoint)
  : illuminant(whitepoint)
{
}

template <class Vec3>
void ColorConverter<Vec3>::convert(const Vec3& input, Vec3& output,
                                   Conversion conversion) const
{
  switch (conversion)
    {
    case Conversion::CIELab_2_XYZ:
      {
        lab2xyz(input, output);
        break;
      }
    case Conversion::XYZ_2_CIELab:
      {
        xyz2lab(input, output);
        break;
      }
    case Conversion::CIELab_2_LCHab:
      {
        lab2LCHab(input, output);
        break;
      }
    case Conversion::LCH_ab_2_CIELab:
      {
        LCHab2lab(input, output);
        break;
      }
    case Conversion::ryb_2_rgb:
      {
        ryb2rgb(input, output);
        break;
      }
    case Conversion::rgb_2_cmy:
      {
        rgb2cmy(input, output);
        break;
      }
    case Conversion::cmy_2_rgb:
      {
        cmy2rgb(input, output);
        break;
      }
    case Conversion::rgb_2_XYZ:
      {
        rgb2xyz(input, output);
        break;
      }
    case Conversion::XYZ_2_rgb:
      {
        xyz2rgb(input, output);
        break;
      }
    case Conversion::srgb_2_rgb:
      {
        srgb2rgb(input, output);
        break;
      }
    case Conversion::rgb_2_srgb:
      {
        rgb2srgb(input, output);
        break;
      }
    case Conversion::CIELab_2_srgb:
      {
        lab2srgb(input, output);
        break;
      }
    case Conversion::srgb_2_CIELab:
      {
        srgb2lab(input, output);
        break;
      }
    case Conversion::CIELab_2_rgb:
      {
        lab2rgb(input, output);
        break;
      }
    case Conversion::rgb_2_CIELab:
      {
        rgb2lab(input, output);
        break;
      }
    case Conversion::XYZ_2_srgb:
      {
        xyz2srgb(input, output);
        break;
      }
    case Conversion::hsv_2_srgb:
      {
        hsv2srgb(input, output);
        break;
      }
    case Conversion::srgb_2_hsv:
      {
        srgb2hsv(input, output);
        break;
      }
    case Conversion::srgb_2_XYZ:
      {
        srgb2xyz(input, output);
        break;
      }
    case Conversion::XYZ_2_xyY:
      {
        xyz2xyY(input, output);
        break;
      }
    case Conversion::xyY_2_XYZ:
      {
        xyY2xyz(input, output);
        break;
      }
    case Conversion::Luv_2_XYZ:
      {
        Luv2XYZ(input, output);
        break;
      }
    case Conversion::XYZ_2_Luv:
      {
        XYZ2Luv(input, output);
        break;
      }
    case Conversion::Luv_2_LCH_uv:
      {
        Luv2LCHuv(input, output);
        break;
      }
    case Conversion::LCH_uv_2_Luv:
      {
        LCHuv2Luv(input, output);
        break;
      }
    case Conversion::Yuv_2_rgb:
      {
        Yuv2rgb(input, output);
        break;
      }
    case Conversion::rgb_2_Yuv:
      {
        rgb2Yuv(input, output);
        break;
      }
    case Conversion::srgb_2_CIELCHab:
      {
        Vec3 v;
        srgb2lab(input, v);
        lab2LCHab(v, output);
        break;
      }
    case Conversion::rgb_2_CIELCHab:
      {
        Vec3 v;
        rgb2lab(input, v);
        lab2LCHab(v, output);
        break;
      }
    case Conversion::CIELCHab_2_srgb:
      {
        Vec3 v;
        LCHab2lab(input, v);
        lab2srgb(v, output);
        break;
      }
    case Conversion::CIELCHab_2_rgb:
      {
        Vec3 v;
        LCHab2lab(input, v);
        lab2rgb(v, output);
        break;
      }
    default: break;
    }
}

// rgb to cmy
template <class Vec3>
void ColorConverter<Vec3>::rgb2cmy(const Vec3& rgb, Vec3& cmy) const
{
  cmy[0] = 1. - rgb[0];
  cmy[1] = 1. - rgb[1];
  cmy[2] = 1. - rgb[2];
}

// cmy to rgb
template <class Vec3>
void ColorConverter<Vec3>::cmy2rgb(const Vec3& cmy, Vec3& rgb) const
{
  rgb[0] = 1. - cmy[0];
  rgb[1] = 1. - cmy[1];
  rgb[2] = 1. - cmy[2];
}

// http://en.wikipedia.org/wiki/RYB_color_model
// http://threekings.tk/mirror/ryb_TR.pdf
template <class Vec3>
void ColorConverter<Vec3>::ryb2rgb(const Vec3& ryb, Vec3& rgb) const
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

template <class T>
static inline T f(T t)
{
  return (t > std::pow<T>(6. / 29., 3.))
           ? std::pow<T>(t, 1. / 3.)
           : (1. / 3.) * std::pow<T>(29. / 6., 2.) * t + (4. / 29.);
}

template <class T>
static inline T fi(T t)
{
  return (t > 6. / 29.) ? std::pow<T>(t, 3.)
                        : 3. * std::pow<T>(6. / 29., 2.) * (t - (4. / 29.));
}

template <class Vec3>
void ColorConverter<Vec3>::xyz2lab(const Vec3& XYZ, Vec3& Lab) const
{
  Lab[0] = 116. * f(XYZ[1] / illuminant[1]) - 16.;
  Lab[1] = 500. * (f(XYZ[0] / illuminant[0]) - f(XYZ[1] / illuminant[1]));
  Lab[2] = 200. * (f(XYZ[1] / illuminant[1]) - f(XYZ[2] / illuminant[2]));
}

template <class Vec3>
void ColorConverter<Vec3>::lab2xyz(const Vec3& Lab, Vec3& XYZ) const
{
  // chromatic adaption, reference white
  XYZ[1] = illuminant[1] * fi((1. / 116.) * (Lab[0] + 16.)); // Y
  XYZ[0] = illuminant[0] *
           fi((1. / 116.) * (Lab[0] + 16.) + (1. / 500.) * Lab[1]); // X
  XYZ[2] = illuminant[2] *
           fi((1. / 116.) * (Lab[0] + 16.) - (1. / 200.) * Lab[2]); // Z
}

template <class Vec3>
void ColorConverter<Vec3>::lab2LCHab(const Vec3& Lab, Vec3& LCHab) const
{
  LCHab[0] = Lab[0];                                       // [0,100]
  LCHab[1] = std::sqrt(Lab[1] * Lab[1] + Lab[2] * Lab[2]); // [0,100]

  LCHab[2] = std::atan2(Lab[2], Lab[1]);
  if (LCHab[2] < 0)
    {
      LCHab[2] += M_PI * 2.; // [0, 2pi]
    }
}

template <class Vec3>
void ColorConverter<Vec3>::LCHab2lab(const Vec3& LCHab, Vec3& Lab) const
{
  Lab[0]   = LCHab[0];
  Scalar h = LCHab[2];
  if (h > M_PI)
    {
      h -= M_PI * 2.; // [0, 2pi]
    }
  Lab[1] = LCHab[1] * std::cos(h);
  Lab[2] = LCHab[1] * std::sin(h);
}

// uses sRGB chromatic adapted matrix
template <class Vec3>
void ColorConverter<Vec3>::rgb2xyz(const Vec3& rgb, Vec3& XYZ) const
{
  XYZ[0] = RGB2XYZ_MATRIX[0][0] * rgb[0] + RGB2XYZ_MATRIX[0][1] * rgb[1] +
           RGB2XYZ_MATRIX[0][2] * rgb[2];
  XYZ[1] = RGB2XYZ_MATRIX[1][0] * rgb[0] + RGB2XYZ_MATRIX[1][1] * rgb[1] +
           RGB2XYZ_MATRIX[1][2] * rgb[2];
  XYZ[2] = RGB2XYZ_MATRIX[2][0] * rgb[0] + RGB2XYZ_MATRIX[2][1] * rgb[1] +
           RGB2XYZ_MATRIX[2][2] * rgb[2];
}

// uses sRGB chromatic adapted matrix
template <class Vec3>
void ColorConverter<Vec3>::xyz2rgb(const Vec3& XYZ, Vec3& rgb) const
{
  rgb[0] = XYZ2RGB_MATRIX[0][0] * XYZ[0] + XYZ2RGB_MATRIX[0][1] * XYZ[1] +
           XYZ2RGB_MATRIX[0][2] * XYZ[2];
  rgb[1] = XYZ2RGB_MATRIX[1][0] * XYZ[0] + XYZ2RGB_MATRIX[1][1] * XYZ[1] +
           XYZ2RGB_MATRIX[1][2] * XYZ[2];
  rgb[2] = XYZ2RGB_MATRIX[2][0] * XYZ[0] + XYZ2RGB_MATRIX[2][1] * XYZ[1] +
           XYZ2RGB_MATRIX[2][2] * XYZ[2];
}

// make linear rgb, no chromatic adaption
template <class Vec3>
void ColorConverter<Vec3>::srgb2rgb(const Vec3& srgb, Vec3& rgb) const
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
template <class Vec3>
void ColorConverter<Vec3>::rgb2srgb(const Vec3& rgb, Vec3& srgb) const
{
  for (auto i = 0U; i < 3U; ++i)
    {
      if (rgb[i] <= 0.0031308)
        srgb[i] = rgb[i] * 12.92;
      else
        srgb[i] = 1.055f * std::pow(rgb[i], 1. / 2.4) - 0.055;
    }
}

// Lab (D50) -> XYZ -> rgb (D65) -> sRGB (D65)
template <class Vec3>
void ColorConverter<Vec3>::lab2srgb(const Vec3& Lab, Vec3& srgb) const
{
  Vec3 XYZ;
  lab2xyz(Lab, XYZ);
  Vec3 rgb;
  xyz2rgb(XYZ, rgb);
  rgb2srgb(rgb, srgb);
}

// sRGB (D65) -> rgb (D65) -> XYZ -> Lab
template <class Vec3>
void ColorConverter<Vec3>::srgb2lab(const Vec3& srgb, Vec3& Lab) const
{
  Vec3 rgb;
  srgb2rgb(srgb, rgb);
  Vec3 XYZ;
  rgb2xyz(rgb, XYZ);
  xyz2lab(XYZ, Lab);
}

// Lab  -> XYZ -> rgb (D65)
template <class Vec3>
void ColorConverter<Vec3>::lab2rgb(const Vec3& Lab, Vec3& rgb) const
{
  Vec3 XYZ;
  lab2xyz(Lab, XYZ);
  xyz2rgb(XYZ, rgb);
}

// rgb (D65) -> XYZ -> Lab
template <class Vec3>
void ColorConverter<Vec3>::rgb2lab(const Vec3& rgb, Vec3& Lab) const
{
  Vec3 XYZ;
  rgb2xyz(rgb, XYZ);
  xyz2lab(XYZ, Lab);
}

// XYZ -> rgb (D65) -> sRGB (D65)
template <class Vec3>
void ColorConverter<Vec3>::xyz2srgb(const Vec3& XYZ, Vec3& srgb) const
{
  Vec3 rgb;
  xyz2rgb(XYZ, rgb);
  rgb2srgb(rgb, srgb);
}

// [0..1] -> [0..1]
template <class Vec3>
void ColorConverter<Vec3>::hsv2srgb(const Vec3& hsv, Vec3& srgb) const
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
template <class Vec3>
void ColorConverter<Vec3>::srgb2hsv(const Vec3& srgb, Vec3& hsv) const
{
  Scalar min;
  Scalar max;
  Scalar delMax;

  min    = std::min<Scalar>(std::min<Scalar>(srgb[0], srgb[1]), srgb[2]);
  max    = std::max<Scalar>(std::max<Scalar>(srgb[0], srgb[1]), srgb[2]);
  delMax = 1. / (max - min);

  const Scalar fa = 1. / 360.0;

  if (fuzzy(max, min))
    hsv[0] = 0;
  else if (fuzzy(max, srgb[0]))
    hsv[0] = 60.0 * (0 + (srgb[1] - srgb[2]) * delMax);
  else if (fuzzy(max, srgb[1]))
    hsv[0] = 60.0 * (2 + (srgb[2] - srgb[0]) * delMax);
  else if (fuzzy(max, srgb[2]))
    hsv[0] = 60.0 * (4 + (srgb[0] - srgb[1]) * delMax);

  if (hsv[0] < 0.0)
    hsv[0] += 360.0;

  if (fuzzy(max, 0.0))
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
template <class Vec3>
void ColorConverter<Vec3>::srgb2xyz(const Vec3& srgb, Vec3& XYZ) const
{
  Vec3 rgb;
  srgb2rgb(srgb, rgb);
  rgb2xyz(rgb, XYZ);
}

template <class Vec3>
void ColorConverter<Vec3>::xyz2xyY(const Vec3& XYZ, Vec3& xyY) const
{
  xyY[0] = XYZ[0] / (XYZ[0] + XYZ[1] + XYZ[2]);
  xyY[1] = XYZ[1] / (XYZ[0] + XYZ[1] + XYZ[2]);
  xyY[2] = XYZ[1];
}

template <class Vec3>
void ColorConverter<Vec3>::xyY2xyz(const Vec3& xyY, Vec3& XYZ) const
{
  XYZ[1] = xyY[2];
  XYZ[0] = (xyY[2] / xyY[1]) * xyY[0];
  XYZ[2] = (xyY[2] / xyY[1]) * (1 - xyY[0] - xyY[1]);
}

template <class Vec3>
void ColorConverter<Vec3>::Luv2XYZ(const Vec3& Luv, Vec3& XYZ) const
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

template <class Vec3>
void ColorConverter<Vec3>::Yuv2rgb(const Vec3& Yuv, Vec3& rgb) const
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

template <class Vec3>
void ColorConverter<Vec3>::rgb2Yuv(const Vec3& rgb, Vec3& Yuv) const
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

template <class Vec3>
void ColorConverter<Vec3>::XYZ2Luv(const Vec3& XYZ, Vec3& Luv) const
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

template <class Vec3>
void ColorConverter<Vec3>::Luv2LCHuv(const Vec3& Luv, Vec3& LCHuv) const
{
  LCHuv[0] = Luv[0];
  LCHuv[1] = std::sqrt((Luv[1] * Luv[1]) + (Luv[2] * Luv[2]));
  LCHuv[2] = atan2(Luv[2], Luv[1]);
}

template <class Vec3>
void ColorConverter<Vec3>::LCHuv2Luv(const Vec3& LCHuv, Vec3& Luv) const
{
  Luv[0] = LCHuv[0];
  Luv[1] = LCHuv[1] * cos(LCHuv[2]);
  Luv[2] = LCHuv[1] * sin(LCHuv[2]);
}

/**
 * @brief Simple fuzzy comparison of floating point values.
 *
 * @tparam T
 * @param a
 * @param b
 * @return true
 * @return false
 */
template <class Vec3>
bool ColorConverter<Vec3>::fuzzy(const typename VecType<Vec3>::Scalar a,
                                 const typename VecType<Vec3>::Scalar b)
{
  return std::fabs(a - b) <
         std::numeric_limits<typename VecType<Vec3>::Scalar>::epsilon();
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
 * @tparam T
 * @param lab1
 * @param lab2
 * @return T
 */
template <class Vec3>
typename VecType<Vec3>::Scalar
ColorDifferenceCIEDE2000(const typename ColorConverter<Vec3>::Vec3& lab1,
                         const typename ColorConverter<Vec3>::Vec3& lab2)
{
  using Scalar = typename VecType<Vec3>::Scalar;
  Scalar Lstd  = lab1[0];
  Scalar astd  = lab1[1];
  Scalar bstd  = lab1[2];

  Scalar Lsample = lab2[0];
  Scalar asample = lab2[1];
  Scalar bsample = lab2[2];

  constexpr Scalar pi = 3.1415926535897932384626433832795;

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
    hpstd += 2. * pi; // rollover ones that come -ve

  Scalar hpsample = std::atan2(bsample, apsample);
  if (hpsample < 0)
    hpsample += 2. * pi;
  if (fuzzy((fabs(apsample) + fabs(bsample)), 0.))
    hpsample = 0.;

  Scalar dL = (Lsample - Lstd);
  Scalar dC = (Cpsample - Cpstd);

  // Computation of hue difference
  Scalar dhp = (hpsample - hpstd);
  if (dhp > pi)
    dhp -= 2. * pi;
  if (dhp < -pi)
    dhp += 2. * pi;
  // set chroma difference to zero if the product of chromas is zero
  if (fuzzy(Cpprod, 0.))
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
  if (fabs(hpstd - hpsample) > pi)
    hp -= pi;
  // rollover ones that come -ve
  if (hp < 0)
    hp += 2. * pi;

  // Check if one of the chroma values is zero, in which case set
  // mean hue to the sum which is equivalent to other value
  if (fuzzy(Cpprod, 0.))
    hp = hpsample + hpstd;

  Scalar Lpm502 = (Lp - 50.) * (Lp - 50.);
  Scalar Sl     = 1. + 0.015f * Lpm502 / std::sqrt(20.0f + Lpm502);
  Scalar Sc     = 1. + 0.045f * Cp;
  Scalar Ta = 1. - 0.17f * std::cos(hp - pi / 6.) + 0.24f * std::cos(2. * hp) +
              0.32f * std::cos(3. * hp + pi / 30.) -
              0.20f * std::cos(4. * hp - 63. * pi / 180.);
  Scalar Sh          = 1. + 0.015f * Cp * Ta;
  Scalar delthetarad = (30. * pi / 180.) *
                       std::exp(-std::pow(((180. / pi * hp - 275.) / 25.), 2.));
  Scalar Rc =
    2. * std::sqrt(std::pow(Cp, 7.) / (std::pow(Cp, 7.) + std::pow(25., 7.)));
  Scalar RT = -std::sin(2.0f * delthetarad) * Rc;

  // The CIE 00 color difference
  return std::sqrt(std::pow((dL / Sl), 2.) + std::pow((dC / Sc), 2.) +
                   std::pow((dH / Sh), 2.) + RT * (dC / Sc) * (dH / Sh));
}

} // namespace color

#endif // COLOR_CONVERTER_HXX
