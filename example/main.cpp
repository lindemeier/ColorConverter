#include <ColorConverter/ColorConverter.hxx>

#include <iostream>

int main()
{
  using vec = color::ColorConverter<double>::vec3;

  color::ColorConverter<double> converter;

  vec linearRgb = {0.2, 0.3, 0.4};
  vec sRgb;

  converter.convert(linearRgb, sRgb,
                    color::ColorConverter<double>::Conversion::rgb_2_srgb);

  vec cieLab;
  converter.convert(sRgb, cieLab,
                    color::ColorConverter<double>::Conversion::srgb_2_CIELab);

  vec convertedLinearRgb;
  converter.convert(cieLab, convertedLinearRgb,
                    color::ColorConverter<double>::Conversion::CIELab_2_rgb);
  const auto differenceEuclid =
    std::sqrt(std::pow(convertedLinearRgb[0U] - linearRgb[0U], 2.0) +
              std::pow(convertedLinearRgb[1U] - linearRgb[1U], 2.0) +
              std::pow(convertedLinearRgb[2U] - linearRgb[2U], 2.0));
  if (differenceEuclid < 0.0000001)
    {
      std::cout << "it worked!";
    }
  else
    {
      std::cout << "it didn't work!";
    }

  return 0;
}
