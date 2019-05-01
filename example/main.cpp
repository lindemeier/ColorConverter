

#include <iostream>

#include <ColorConverter/ColorConverter.hxx>

int main()
{
  // define custom vec type that is defined by our type trait above
  using vec = std::array<double, 3U>;
  // Save some typing
  using Converter = color::ColorConverter<double, std::array>;

  // initialize the converter to standard sRGB illuminant (which is also the
  // default argument).
  Converter converter(Converter::IM_ILLUMINANT_D65);

  // stating values in linear RGB
  vec linearRgb = {0.2, 0.3, 0.4};
  vec sRgb;

  // convert linear to sRGB
  converter.rgb2srgb(linearRgb, sRgb);

  // convert the sRGB value to CIELab
  vec cieLab;
  converter.srgb2lab(sRgb, cieLab);

  // convert back to linear RGB from CIELab
  vec convertedLinearRgb;
  converter.lab2rgb(cieLab, convertedLinearRgb);
  // cpmpute the difference of input and converted linear RGB
  const auto differenceEuclid =
    std::sqrt(std::pow(convertedLinearRgb[0U] - linearRgb[0U], 2.0) +
              std::pow(convertedLinearRgb[1U] - linearRgb[1U], 2.0) +
              std::pow(convertedLinearRgb[2U] - linearRgb[2U], 2.0));
  // check if our converter did a good job
  if (differenceEuclid < 0.0000001)
    {
      std::cout << "it worked!";
    }
  else
    {
      std::cout << "it didn't work!";
    }
  std::cout << std::endl;

  return 0;
}
