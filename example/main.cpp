

#include <iostream>

#include <ColorConverter/ColorConverter.hxx>

namespace color
{
/**
 * @brief Example type trait for std::array
 *
 * @tparam T the scalar floating point value
 */
template <typename T>
class VecType<std::array<T, 3U>>
{
public:
  /**
   * @brief This is needed by the converter to know what the scalr type of your
   * vec type is.
   */
  using Scalar = T;
  static_assert(std::is_floating_point<T>::value,
                "only floating point scalar types supported.");

  static std::string getName() { return typeid(T).name(); }
};

} // namespace color

int main()
{
  // define custom vec type that is defined by our type trait above
  using vec = std::array<double, 3U>;
  // Save some typing
  using Converter = color::ColorConverter<vec>;

  // initialize the converter to standard sRGB illuminant (which is also the
  // default argument).
  Converter converter(Converter::IM_ILLUMINANT_D65);

  // stating values in linear RGB
  vec linearRgb = {0.2, 0.3, 0.4};
  vec sRgb;

  // convert linear to sRGB
  converter.convert(linearRgb, sRgb, Converter::Conversion::rgb_2_srgb);

  // convert the sRGB value to CIELab
  vec cieLab;
  converter.convert(sRgb, cieLab, Converter::Conversion::srgb_2_CIELab);

  // convert back to linear RGB from CIELab
  vec convertedLinearRgb;
  converter.convert(cieLab, convertedLinearRgb,
                    Converter::Conversion::CIELab_2_rgb);
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

  return 0;
}
