#include "polynomial_nttp.cpp"
#include <filesystem>
#include <fstream>
#include <thread>
#include <ostream>
#include <print>
#include <string>
#include <tuple>
#include <type_traits>
#include <variant>

namespace fs = std::filesystem;

const std::string filename{"thousand_divisions"};
const std::string quotient_file_ext{".quotient"};
const std::string remainder_file_ext{".remainder"};

template<int N>
constexpr auto division_test()
{
  constexpr std::array<double, 5> coefficients_a{1, 0, N, 0, 1};
  constexpr math::polynomial_nttp<double, 4> a(std::move(coefficients_a));
  constexpr std::array<double, 3> coefficients_b{1, 0, 1};
  constexpr math::polynomial_nttp<double, 2> b(std::move(coefficients_b));
  constexpr auto a_divided_by_b = math::
                                  division_prototype<double, 4, a, 2, b>();
  return a_divided_by_b;
}

consteval auto divide_thousand_polynomials()
{
  return [&]() {
    using deg_0_poly = math::polynomial_nttp<double, 0>;
    using deg_0_poly_pair = std::pair<int, deg_0_poly>;
    std::array<deg_0_poly_pair, 1000 + 1> possible_deg_0_remainders{};

    using deg_1_poly = math::polynomial_nttp<double, 1>;
    using deg_1_poly_pair = std::pair<int, deg_1_poly>;
    std::array<deg_1_poly_pair, 1000 + 1> possible_deg_1_remainders{};

    using deg_2_poly = math::polynomial_nttp<double, 2>;
    using deg_2_poly_pair = std::pair<int, deg_2_poly>;
    std::array<deg_2_poly_pair, 1000 + 1> possible_deg_2_quotients{};
    std::array<deg_2_poly_pair, 1000 + 1> possible_deg_2_remainders{};

    using deg_3_poly = math::polynomial_nttp<double, 3>;
    using deg_3_poly_pair = std::pair<int, deg_3_poly>;
    std::array<deg_3_poly_pair, 1000 + 1> possible_deg_3_quotients{};

    using deg_4_poly = math::polynomial_nttp<double, 4>;
    using deg_4_poly_pair = std::pair<int, deg_4_poly>;
    std::array<deg_4_poly_pair, 1000 + 1> possible_deg_4_quotients{};

    using deg_unknown_poly = std::variant<deg_0_poly,
                                          deg_1_poly,
                                          deg_2_poly,
                                          deg_3_poly,
                                          deg_4_poly>;

    std::array<int, 1000 + 1> quotient_errors{};
    auto filter_quotients = [&] (int i, deg_unknown_poly q) {
      std::size_t index = static_cast<std::size_t>(i);
      if (q.index() == 2)
        possible_deg_2_quotients[index] = std::make_pair(i, std::get<2>(q));
      else if (q.index() == 3)
        possible_deg_3_quotients[index] = std::make_pair(i, std::get<3>(q));
      else if (q.index() == 4)
        possible_deg_4_quotients[index] = std::make_pair(i, std::get<4>(q));
      else
        quotient_errors[index] = i;
    };

    std::array<int, 1000 + 1> remainder_errors{};
    auto filter_remainders = [&] (int i, deg_unknown_poly r) {
      std::size_t index = static_cast<std::size_t>(i);
      if (r.index() == 0)
        possible_deg_0_remainders[index] = std::make_pair(i, std::get<0>(r));
      else if (r.index() == 1)
        possible_deg_1_remainders[index] = std::make_pair(i, std::get<1>(r));
      else if (r.index() == 2)
        possible_deg_2_remainders[index] = std::make_pair(i, std::get<2>(r));
      else
        remainder_errors[index] = i;
    };

    auto monster = [&] <int... ints> (std::integer_sequence<int, ints...>) {
      ([&]
        {
          auto a_divided_by_b = division_test<ints>();
          auto ith_quotient = a_divided_by_b.first;
          filter_quotients(ints, ith_quotient);
          auto ith_remainder = a_divided_by_b.second;
          filter_remainders(ints, ith_remainder);
        }(),
      ...);
    };
    // clang is being weird for some reason and will not work for N > ~400???
    monster(std::make_integer_sequence<int, 1000>());

    return std::make_pair(std::make_tuple(possible_deg_2_quotients,
                                          possible_deg_3_quotients,
                                          possible_deg_4_quotients,
                                          quotient_errors),
                          std::make_tuple(possible_deg_0_remainders,
                                          possible_deg_1_remainders,
                                          possible_deg_2_remainders,
                                          remainder_errors)
    );
  }();
}

int main()
{
  const std::string quotient_deg_2_file_string = filename + "2"
                                                + quotient_file_ext;
  const std::string quotient_deg_3_file_string = filename + "3"
                                                + quotient_file_ext;
  const std::string quotient_deg_4_file_string = filename + "4"
                                                + quotient_file_ext;
  const std::string remainder_deg_0_file_string = filename + "0"
                                                  + remainder_file_ext;
  const std::string remainder_deg_1_file_string = filename + "1"
                                                  + remainder_file_ext;
  const std::string remainder_deg_2_file_string = filename + "2"
                                                  + remainder_file_ext;
  /*std::ofstream quotient_deg_2_file;
  std::ofstream quotient_deg_3_file;
  std::ofstream quotient_deg_4_file;
  quotient_deg_2_file.open(quotient_deg_2_file_string);
  quotient_deg_3_file.open(quotient_deg_3_file_string);
  quotient_deg_4_file.open(quotient_deg_4_file_string);*/
  std::ofstream remainder_deg_0_file;
  //std::ofstream remainder_deg_1_file;
  //std::ofstream remainder_deg_2_file;
  remainder_deg_0_file.open(remainder_deg_0_file_string);
  //remainder_deg_1_file.open(remainder_deg_1_file_string);
  //remainder_deg_2_file.open(remainder_deg_2_file_string);
  constexpr auto quotients_and_remainders = divide_thousand_polynomials();
  //constexpr auto quotients = quotients_and_remainders.first;
  //constexpr auto deg_2_quotients = std::get<0>(quotients);
  constexpr auto remainders = quotients_and_remainders.second;
  constexpr auto deg_0_remainders = std::get<0>(remainders);
  if (remainder_deg_0_file.is_open())
  {
    std::println("succeeded to open thousand_divisions");
    for (auto pair : deg_0_remainders)
    {
      std::print(remainder_deg_0_file, "{}\t", pair.first);
      auto r = pair.second;
      for (std::size_t i = 0; i <= math::norm(r); ++i)
      {
        if (i != math::norm(r))
          std::print(remainder_deg_0_file, "{} x^{} + ", r[i], i);
        else
          std::print(remainder_deg_0_file, "{} x^{}", r[i], i);
      }
      std::print(remainder_deg_0_file, "\n");
    }
  } else
    std::println("failed to open thousand_divisions");
}
