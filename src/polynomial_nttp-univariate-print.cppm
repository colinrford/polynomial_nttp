/*
 *  polynomial_nttp-univariate-print.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

export module lam.polynomial_nttp:univariate.print;

import std;
import :univariate.structure;

template<lam::polynomial::univariate::ring_element_c_weak R, std::size_t N>
struct std::formatter<lam::polynomial_nttp<R, N>>
{
  char variable = 'x';
  int precision = -1;
  bool show_plus_sign = false;

  constexpr auto parse(std::format_parse_context& ctx)
  {
    auto it = ctx.begin();
    auto end = ctx.end();

    while (it != end && *it != '}')
    {
      if (*it == '.')
      {
        ++it;
        precision = 0;
        while (it != end && *it >= '0' && *it <= '9')
        {
          precision = precision * 10 + (*it - '0');
          ++it;
        }
      } else if (*it == '+')
      {
        show_plus_sign = true;
        ++it;
      } else if (*it >= 'a' && *it <= 'z')
      {
        variable = *it;
        ++it;
      } else
        ++it;
    }
    return it;
  }

  template<typename FormatContext>
  auto format(const lam::polynomial_nttp<R, N>& p, FormatContext& ctx) const
  {
    auto out = ctx.out();
    bool first = true;

    for (std::size_t i = 0; i <= N; ++i)
    {
      R coeff = p.coefficients[i];
      // Skip negligible coefficients (except for constant polynomial)
      if (lam::is_negligible(coeff) && N > 0)
        continue;
      // Handle sign and spacing
      if (!first)
      {
        if (coeff >= R(0))
          out = std::format_to(out, " + ");
        else
        {
          out = std::format_to(out, " - ");
          coeff = -coeff;
        }
      } else if (coeff < R(0))
      {
        out = std::format_to(out, "-");
        coeff = -coeff;
      }
      first = false;

      // Format coefficient with optional precision
      auto format_coeff = [&](R c) {
        if (precision >= 0)
          return std::format("{:.{}f}", static_cast<double>(c), precision);
        else
          return std::format("{}", c);
      };

      bool coeff_is_one = lam::is_approx_equal(coeff, R(1));
      // Constant term - always show coefficient
      if (i == 0)
        out = std::format_to(out, "{}", format_coeff(coeff));
      else if (i == 1)
      {
        if (coeff_is_one) // Linear term
          out = std::format_to(out, "{}", variable);
        else
          out = std::format_to(out, "{}{}", format_coeff(coeff), variable);
      } else
      {
        if (coeff_is_one) // Higher degree terms
          out = std::format_to(out, "{}^{}", variable, i);
        else
          out = std::format_to(out, "{}{}^{}", format_coeff(coeff), variable, i);
      }
    }
    // Handle zero polynomial
    if (first)
    {
      if (precision >= 0)
        out = std::format_to(out, "{:.{}f}", 0.0, precision);
      else
        out = std::format_to(out, "0");
    }
    return out;
  }
};
