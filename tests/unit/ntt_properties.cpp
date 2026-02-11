/*
 *  ntt_properties.cpp
 *  Verifies fundamental mathematical properties of the Number Theoretic Transform.
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::cbn::literals;

// Solinas-like Prime: 4179340454199820289 = 29 * 2^57 + 1
constexpr auto ntt_prime = 4179340454199820289_Z; 
using Field = decltype(lam::cbn::Zq(ntt_prime));

// Must specialize finite_field_traits so NTT engine recognizes ZqElement
namespace lam::polynomial::univariate
{
  template<typename T, T... Modulus>
  struct finite_field_traits<lam::cbn::ZqElement<T, Modulus...>>
  {
    static constexpr bool is_finite_field = true;
    static constexpr T modulus = []() { 
        constexpr T mods[] = {Modulus...};
        using K = lam::cbn::ZqElement<T, Modulus...>;
        return mods[0]; 
    }();
    using K = lam::cbn::ZqElement<T, Modulus...>;
    
    static constexpr K mul(const K& a, const K& b) { return a * b; }
    static constexpr K add(const K& a, const K& b) { return a + b; }
    static constexpr K sub(const K& a, const K& b) { return a - b; }
  };
}

using namespace lam::polynomial::univariate;
using namespace lam::polynomial::univariate::ntt;

void test_round_trip() {
    std::println("  Verifying Round-Trip Property: inv(ntt(x)) == x");
    
    std::vector<Field> data = { Field(1), Field(2), Field(3), Field(4), Field(0), Field(0), Field(0), Field(0) };
    auto original = data;

    // Forward
    ntt_transform(data, false);
    
    // Inverse
    ntt_transform(data, true);

    for(std::size_t i=0; i<data.size(); ++i) {
        if (data[i] != original[i]) {
            std::println("Round trip failed at index {}! Expected {}, got {}", i, original[i], data[i]);
            std::exit(1);
        }
    }
    std::println("    [Pass]");
}

void test_linearity() {
    std::println("  Verifying Linearity: ntt(a + b) == ntt(a) + ntt(b)");
    
    std::vector<Field> a = { Field(10), Field(20), Field(30), Field(40) };
    std::vector<Field> b = { Field(5), Field(5), Field(5), Field(5) };
    std::vector<Field> sum_ab(4);
    
    for(std::size_t i=0; i<4; ++i) sum_ab[i] = a[i] + b[i];

    ntt_transform(a, false);
    ntt_transform(b, false);
    ntt_transform(sum_ab, false); // NTT(a+b)

    for(std::size_t i=0; i<4; ++i) {
        if (sum_ab[i] != (a[i] + b[i])) { // Should equal NTT(a) + NTT(b)
            std::println("Linearity failed at index {}", i);
            std::exit(1);
        }
    }
    std::println("    [Pass]");
}

// Convolution Theorem: NTT(a * b) = NTT(a) . NTT(b)
// Where '*' is polynomial convolution (cyclic), and '.' is pointwise
void test_convolution_theorem() {
    std::println("  Verifying Convolution Theorem");
    
    // N=4.
    // a = 1 + x
    // b = 1 - x
    // a*b (linear) = 1 - x^2. 
    // Cyclic convolution mod 4: 1 - x^2. (Coeffs: 1, 0, -1, 0)
    
    std::vector<Field> a = { Field(1), Field(1), Field(0), Field(0) };
    std::vector<Field> b = { Field(1), Field(4179340454199820288ULL), Field(0), Field(0) }; // 1, -1, 0, 0
    
    // Expected Result: 1, 0, -1, 0
    // 4179340454199820289 - 1 = 4179340454199820288
    std::vector<Field> expected = { Field(1), Field(0), Field(4179340454199820288ULL), Field(0) };

    ntt_transform(a, false);
    ntt_transform(b, false);
    
    // Pointwise Multiply
    std::vector<Field> c(4);
    for(std::size_t i=0; i<4; ++i) c[i] = a[i] * b[i];
    
    // Inverse Transform
    ntt_transform(c, true);
    
    for(std::size_t i=0; i<4; ++i) {
        if (c[i] != expected[i]) {
            std::println("Convolution failed at index {}! Expected {}, got {}", i, expected[i], c[i]);
            std::exit(1);
        }
    }
    std::println("    [Pass]");
}

int main() {
    std::println("Running NTT Property Tests...");
    test_round_trip();
    test_linearity();
    test_convolution_theorem();
    std::println("All NTT Property Tests Passed.");
    return 0;
}
