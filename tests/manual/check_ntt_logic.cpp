
import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate::ntt;

void check_bit_reversal() {
    std::println("Checking Bit Reversal for N=16...");
    std::size_t n = 16;
    std::vector<int> data(n);
    for(std::size_t i=0; i<n; ++i) data[i] = i;

    // Direct copy of bit-reversal logic from ntt_transform
    std::size_t j = 0;
    for (std::size_t i = 1; i < n; i++)
    {
      std::size_t bit = n >> 1;
      while (j & bit)
      {
        j ^= bit;
        bit >>= 1;
      }
      j ^= bit;
      if (i < j)
      {
        std::swap(data[i], data[j]);
      }
    }

    // Expected: 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    std::vector<int> expected = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
    bool pass = true;
    for(std::size_t i=0; i<n; ++i) {
        if (data[i] != expected[i]) {
            std::println("Mismatch at {}: Got {}, Expected {}", i, data[i], expected[i]);
            pass = false;
        }
    }
    std::println("Bit Reversal N=16: {}", pass ? "PASS" : "FAIL");
}

void check_root_properties() {
    std::println("Checking Root Properties for N=4096...");
    std::uint64_t P = 18446744069414584321ULL;
    std::uint64_t N = 4096;
    
    // Call find_nth_root_of_unity
    std::uint64_t root = find_nth_root_of_unity<std::uint64_t>(N, P);
    std::println("Root found: {}", root);
    
    // Check order
    auto w_n = power(root, N, P);
    std::println("root^N (mod P) = {} (Expected 1)", w_n);
    
    auto w_half = power(root, N/2, P);
    std::println("root^(N/2) (mod P) = {} (Expected P-1)", w_half); // Should be -1
    
    if (w_n == 1 && w_half != 1) {
         std::println("Root Order: PASS");
    } else {
         std::println("Root Order: FAIL");
    }
}

void check_modular_inverse() {
    std::println("Checking Modular Inverse for N=4096...");
    std::uint64_t P = 18446744069414584321ULL;
    std::uint64_t N = 4096;

    std::uint64_t inv_N = mod_inverse(N, P);
    unsigned __int128 prod = static_cast<unsigned __int128>(inv_N) * N;
    std::uint64_t check = static_cast<std::uint64_t>(prod % P);
    
    std::println("inv(N) = {}", inv_N);
    std::println("inv(N) * N (mod P) = {} (Expected 1)", check);
    
    if (check == 1) {
        std::println("Modular Inverse: PASS");
    } else {
        std::println("Modular Inverse: FAIL");
    }
}

import lam.ctbignum;

constexpr std::uint64_t P = 18446744069414584321ULL;
using Zq = lam::cbn::ZqElement<std::uint64_t, P>;

// Replicate benchmark traits EXACTLY
namespace lam::polynomial::univariate {
    template<>
    struct finite_field_traits<Zq> {
        using K = Zq;
        static constexpr bool is_finite_field = true;
        static constexpr std::uint64_t modulus = P;

        static constexpr K mul(const K& a, const K& b) {
            unsigned __int128 prod = static_cast<unsigned __int128>(a.data[0]) * b.data[0];
            std::uint64_t solinas_res = lam::polynomial::univariate::ntt::reduce_solinas(prod);
            
            // Cross-check
            if (!std::is_constant_evaluated()) {
                 unsigned __int128 safe_res = prod % P;
                 if (solinas_res != static_cast<std::uint64_t>(safe_res)) {
                      std::println("MUL MISMATCH! a={:x} b={:x} prod={:x}", a.data[0], b.data[0], (std::uint64_t)prod);
                      std::println("  Solinas: {:x}", solinas_res);
                      std::println("  Safe:    {:x}", (std::uint64_t)safe_res);
                      std::exit(1); 
                 }
            }
            return K(solinas_res);
        }

        static constexpr K add(const K& a, const K& b) { 
            K res = a + b;
            if (!std::is_constant_evaluated()) {
                // Trap the specific failure case: 1 + (P-1)
                if (a.data[0] == 1 && b.data[0] == P-1) {
                    std::println("TRAIT ADD(1, P-1) -> {:x}", res.data[0]);
                }
            }
            return res;
        }
        static constexpr K sub(const K& a, const K& b) { 
            K res = a - b;
            if (!std::is_constant_evaluated()) {
                // Trap the specific failure case: 1 - (P-1)
                if (a.data[0] == 1 && b.data[0] == P-1) {
                    std::println("TRAIT SUB(1, P-1) -> {:x}", res.data[0]);
                }
            }
            return res; 
        }
    };
}

using Wrapper = Zq;

// Naive DFT for verification
std::vector<std::uint64_t> naive_dft(const std::vector<std::uint64_t>& input, std::uint64_t root) {
    std::size_t N = input.size();
    std::vector<std::uint64_t> output(N);
    for (std::size_t k = 0; k < N; ++k) {
        unsigned __int128 sum = 0;
        unsigned __int128 w_k = 1; // root^0
        unsigned __int128 w_step = power(root, k, P);
        
        for (std::size_t n = 0; n < N; ++n) {
            unsigned __int128 term = mul_mod(input[n], static_cast<std::uint64_t>(w_k), P);
            sum += term;
            if (sum >= P) sum -= P;
            w_k = mul_mod(static_cast<std::uint64_t>(w_k), static_cast<std::uint64_t>(w_step), P);
        }
        output[k] = static_cast<std::uint64_t>(sum % P);
    }
    return output;
}

void check_ntt_transform_logic() {
    std::println("Checking Full NTT Logic (Local vs Naive DFT)...");
    
    std::size_t N = 1024; // Test blocking size
    std::uint64_t root = find_nth_root_of_unity<std::uint64_t>(N, P);
    
    std::vector<Wrapper> ntt_data(N);
    std::vector<std::uint64_t> raw_data(N);
    
    for(std::size_t i=0; i<N; ++i) {
        ntt_data[i] = Wrapper(i + 1);
        raw_data[i] = i + 1;
    }
    
    // Run NTT
    ntt_transform(ntt_data, false);
    
    // Run Naive
    auto expected = naive_dft(raw_data, root);
    
    bool pass = true;
    for(std::size_t i=0; i<N; ++i) {
        if (ntt_data[i].data[0] != expected[i]) {
            std::println("Mismatch at {}: Naive={} NTT={}", i, expected[i], ntt_data[i].data[0]);
            pass = false;
            break; // Stop spam
        }
    }
    
    std::println("NTT Transform N={}: {}", N, pass ? "PASS" : "FAIL");
}

void check_convolution_logic() {
    std::println("Checking Convolution Logic (Forward -> Pointwise -> Inverse)...");
    
    std::size_t N = 4096; // Pad size for degree 1024
    
    // Create A and B
    std::vector<Wrapper> a(N, Wrapper(0));
    std::vector<Wrapper> b(N, Wrapper(0));
    
    // A = 1 => {1, 0, ...}
    // B = 1 => {1, 0, ...}
    // A*B = 1 => {1, 0, ...}
    a[0] = Wrapper(1);
    b[0] = Wrapper(1);
    
    // Forward
    ntt_transform(a, false);
    ntt_transform(b, false);
    
    // Pointwise
    std::vector<Wrapper> c(N);
    for(std::size_t i=0; i<N; ++i) {
        c[i] = lam::polynomial::univariate::finite_field_traits<Wrapper>::mul(a[i], b[i]);
    }
    
    // Inverse
    ntt_transform(c, true);
    
    // Check
    bool pass = true;
    if (c[0].data[0] != 1) pass = false;
    for(std::size_t i=1; i<N; ++i) {
        if (c[i].data[0] != 0) pass = false;
    }
    
    if (c[0].data[0] != 1) 
        std::println("Mismatch at 0: Expected 1, Got {}", c[0].data[0]);

    std::println("Convolution Identity (1*1): {}", pass ? "PASS" : "FAIL");
    
    // Check non-trivial
    // A = {1, 1, 0...} (1+x)
    // B = {1, -1, 0...} (1-x) => -1 is P-1
    // C = {1, 0, -1, 0...} (1-x^2)
    
    std::fill(a.begin(), a.end(), Wrapper(0));
    std::fill(b.begin(), b.end(), Wrapper(0));
    
    a[0] = 1; a[1] = 1;
    b[0] = 1; b[1] = P - 1;
    
    ntt_transform(a, false);
    ntt_transform(b, false);
    
    for(std::size_t i=0; i<N; ++i) {
         c[i] = lam::polynomial::univariate::finite_field_traits<Wrapper>::mul(a[i], b[i]);
    }
    
    ntt_transform(c, true);
    
    Wrapper exp0(1);
    Wrapper exp1(0);
    Wrapper exp2(P-1);
    
    if (c[0].data[0] == exp0.data[0] && c[1].data[0] == exp1.data[0] && c[2].data[0] == exp2.data[0]) {
        std::println("Convolution Non-Trivial: PASS");
    } else {
        std::println("Convolution Non-Trivial: FAIL");
        std::println("0: {} (Exp 1)", c[0].data[0]);
        std::println("1: {} (Exp 0)", c[1].data[0]);
        std::println("2: {} (Exp P-1)", c[2].data[0]);
    }
}

void check_field_arithmetic() {
    std::println("Checking ZqElement Arithmetic...");
    
    // ZqElement is aliased as Wrapper aka Zq
    Zq a(P - 1);
    Zq b(1);
    std::println("a: {:x}", a.data[0]);
    std::println("b: {:x}", b.data[0]);
    
    Zq c = a + b; // Should be 0
    std::println("c: {:x}", c.data[0]);
    
    // Manual 128-bit check
    unsigned __int128 sum = (unsigned __int128)a.data[0] + b.data[0];
    std::uint64_t res = (sum >= P) ? (std::uint64_t)(sum - P) : (std::uint64_t)sum;
    std::println("Manual 128-bit sum: {:x}", res);

    if (c.data[0] == 0) std::println("Add Wrap (P-1 + 1): PASS");
    else std::println("Add Wrap (P-1 + 1): FAIL. Got {:x} (Expected 0)", c.data[0]);

    Zq d(0);
    Zq e(1);
    Zq f = d - e; // Should be P-1
    std::println("f: {:x}", f.data[0]);
    
    if (f.data[0] == P - 1) std::println("Sub Wrap (0 - 1): PASS");
    else std::println("Sub Wrap (0 - 1): FAIL. Got {:x}", f.data[0]);
    
    // Check associativity
    // (P-5) + 10 -> 5
    Zq g(P - 5);
    Zq h(10);
    Zq i = g + h;
    std::println("i: {:x}", i.data[0]);
    if (i.data[0] == 5) std::println("Add Assoc: PASS");
    else std::println("Add Assoc: FAIL. Got {:x}", i.data[0]);

    // Commutativity: 1 + (P-1)
    Zq l(1);
    Zq m(P - 1);
    Zq n = l + m;
    std::println("1 + (P-1): {:x} (Expected 0)", n.data[0]);
    if (n.data[0] == 0) std::println("Add Commute: PASS");
    else std::println("Add Commute: FAIL");

    // Subtraction: 1 - (P-1) = 2
    Zq o = l - m;
    std::println("1 - (P-1): {:x} (Expected 2)", o.data[0]);
    if (o.data[0] == 2) std::println("Sub P-1: PASS");
    else std::println("Sub P-1: FAIL");

    // Check (P-1) + (P-1) = P-2
    Zq j(P - 1);
    Zq k = j + j;
    std::println("k: {:x} (Expected {:x})", k.data[0], P-2);
    std::println("k: {:x} (Expected {:x})", k.data[0], P-2);
    if (k.data[0] == P - 2) std::println("Add Large: PASS");
    else std::println("Add Large: FAIL. Got {:x}", k.data[0]);
}

void check_multiplication() {
    std::println("Checking ZqElement Multiplication...");
    
    // (P-1) * (P-1) = (-1)*(-1) = 1
    Zq a(P - 1);
    Zq b = a * a;
    std::println("(P-1)*(P-1): {:x} (Expected 1)", b.data[0]);
    
    if (b.data[0] == 1) std::println("Mul NegOne: PASS");
    else std::println("Mul NegOne: FAIL");

    // (P-1) * 1 = P-1
    Zq c = a * Zq(1);
    if (c.data[0] == P-1) std::println("Mul Identity: PASS");
    else std::println("Mul Identity: FAIL. Got {:x}", c.data[0]);
    
    // 2^32 * 2^32 = 2^64
    // 2^64 mod P = 2^32 - 1
    std::uint64_t val = 1ULL << 32;
    Zq d(val);
    Zq e = d * d;
    std::println("2^32 * 2^32: {:x} (Expected {:x})", e.data[0], val - 1);
    
    if (e.data[0] == val - 1) std::println("Mul 2^64: PASS");
    else std::println("Mul 2^64: FAIL");
}

void check_ntt_edge_cases() {
    std::println("Checking NTT Edge Cases (Input P-1)...");
    
    // Test NTT on [1, P-1] (conceptually 1 + (P-1)x)
    // N=2
    std::size_t N = 2;
    std::vector<Zq> data(N);
    data[0] = 1;
    data[1] = P - 1;
    
    // Expected Forward NTT:
    // val[0] = 1 + (P-1) = P = 0
    // val[1] = 1 - (P-1) = 1 - (-1) = 2
    
    ntt_transform(data, false);
    
    std::println("NTT([1, P-1])[0]: {:x} (Expected 0)", data[0].data[0]);
    std::println("NTT([1, P-1])[1]: {:x} (Expected 2)", data[1].data[0]);
    
    if (data[0].data[0] == 0 && data[1].data[0] == 2) {
        std::println("NTT P-1: PASS");
    } else {
        std::println("NTT P-1: FAIL");
    }
}

int main() {
    check_bit_reversal();
    check_root_properties();
    check_modular_inverse();
    check_field_arithmetic();
    check_multiplication();
    check_ntt_edge_cases();
    check_ntt_transform_logic();
    check_convolution_logic();
    return 0;
}
