# `lam.polynomial.nttp`
`lam.polynomial.nttp` is a `c++` module and a part of [LAM](https://www.github.com/colinrford/lam). If instead you prefer a header
version, and / or you would like to give it a try, please see
[here](https://godbolt.org/z/s58zEqKeY) for a Compiler Explorer implementation,
which uses `#include`s rather than `import`s, but it is relatively outdated at this point.

In `v0.2.260531`
- module name change `lam.polynomial_nttp` -> `lam.polynomial.nttp`,
- namespace name change `lam::polynomial_nttp` -> `lam::polynomial::nttp`.

## polynomial.nttp, a [LAM](https://www.github.com/colinrford/lam) library
`lam.polynomial.nttp` implements a model of a polynomial algebra in one
indeterminate in `c++` which, as the name suggests, achieves this via
Non-Type Template Parameters (NTTPs) (it's a wrapped `std::array`; thusly I
could only figure out the Division Algorithm / Euclidean Algorithm using NTTPs
– more on this in a moment). Thanks to NTTPs, combined with recent enhancements
for lambdas in `c++20`, a model of a polynomial ring $k[X]$ over a field $k$
(hopefully, of any characteristic, but at least a reasonable subset of
computable rationals) works in `constexpr` and `consteval` contexts. Such a
ring $k[X]$ is a Euclidean Domain, in which case we have a Euclidean Algorithm,
and indeed, in addition to addition, overloading `operator+`, subtraction
`operator-`, and multiplication `operator*`, all of this can be done, if one
*truly* insists, at compile time.

| operation      | implemented                                                                                                    | unit-tested                                                                                                                                                                                                                                                                                                                           |
|----------------|----------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| addition       | [`operator+`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp-univariate-algebra.cppm#L34)            | [`/tests/unit/algebra/binary_add_polynomials.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/binary_add_polynomials.cpp) and [`/tests/unit/algebra/binary_add_polynomial_with_constant.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/binary_add_polynomial_with_constant.cpp)                     |
| subtraction    | [`operator-`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp-univariate-algebra.cppm#L205)            | [`/tests/unit/algebra/binary_subtract_polynomials.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/binary_subtract_polynomials.cpp) and [`/tests/unit/algebra/binary_subtract_polynomial_with_constant.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/binary_subtract_polynomial_with_constant.cpp) |
| multiplication | [`operator*`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp-univariate-algebra.cppm#L417)            | [`/tests/unit/algebra/binary_multiply_polynomials.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/binary_multiply_polynomials.cpp) and [`/tests/unit/algebra/binary_multiply_polynomial_with_constant.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/binary_multiply_polynomial_with_constant.cpp) |
| division       | [`division_prototype<>`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp-univariate-algebra.cppm#L644) | [`/tests/unit/algebra/division.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/division.cpp)                                                                                                                                                                                                                         |
| derivative     | [`derivative`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp-univariate-algebra.cppm#L726)           | [`/tests/unit/algebra/derivative.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/derivative.cpp)                                                                                                                                                                                                                     |
| antiderivative | [`antiderivative`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp-univariate-algebra.cppm#L746)       | [`/tests/unit/algebra/antiderivative.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/algebra/antiderivative.cpp)                                                                                                                                                                                                             |
| root-finding   | [`roots<>`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp-univariate-roots.cppm#L762) | [`/tests/unit/roots/roots.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/unit/roots/roots.cpp) |

Syntax such as `3 * p + q`, `p - 1.5 * q`, `6 * p * q`, is all possible.
Dividing a polynomial `p` of degree `M` by a polynomial `q` of degree `N`, both
with coefficients represented as `double`s, is achieved by
`division_prototype<double, M, p, N, q>()`.

The biggest flaw of this implementation is division - because it is not
`operator/`... This author has not (yet?) figured out how to use the choice of
container / data structure, `std::array` (where the entries are coefficients
$a_i$, $i = 0, \ldots, n$) in a Euclidean Algorithm which both i) overloads
`operator/` and ii) compiles inside `constexpr`, `consteval` contexts. This
deeply saddens the author, but at least its still possible to achieve the
second point ii).

## `division_prototype()` (find it in [`src/polynomial_nttp-univariate-algebra.cppm`, line `644`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp-univariate-algebra.cppm#L644))
So, why was this implementation ~~doomed~~ forced from the outset (i.e. the
author's choice to use `std::array`) to rely on NTTPs to achieve simple
polynomial division at compile time? Well, that's just it, apparently, since
we cannot mark function arguments `constexpr` (still true as of `c++23`),
**it is this author's current conclusion that the only way to return
appropriately-sized quotients and remainders (at compile time) is to use a
dynamic data structure (like `std::vector`) for the Euclidean Algorithm while
temporarily storing the quotient and remainder in oversized arrays, and
subsequently copying the quotient and remainder into "right-sized" arrays.**
This technique, coined ["The constexpr 2-step" by Jason Turner](https://youtu.be/_AefJX66io8),
is useful, as it is the only (first?) way this author found polynomial division
was attainable... (at compile time... (using this implementation...))

Now the astute reader may already realize that a chain reaction has occurred.
Having been thrown into the land of template metaprogramming, as
`division_prototype()` takes 5 template parameters, 4 of which are NTTPs,
performing large numbers of tests of such a function requires something like
the use of `std::integer_sequence<int, ints...>` and fold expressions.
One could use `std::index_sequence` instead, of course.

### implementation
I have partially documented my approach in the
[wiki for this repository.](https://github.com/colinrford/polynomial_nttp/wiki/Implementation-of-%60division_prototype()%60)

## building

### preliminaries
To build with `cmake` you will need
* [`lam.concepts`](https://www.github.com/colinrford/concepts)
* `cmake`
* `ninja`
* compiler with support for `import std;`

You may also consume with `conan`, be sure to set up a profile with
support for `import std`.

For the `macOS` user, please see [here](https://www.github.com/colinrford/polynomial_nttp/wiki/Troubles-with-macOS,-Apple-Clang,-and-LLVM-Clang)
for details on getting `cmake` to avoid using Apple `clang` and `ld64` (Apple's
`clang` doesn't support `c++` modules as of circa October 1, 2025;
`ld64` still linked in my experience).

### `cmake`
1. Grab a copy
   ```bash
   git clone https://github.com/colinrford/polynomial_nttp.git
   ```
2. Inside a build directory, run `cmake`
   ```bash
   cd /to/wherever/you/placed/polynomial_nttp
   cmake -S . -B build -G Ninja
   ```
   in my case I point `cmake` to a toolchain file `-DCMAKE_TOOLCHAIN_FILE=/path/to/your/cmake/toolchain.cmake`,
   so my command is
   ```bash
   cmake -S . -B build -G Ninja -DCMAKE_TOOLCHAIN_FILE=/path/to/your/cmake/toolchain.cmake
   ```

   If this does not execute properly, and you are on `macOS`,
   please take a look at
   [this wiki entry.](https://github.com/colinrford/polynomial_nttp/wiki/Troubles-with-macOS,-Apple-Clang,-and-LLVM-Clang).
3. Run `cmake --build build` or `ninja` to build the project, which will also build tests and examples.
   ```bash
   cd build
   ninja
   ```
   If you made it past step 2. successfully, step 3 should hopefully work out
   okay. If not, and you are on `macOS`, again, take a look
   [here.](https://github.com/colinrford/polynomial_nttp/wiki/Troubles-with-macOS,-Apple-Clang,-and-LLVM-Clang)
4. (optional) in the build folder, run `ctest` to run all the unit tests
   ```bash
   ctest
   ```
### `conan`
1. Consume [`lam.concepts`](https://www.github.com/colinrford/concepts) with `conan`, i.e. ensure `lam.concepts` is in `conan` cache.
2. Consume `lam.polynomial.nttp`
   ```bash
   conan create . --profile <my-import-std-profile>
   ```

## thank you for checking out `lam.polynomial.nttp`
I do not claim this is the first of its kind, but the work here is my own.
There are possible blind spots of the author. In its current state, please
expect it to work, albeit with an occasional quirk. If you find it does not
work for some reason, please alert me somehow, either via GitHub or via
email. Any feedback is greatly appreciated.
