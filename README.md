# `polynomial_nttp`
`polynomial_nttp` is now a `c++` module. If instead you prefer a header 
version, and / or you would like to give it a try, please see 
[here](https://godbolt.org/z/s58zEqKeY) for a Compiler Explorer implementation,
which uses `#include`s rather than `import`s, but also only builds 
`test_polynomial_nttp.cpp` at this time and not `thousand_divisions.cpp`. More 
build information is found below.

## Overview
`polynomial_nttp` implements a model of a polynomial algebra in `c++`
which, as the name suggests, achieves this via  Non-Type Template Parameters 
(NTTPs) (it's a wrapped `std::array`; thusly I could only figure out the 
Division Algorithm / Euclidean Algorithm using NTTPs – more on this in a 
moment). Thanks to NTTPs, combined with recent enhancements for lambdas in 
`c++20`, a model of a polynomial ring $k[X]$ over a field $k$ (hopefully, of 
any characteristic, but at least a reasonable subset of computable rationals) 
works in `constexpr` and `consteval` contexts. Such a ring $k[X]$ is a 
Euclidean Domain, in which case we have a Euclidean Algorithm, and indeed, in 
addition to addition, overloading `operator+`, subtraction `operator-`, and 
multiplication `operator*`, all of this can be done, if one *truly* insists, 
at compile time.

I do not claim this is the first of its kind, but the work here is my own. 
There is also quite a bit of testing and refinements left to do, and possible 
blind spots of the author. In its current state, please expect it to work, 
albeit with an occasional quirk (early days). If you find it does not work for 
some reason, please alert me somehow, either via GitHub or via email. Any 
feedback is greatly appreciated.

The biggest flaw of this implementation is division (because it's not 
`operator/`). The author is not an expert in `c++` or the kinds of things this 
implentation found thrust upon itself at its onset – all this is to say, this 
author has not (yet?) figured out how to use the choice of container / data 
structure, `std::array` (where the entries are coefficients $a_i$, $i = 0, 
\ldots, n$) in a Euclidean Algorithm which both i) overloads `operator/` and 
ii) compiles inside `constexpr`, `consteval` contexts. This deeply saddens the 
author, but at least its still possible to achieve the second point ii). 

## `division_prototype()` (find it in [`src/polynomial_nttp.cpp`, line `176`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp.cpp#L176))
So, why was this implementation ~~doomed~~ forced from the outset (i.e. the 
author's choice to use `std::array`) to rely on NTTPs to achieve simple 
polynomial division at compile time? Well, that's just it, apparently, since 
we cannot mark function arguments `constexpr` (still true as of `c++23`), 
**it is this author's current conclusion that the only way to return 
appropriately-sized quotients and remainders (at compile time) is to use a 
dynamic data structure (like `std::vector`) for the Euclidean Algorithm while 
temporarily storing the quotient and remainder in oversized arrays, and
subsequently copying the quotient and remainder into "right-sized" arrays.** 
This technique, coined ["The constexpr 2-step" by Jason Turner](https://youtu.be/_AefJX66io8?si=6oBkCYUXy5VfIaOQ),
is useful, as it is the only (first?) way this author found polynomial division
was attainable... (at compile time... (using this implementation...))

Now the astute reader may already realize that a chain reaction has occured. 
Having been thrown into the land of template metaprogramming, as 
`division_prototype()` takes 5 template parameters, 4 of which are NTTPs, 
performing large numbers of tests of such a function requires (as far as I know)
something like the use of `std::integer_sequence<int, ints...>` and fold 
expressions. One could use `std::index_sequence` instead, of course.

### implementation
I documented my approach in the 
[wiki for this repository.](https://github.com/colinrford/polynomial_nttp/wiki/Implementation-of-%60division_prototype()%60)

### testing
Unit tests for the contents in `polynomial_nttp.cppm` are found in in 
`/tests/unit`. Currently they will build with a proper `Cmake` command, and 
run with `ctest`. The unit tests all have `if constexpr` statements so that
even when asserts are turned off the program will `return 0` for passing and 
`return 1` otherwise, and of course this switch is determined at compile time. 
There is a chance I incorporate a different or additional test mechanism 
(Catch2? gtest?) in the future. 

this readme section is wip, for now check out [`tests/thousand_divisions.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/thousand_divisions.cpp)
take note that this may require some compiler flags to increase the number of 
iterations, and `clang` may not want to do 1000 at a time (in my sad 
experience). currently `thousand_divisions.cpp` is a little messy, the 
chosen filenames are kind of lame, and there is some redundancy that I plan to 
reduce into a function soon enough. 

## building

### some preliminaries
Despite limited support, still, as of (circa) October 1, 2025, for `c++` 
modules, the project will build as a `c++` module. The primary build route 
currently provided builds the project using `cmake`, so you will need `cmake` 
at least version `3.30` (use `3.31.6` if you can, otherwise you'll need to 
adjust `CMakeLists.txt`), as well as at least `llvm` version `18.1.1` (newer is 
better, given status of module `std`) OR `MSVC` version `14.35`. Minimum 
requirements can be found 
[on Kitware's announcement for `import std`](https://www.kitware.com/import-std-in-cmake-3-30/).
Sadly `cmake` does not yet support `g++` or `stdlibc++` with `import std`, so I 
have included a separate `Makefile` intended exclusively for compilation with 
`g++` and `stdlibc++`; `cmake` also does not support `Makefile`s when using 
this feature, and instead solely relies on / supports `ninja` at this time 
(circa Oct. 1, 2025). In other words, to build with `cmake` you will need
* `cmake`
* `ninja`
* `clang` or `msvc`

and to build with the `Makefile`, have at least `g++` version `15.1.0` and be 
willing to make changes based upon your machine or distribution. I will try to 
elaborate on this in the near future.

For the `macOS` user, please see [here](https://github.com/colinrford/polynomial_nttp/wiki/Troubles-with-macOS,-Apple-Clang,-and-LLVM-Clang)
for details on getting `cmake` to avoid using Apple `clang` and `ld64` (Apple's 
`clang` doesn't support `c++` modules as of circa October 1, 2025; and `ld64` 
will do more than it needs to).

### building (for real this time)
1. Grab a copy
   ```bash
   git clone https://github.com/colinrford/polynomial_nttp.git
   ```
2. Inside a build directory, run `cmake`
   ```bash
   cd /to/wherever/you/placed/polynomial_nttp
   mkdir build && cd build
   cmake ..
   ```
   where it is important to note that the `cmake` command likely won't be this 
   pretty, and you will likely need to adjust the command by adding `-G Ninja` 
   and/or `-DCMAKE-TOOLCHAIN-FILE=/path/to/your/cmake/toolchain.cmake` as I 
   have had to do:
   ```bash
   cmake -DCMAKE-TOOLCHAIN-FILE=/path/to/your/cmake/toolchain.cmake .. -G Ninja
   ```
   You can alternatively add Cmake presets to shorten this command. If this\
   does not execute properly, and you are on `macOS`, please take a look at
   [this wiki entry.](https://github.com/colinrford/polynomial_nttp/wiki/Troubles-with-macOS,-Apple-Clang,-and-LLVM-Clang)
3. Run `ninja` to build the project, which will also build tests and examples 
   at this time. 
   ```bash 
   ninja
   ```
   If you made it past step 2. successfully, step 3 should hopefully work out 
   okay. If not, and you are on `macOS`, again, take a look 
   [here.](https://github.com/colinrford/polynomial_nttp/wiki/Troubles-with-macOS,-Apple-Clang,-and-LLVM-Clang)
4. (optional) run `ctest` to run all the unit tests
   ```bash
   ctest
   ```
