# `polynomial_nttp`
`polynomial_nttp` is a fun little implementation of polynomial algebra in `c++`
which, as the name suggests, achieves this via  Non-Type Template Parameters 
(NTTPs) (it's a wrapped `std::array`; however, I could only figure out the 
Division Algorithm / Euclidean Algorithm using NTTPs – more on this in a 
moment). Thanks to NTTPs, combined with recent enhancements for lambdas in 
`c++20`, a model of a polynomial ring $k[X]$ over a field $k$ (hopefully, of 
any characteristic, but at least a reasonable subset of computable rationals) 
works in `constexpr` and `consteval` contexts. Such a ring $k[X]$ is a 
Euclidean Domain, in which case we have a Euclidean Algorithm, and indeed, in 
addition to addition, overloading `operator+`, subtraction `operator-`, and 
multiplication `operator*`, all of this can be done, if one *truly* insists, 
at compile time.

I do not claim this is the first of its kind, it is simply something I found 
quite stimulating for a little while. There is also quite a bit of testing and 
refinements left to do, and possible blind spots of the author. In its current 
state, please expect it to work, albeit with an occasional quirk (early days). 
If you find it does not work for some reason, please alert me somehow, either 
via GitHub or via email. Any feedback is greatly appreciated.

The biggest flaw of this implementation is division (because its not 
`operator/`). The author is not an expert in `c++` or the kinds of things this 
implentation found thrust upon itself at its onset – all this is to say, this 
author has not (yet?) figured out how to use the choice of container / data 
structure, `std::array` (where the entries are coefficients $a_i$, $i = 0, 
\ldots, n$) in a Euclidean Algorithm which both i) overloads `operator/` and 
ii) compiles inside `constexpr`, `consteval` contexts. This deeply saddens the 
author, but at least its still possible to achieve the second point ii). 

## `division_prototype()` (find it in [`src/polynomial_nttp.cpp`, line `173`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp.cpp#L173))
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
performing large number of tests of such a function requires (as far as I know)
something like the use of `std::integer_sequence<int, ints...>` and fold 
expressions. One could use `std::index_sequence` instead, of course.

### implementation
this readme section is wip, for now check out
[`src/polynomial_nttp.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/src/polynomial_nttp.cpp)
### testing
this readme section is wip, for now check out [`tests/thousand_divisions.cpp`](https://github.com/colinrford/polynomial_nttp/blob/main/tests/thousand_divisions.cpp)
take note that this may require some compiler flags to increase the number of 
iterations, and `clang` may not want to do 1000 at a time (in my sad 
experience). currently `thousand_divisions.cpp` is a little messy, the 
chosen filenames are kind of lame, and there is some redundancy that I plan to 
reduce into a function soon enough. 

## building
this readme section is wip, for now e.g. in the `tests` directory
```c++
g++ -std=c++23 -I../src test_polynomial_nttp.cpp -o test_polynomial
g++ -std=c++23 -I../src thousand_divisions.cpp -o thousand_divisions
```

