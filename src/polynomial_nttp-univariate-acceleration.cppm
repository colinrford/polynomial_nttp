/*
 *  polynomial_nttp-univariate-acceleration.cppm
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */

module;

#ifdef LAM_USE_ACCELERATE
// Avoid including <Accelerate/Accelerate.h> in module interfaces to prevent macro pollution
extern "C"
{
  typedef unsigned long vDSP_Length;
  typedef long vDSP_Stride;

  // Split complex for double precision
  struct DSP_double_split_complex
  {
    double* realp;
    double* imagp;
  };

  // Split complex for single precision
  struct DSP_float_split_complex
  {
    float* realp;
    float* imagp;
  };

  // Double precision
  void vDSP_vsmaD(const double* __A, vDSP_Stride __IA, const double* __S, const double* __C, vDSP_Stride __IC,
                  double* __D, vDSP_Stride __ID, vDSP_Length __N);

  void vDSP_zvmaD(const DSP_double_split_complex* __A, vDSP_Stride __IA, const DSP_double_split_complex* __B,
                  vDSP_Stride __IB, const DSP_double_split_complex* __C, vDSP_Stride __IC,
                  const DSP_double_split_complex* __D, vDSP_Stride __ID, vDSP_Length __N);

  void vDSP_vaddD(const double* __A, vDSP_Stride __IA, const double* __B, vDSP_Stride __IB, double* __C,
                  vDSP_Stride __IC, vDSP_Length __N);

  void vDSP_vsubD(const double* __A, vDSP_Stride __IA, const double* __B, vDSP_Stride __IB, double* __C,
                  vDSP_Stride __IC, vDSP_Length __N);

  // Single precision
  void vDSP_vsma(const float* __A, vDSP_Stride __IA, const float* __S, const float* __C, vDSP_Stride __IC, float* __D,
                 vDSP_Stride __ID, vDSP_Length __N);

  void vDSP_zvma(const DSP_float_split_complex* __A, vDSP_Stride __IA, const DSP_float_split_complex* __B,
                 vDSP_Stride __IB, const DSP_float_split_complex* __C, vDSP_Stride __IC,
                 const DSP_float_split_complex* __D, vDSP_Stride __ID, vDSP_Length __N);

  void vDSP_vadd(const float* __A, vDSP_Stride __IA, const float* __B, vDSP_Stride __IB, float* __C, vDSP_Stride __IC,
                 vDSP_Length __N);

  void vDSP_vsub(const float* __A, vDSP_Stride __IA, const float* __B, vDSP_Stride __IB, float* __C, vDSP_Stride __IC,
                 vDSP_Length __N);
}
#endif

#ifdef LAM_USE_BLAS
extern "C"
{
  // Double precision
  void cblas_daxpy(const int N, const double alpha, const double* X, const int incX, double* Y, const int incY);
  void cblas_zaxpy(const int N, const void* alpha, const void* X, const int incX, void* Y, const int incY);
  // Single precision
  void cblas_saxpy(const int N, const float alpha, const float* X, const int incX, float* Y, const int incY);
  void cblas_caxpy(const int N, const void* alpha, const void* X, const int incX, void* Y, const int incY);
}
#endif

export module lam.polynomial_nttp:univariate.acceleration;

export namespace lam::polynomial::univariate::acceleration
{
#ifdef LAM_USE_ACCELERATE
using ::DSP_double_split_complex;
using ::DSP_float_split_complex;
using ::vDSP_Length;
using ::vDSP_Stride;
using ::vDSP_vadd;
using ::vDSP_vaddD;
using ::vDSP_vsma;
using ::vDSP_vsmaD;
using ::vDSP_vsub;
using ::vDSP_vsubD;
using ::vDSP_zvma;
using ::vDSP_zvmaD;
#endif

#ifdef LAM_USE_BLAS
using ::cblas_caxpy;
using ::cblas_daxpy;
using ::cblas_saxpy;
using ::cblas_zaxpy;
#endif
} // namespace lam::polynomial::univariate::acceleration
