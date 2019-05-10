// Various helper stuff

# if !defined Helpers_h
# define Helpers_h

#ifndef NO_SIMD

#if defined(__SSE2__) || defined(__ARM_NEON)
#define HAS_SIMD 1
#endif

#ifdef __AVX512F__
#define ALIGNED __attribute__((aligned(64)))
#elif defined (__SSE2__) || defined(__ARM_NEON)
#define ALIGNED __attribute__((aligned(32)))
#else
#define ALIGNED
#endif

#else // !NO_SIMD
#define ALIGNED
#endif // !NO_SIMD

# endif
