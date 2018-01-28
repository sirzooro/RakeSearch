// Various helper stuff

# if !defined Helpers_h
# define Helpers_h

#ifdef __AVX512F__
#define ALIGNED __attribute__((aligned(64)))
#elif defined (__SSE2__) || defined(__ARM_NEON)
#define ALIGNED __attribute__((aligned(32)))
#else
#define ALIGNED
#endif

# endif
