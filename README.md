# RakeSearch
Rake search of Diagonal Latin Squares

## Compilation

Application is compiled using gcc and make. Windows version is compiled using MinGW gcc crosscompiler.

To compile everything, you must have BOINC client libraries. Make sure you compile them using the same gcc version as this app, otherwise you may get link errors.

When compiling app for x86/x86_64, you will need gcc 8.x (I used gcc 8.2). Older gcc versions may not support some options used in Makefile.

To compile app, enter RakeSearch/RakeDiagSearch/RakeDiagSearch directory first. Open Makefile and update `BOINC_DIR` variable, so it will point to place where you have BOINC client library and its include files. After doing so, type `make` to start compilation.

Makefile supports number of extra parameters. Here are ones used for x86 and x86_64:

- `SSE2=1` - enable SSE2 instructions (x86 and x86_64)
- `SSSE3=1` - enable SSSE3 instructions (x86 and x86_64)
- `AVX=1` - enable AVX instructions (x86_64 only)
- `AVX2=1` - enable AVX2 and BMI1/2 instructions (x86_64 only)
- `AVX512=1` - enable AVX512 instructions (x86_64 only)

Note: SSE2 is always enabled on x86_64, support for it is part of AMD64 specification.

When compiling for ARM, you can enable NEON instruction by using `NEON=1` parameter. AARCH64 supports it by default, so no need to explicitly enable it.

You can also specify target platform:
- `M32=1` - compile 32-bit app version (used for Linux app)
- `MinGW64=1` - compile 64-bit app for Windows using MinGW crosscompiler
- `MinGW32=1` - compile 32-bit app for Windows using MinGW crosscompiler
- `AARCH64=1` - compile 64-bit app for Linux using AARCH64 crosscompiler
- `ARM=1` - compile 32-bit app for Linux using ARM crosscompiler

When compiling for ARM or AARCH64 using native compiler, no extra options are needed. You can only add `NEON=1` to enable NEON for ARM.

It is also possible to disable SIMD instructions with `NOSIMD=1` parameter. Thismay be handy if you want to disable SIMD instructions on platform which always have them enabled - namely SSE2 on x86_s64, NEON on ARM64.
