# <u>propagate</u>
1D TDSE FFT wavepacket propagation and imaginary time propagation for finding ground states.

---
### Details
**Language:** C++23

**Machine requirements:** ARM64 or Intel x86_64 (use $ cmake -B build -DBUILD_INTEL=ON). RAM usage varies with run settings.

**External dependencies:** ArmPL or imkl (specifically for fftw). 

**Compiler used and version:** gcc (Apple clang version 17.0.0 (clang-1700.0.13.5)) or g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0

---

To configure build on x86_64 systems (Intel), use

```bash
cmake -B build -DBUILD_INTEL=ON
cmake --build build
```
*NOTE: OpenMP flags may need to be changed*

Static linking may be needed to compile on some systems (like HPCs). For example,

```bash
g++-14 -static -std=c++23 propagate_td.cpp propagate.cpp filetools.cpp fftw_complex_tools.cpp console_tools.cpp -o propagate_td -lfftw3 -lm
```

For ARM64 processors,

```bash
cmake -B build
cmake --build build
```