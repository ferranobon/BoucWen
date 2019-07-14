# Bouc-Wen implementation

This is a C++ Bouc-Wen implementation supporting the following cases:

- Bouc-Wen model
- Bouc-Wen model with material degradation
- Bouc-Wen with pitching and degradation (Baber-Noori)

It reads displacements stored in a txt file, and outputs input displacement and force.

## Compilation

The following options are available with CMAKE

### Build type

Can be either Debug or Release.

### Compiler

Supported compilers are: g++, clang, Intel C++ and PGI C++ compilers

### Compiling

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ ..
```
