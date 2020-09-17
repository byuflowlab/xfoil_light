# xfoil-light

This package provides the shared libraries `xfoil_light` and `xfoil_light_cs`, which are both versions of Mark Drela's [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) code, with the GUI features removed.  These shared libraries may be used to directly access the functions and variables in XFOIL from various programming languages.  The functions and variables found in the shared libraries `xfoil_light` and `xfoil_light_cs` are identical, except `xfoil_light_cs` provides a complex implementation of XFOIL which may be used to determine gradients of XFOIL functions via the complex step method.

The code in this package was originally forked from the [pyXLIGHT](https://github.com/mdolab/pyXLIGHT) project.  

## Usage

### Julia

We recommend using the [Xfoil.jl](https://github.com/byuflowlab/Xfoil.jl) Julia package

### Python

Since at this point in time the code in this package is essentially the same as that found in th [pyXLIGHT](https://github.com/mdolab/pyXLIGHT) project we recommend using the [pyXLIGHT](https://github.com/mdolab/pyXLIGHT) project instead.

### Other Languages

The functions and variables in the shared libraries provided by this package may be accessed using any language's foreign function interface.  Creating such an interface may or may not require a lot of work, depending on the programming language.

Assuming that no special boilerplate code is needed, the shared libraries in this package can be built using CMake. For a Unix-like system the following commands may be used:
```
cd [libxfoil directory]
mkdir build && cd build
cmake ..
make
make install
```

By default these commands will build and install the libraries `libxfoil` and `libxfoil_cs` to `/usr/local/lib`  Note that root privileges (`sudo`) may be necessary for the installation command.

For examples of how to access these shared libraries see the [Xfoil.jl](https://github.com/byuflowlab/Xfoil.jl) and [pyXLIGHT](https://github.com/mdolab/pyXLIGHT) packages.

## Alternatives

Another shared library impelementation of XFOIL (with optional Python bindings) is provided by the [libxfoil](https://github.com/montagdude/libxfoil) project.
