# xfoil_light

This package provides the shared libraries `xfoil_light` and `xfoil_light_cs`, which are both versions of Mark Drela's [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) code, with the GUI features removed.  These shared libraries may be used to directly access the functions and variables in [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) from various programming languages.  The functions and variables found in the shared libraries `xfoil_light` and `xfoil_light_cs` are identical, except `xfoil_light_cs` provides a complex implementation of [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) which may be used to determine gradients of [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) functions via the complex step method.

The code in this package was originally forked from [pyXLIGHT](https://github.com/mdolab/pyXLIGHT) (which has since been renamed to [CMPLXFOIL](https://github.com/mdolab/CMPLXFOIL)).  

## Usage

### Julia

We recommend using the [Xfoil.jl](https://github.com/byuflowlab/Xfoil.jl) Julia package, which provides a wrapper for the shared libraries provided by this package.  See the installation and usage instructions [there](https://github.com/byuflowlab/Xfoil.jl).

### Python

Since at this point in time the code in this package is essentially the same as that found in the [CMPLXFOIL](https://github.com/mdolab/CMPLXFOIL) project, we recommend using the [CMPLXFOIL](https://github.com/mdolab/CMPLXFOIL) project instead.

### Other Languages

The functions and variables in the shared libraries provided by this package may be accessed using any language's foreign function interface.  Creating such an interface may or may not require a lot of work, depending on the programming language.

Assuming that no special boilerplate code is needed, the shared libraries in this package can be built using CMake. For a Unix-like system the following commands may be used:
```
cd [xfoil_light directory]
mkdir build && cd build
cmake ..
make
make install
```

By default these commands will build and install the libraries `xfoil_light` and `xfoil_light_cs` to `/usr/local/lib`  Note that root privileges (`sudo`) may be necessary for the installation command.

For examples of how to access these shared libraries see the [Xfoil.jl](https://github.com/byuflowlab/Xfoil.jl) and [CMPLXFOIL](https://github.com/mdolab/CMPLXFOIL) packages.

## Alternatives

Another shared library impelementation of XFOIL (with optional Python bindings) is provided by the [libxfoil](https://github.com/montagdude/libxfoil) project.
