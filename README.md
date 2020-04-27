# Cellogram

[![Build Status](https://travis-ci.com/cellogram/cellogram.svg?branch=master)](https://travis-ci.com/cellogram/cellogram)
[![Build Status](https://dev.azure.com/cellogram/cellogram/_apis/build/status/cellogram.cellogram?branchName=master)](https://dev.azure.com/cellogram/cellogram/_build/latest?definitionId=1&branchName=master)

Automatic registration of deformed cell arrays.

## Building Cellogram

After cloning this git repository, do the following (on macOS/Linux):

```
mkdir build
cd build
cmake ..
make -j 8
```

## License

The code of Cellogram itself is licensed under [MIT License](LICENSE). However, please be mindful of third-party libraries which are used by Cellogram, and may be available under a different license. In particular, the `src/tcdf` module follows the Academic License (see header). The files in the following folders are available under MIT:
```
app/
src/cellogram/
src/points_untangler/
```
