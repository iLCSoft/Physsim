# Physsim
[![Build Status](https://travis-ci.org/iLCSoft/Physsim.svg?branch=master)](https://travis-ci.org/iLCSoft/Physsim)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/12399/badge.svg)](https://scan.coverity.com/projects/ilcsoft-physsim)

Physsim is a matrix element package

Physsim is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Building Physsim and the example processors

In order to build Physsim the easiest way to satisfy all dependencies is to either use a Key4hep or a recent iLCSoft release. The instructions are very similar for both, the major difference is in the invokation of the `cmake` command, where `-C ${ILCSOFT}/ILCSoft.cmake` has to be replaced by `-DCMAKE_CXX_STANDARD=17` when using a Key4hep release.

- Setup a recent iLCSoft release, e.g.
```bash
source /cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-01/init_ilcsoft.sh
```

- Clone the Physsim repository
```bash
git clone https://github.com/iLCSoft/Physsim
cd Physsim
```

- Build and install Physsim
```bash
mkdir build
cd build
cmake -C ${ILCSOFT}/ILCSoft.cmake -DCMAKE_INSTALL_PREFIX=../install ..
make install  # use -j<N> to build with N proceses in parallel
```

- Setup the environment to use the Physsim that you have just built (assuming you are still in the build directory). **Note that the `lib64` folder might be called `lib` depending on the OS you are using, please change accordingly**
```bash
cd ../install
export CMAKE_PREFIX_PATH=$(pwd):$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=$(pwd)/lib64:$LD_LIBRARY_PATH
cd ..
```

- Build the example processors (assuming that you start from the `Physsim` top-level directory)
```bash
cd example_stdhep
mkdir build
cd build
cmake -C ${ILCSOFT}/ILCSoft.cmake -DCMAKE_INSTALL_PREFIX=../install ..
make install  # use -j<N> to build with N processes in parallel
```

- Afterwards, the usual steps of adding the processors to `MARLIN_DLL` are necessary in order for Marlin to be able to find them

## License and Copyright
Copyright (C) 2005-2017, Physsim Authors

Physsim is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.
