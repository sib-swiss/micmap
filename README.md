# micmap

	Mapping of short reads in fastq format onto a reference

## License

    Copyright (C) SIB  - Swiss Institute of Bioinformatics,   2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
    Copyright (C) UNIL - University of Lausanne, Switzerland  2019-2020 Nicolas Guex and Christian Iseli
    Copyright (C) EPFL - EPFL, Lausanne, Switzerland          2020-2022 Christian Iseli


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


## Authors

	Code:       Nicolas Guex, Thierry Schuepbach and Christian Iseli
	Contacts:   Nicolas.Guex@unil.ch and Christian.Iseli@epfl.ch
	Repository: https://github.com/sib-swiss/micmap

## Info

	Machine :	Unix
	Language:	C
	Requires:	AVX512, pthread

	Version information

	Version:	2.4  November 2022 Public release of code under GPL2+ license



## Build instructions

- need to manually setup the json-c library in third-party subdirectory

```bash
mkdir /tmp/tst
cd /tmp/tst
git clone --recursive git@github.com:sib-swiss/micmap.git
cd micmap/third-party/json-c/
mkdir build
cd build/
../cmake-configure
make
mv * .. # say n to move the existing directories
cd ..
rm -rf build
cd ../..
```

- create a build directory and build the code

```bash
mkdir build
cd build
cmake ../ -DGTL_TESTS_PATH=/tmp/GTL -DCMAKE_INSTALL_PREFIX=/data6/tools -DUSE_RPATH=0 -DCMAKE_INSTALL_RPATH=\$ORIGIN/../lib
make
```

- it is also possible to compile using intel's compiler

```bash
source /opt/intel/bin/compilervars.sh intel64
mkdir build
cd build
CC=icc cmake ../ -DGTL_TESTS_PATH=/tmp/GTL -DCMAKE_INSTALL_PREFIX=/data6/tools -DUSE_RPATH=0 -DCMAKE_INSTALL_RPATH=\$ORIGIN/../lib
ccmake .. # add -static-intel to the CMAKE_C_FLAGS_* flags (need to toggle in expert mode, edit, then configure and generate)
cmake .. # to check the CFLAGS contain -static-intel
make
```

- can also build a package and/or do testing

```bash
make package
make test
```
