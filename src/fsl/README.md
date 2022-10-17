
# WASM port of FSL 3.3.11 BET2 and FLIRT

This repository contains the FSL 3.3.11 source code modified to compile with the Emscripten Compiler Frontend (emcc).
The `src` directory only contains projects needed to produce `bet2` and `flirt` binaries.

Main fixes on the codebase to compile with `emcc`:
- Due to `round` ambiguity errors `MISCMATHS::round()` is used with namespace prefix
- The `myFloat * myColumnVector` multiplications don't compile with `emcc`. The intention is to use the `operator*` overload in `NEWMAT::BaseMatrix`. The correct usage is to have the `BaseMatrix` on the lhs and the `float`, `double` or `Real` on the rhs:
`myColumnVector * myFloat`
- `SPMatrix SP(const BaseMatrix& bm1,const BaseMatrix& bm2);` has been added to `newmat.h`, and the usages of `SP` in `miscmaths.cc` has been changed to `NEWMAT::SP`
- Replace `<iostream.h>` with `<iostream>`
- Replace `<fstream.h>` with `<fstream>`
- Fix redefinitions of default arguments errors

## How to build

Set the main FSL environment variables

```bash
export FSLDIR=`pwd`/fsl
. ${FSLDIR}/etc/fslconf/fsl.sh
```

Then check if your machine/compiler is supported by default:

```bash
ls $FSLDIR/config/$FSLMACHTYPE
```

If the above directory does not exist (the ls returns an error), Select the closest match from the directories in $FSLDIR/config and do the following:

```bash
cp -r $FSLDIR/config/closestmatch $FSLDIR/config/$FSLMACHTYPE
```

Go to  `$FSLDIR/config/$FSLMACHTYPE` and change the following variables in `systemvars.mk`:

```bash
# Compiler dependent variables

CC = emcc
CXX = emcc
AR = emar
CSTATICFLAGS = -static
CXXSTATICFLAGS = -static

ARCHFLAGS = -m32

DEPENDFLAGS = -MM
```
The important thing is to to replace the default compilers to `emcc`, and the archiver to `emar`.

Build `bet2.wasm` and `flirt.wasm`

```bash
cd $FSLDIR
./build
```

After the build is finished you can find the `wasm` binaries in `/emscripten` folder.
