#!/bin/sh

echo 'Compiling using Cray intel....'
echo 'Creating lib and include folders....'

mkdir results
cd src

mkdir lib
mkdir include

# untar the dependent library 2decomp and fftw
echo 'Untaring....'
rm -rf fftw-3.3.4
rm -rf 2decomp_fft
tar -xvf 2decomp_fft-1.5.847.tar.gz
tar -xvf fftw-3.3.4.tar.gz

echo 'Compiling 2decomp....'
sleep 2
cd 2decomp_fft
mv ./src/Makefile.inc.Cray_XE ./src/Makefile.inc
sed -i '/FFT=fftw3/c\FFT=generic' ./src/Makefile.inc
sed -i '/COMPILER=PGI/c\COMPILER=Intel' ./src/Makefile.inc
make
cp ./lib/lib2decomp_fft.a ../lib/
cp ./include/* ../include/
echo 'Compiling 2decomp Done!'
sleep 5

echo 'Compiling fftw.....'
sleep 5
cd ../fftw-3.3.4
CWD=$(pwd)
./configure CC=icc CFLAGS=-gcc F77=ifort --with-our-malloc16 --enable-threads --with-combined-threads --enable-sse2 --prefix=${CWD}
make
make install
cp ./lib/libfftw3.a ../lib/
cp ./include/fftw3.f03 ../include/
echo 'Compiling fftw Done!'
sleep 5

cd ..
echo 'Compiling HERCULES.....'
cp Makefile.Cray_Intel Makefile
make clean
make
cp hercules.exe ../
echo 'Compiling HERCULES Done!'
sleep 5

