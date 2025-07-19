mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release //make sure nvcc and cmake are present
make
./sph
