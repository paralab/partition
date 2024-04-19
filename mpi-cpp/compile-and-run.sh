export OMP_NUM_THREADS=8

mkdir -p build \
&& \
cmake -G Ninja -S . -B build -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
&& \
ninja -C ./build \
&& \
mpirun -np 1 --bind-to none  ./build/main-new