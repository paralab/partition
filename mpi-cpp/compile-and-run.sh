mkdir -p build \
&& \
cmake -G Ninja -S . -B build -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
&& \
ninja -C ./build \
&& \
mpirun -np 8 --oversubscribe  ./build/main-dist