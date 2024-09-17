rm -rf ./build
mkdir build
g++ -Wall -c -shared -fPIC ./src/*.cpp -o build/libngineer.so