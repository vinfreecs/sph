# assuming GLAD files are present as we set up earlier
clang -c deps/glad/src/glad.c -Ideps/glad/include -o build/glad.o

clang++ -std=c++17 test.cpp build/glad.o \
  -Ideps/glad/include -I"$(brew --prefix)/include" \
  -L"$(brew --prefix)/lib" -lglfw \
  -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo \
  -o sph_shaders
