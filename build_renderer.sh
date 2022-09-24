set -e
mkdir -p bin
g++ -g textlib/text.cpp -DUSE_FREETYPE_AND_HARFBUZZ -I/usr/include/freetype2 -o bin/text.o -c 
g++ renderer.cpp bin/text.o -obin/renderer -O2 -g -Wall -Wextra -Wno-missing-field-initializers -Wno-sequence-point -lfreetype -lharfbuzz
