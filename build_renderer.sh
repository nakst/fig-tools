set -e
mkdir -p bin
g++ textlib/text.cpp -obin/text.o -c -fsanitize=address -g -DUSE_FREETYPE_AND_HARFBUZZ -I/usr/include/freetype2 
g++ renderer.cpp bin/text.o -obin/renderer -fsanitize=address -g -Wall -Wextra -Wno-missing-field-initializers -Wno-sequence-point -lfreetype -lharfbuzz
