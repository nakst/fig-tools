#define BITMAP_FONT_SIGNATURE (0xF67DD870) // (Random number.)

typedef struct BitmapFontHeader {
	uint32_t signature;
	uint16_t glyphCount;
	uint8_t headerBytes, glyphBytes;
	uint16_t yAscent, yDescent;
	uint16_t xEmWidth;
	uint16_t _unused0;
	// Followed by glyphCount copies of BitmapFontGlyph, sorted by codepoint.
} BitmapFontHeader;

typedef struct BitmapFontKerningEntry {
	uint32_t rightCodepoint;
	int16_t xOffset;
	uint16_t _unused0;
} BitmapFontKerningEntry;

typedef struct BitmapFontGlyph {
	uint32_t bitsOffset; // Stored one row after another; each row is padded to a multiple of 8 bits.
	uint32_t codepoint;
	int16_t xOrigin, yOrigin, xAdvance;
	uint16_t kerningEntryCount; // Stored after the bits. Not necessarily aligned!
	uint16_t bitsWidth, bitsHeight;
} BitmapFontGlyph;
