// This file is part of the Essence operating system.
// It is released under the terms of the MIT license -- see LICENSE.md.
// Written by: nakst.

// TODO If the font size is sufficiently large disable subpixel anti-aliasing.
// TODO Variable font support.

/////////////////////////////////////

#include "text.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#define EsCRTexpf expf

typedef uint16_t EsFontFamily;
typedef void *EsHeap;
typedef uint32_t EsStyleID;

struct EsFileStore {
#define FILE_STORE_PATH (2)
	uint8_t type;
	bool operationComplete;
	uint32_t handles;
	uint32_t error;
	const char *cPath;
};

struct EsPoint {
	int32_t x, y;
};

#define EsMessageMutexCheck()

void *EsHeapAllocate(size_t size, bool zeroMemory, EsHeap *heap = nullptr) {
	(void) heap;
	return zeroMemory ? calloc(1, size) : malloc(size);
}

void EsHeapFree(void *address, size_t expectedSize = 0, EsHeap *heap = nullptr) {
	(void) expectedSize;
	(void) heap;
	free(address);
}

void *EsHeapReallocate(void *oldAddress, size_t newAllocationSize, bool zeroNewSpace, EsHeap *heap = nullptr) {
	(void) heap;
	assert(!zeroNewSpace);
	return realloc(oldAddress, newAllocationSize);
}

#define ES_MEMORY_MOVE_BACKWARDS -
#define EsMemoryCopy memcpy
#define EsMemoryCopyReverse memmove
#define EsMemoryZero(destination, bytes) memset(destination, 0, bytes)
#define EsMemoryCompare memcmp

void EsMemoryMove(void *_start, void *_end, intptr_t amount, bool zeroEmptySpace) {
	uint8_t *start = (uint8_t *) _start;
	uint8_t *end = (uint8_t *) _end;
	assert(start <= end);

	if (amount > 0) {
		EsMemoryCopyReverse(start + amount, start, end - start);

		if (zeroEmptySpace) {
			EsMemoryZero(start, amount);
		}
	} else if (amount < 0) {
		EsMemoryCopy(start + amount, start, end - start);

		if (zeroEmptySpace) {
			EsMemoryZero(end + amount, -amount);
		}
	}
}

size_t EsCStringLength(EsCString string) {
	return string ? strlen(string) : 0;
}

int EsStringCompareRaw(const char *s1, ptrdiff_t length1, const char *s2, ptrdiff_t length2) {
	if (length1 == -1) length1 = EsCStringLength(s1);
	if (length2 == -1) length2 = EsCStringLength(s2);
	if (s1 == s2 && length1 == length2) return 0;

	while (length1 || length2) {
		if (!length1) return -1;
		if (!length2) return 1;

		char c1 = *s1;
		char c2 = *s2;

		if (c1 != c2) {
			return c1 - c2;
		}

		length1--;
		length2--;
		s1++;
		s2++;
	}

	return 0;
}

char *EsSystemConfigurationReadString(const char *section, ptrdiff_t sectionBytes, const char *key, ptrdiff_t keyBytes, size_t *valueBytes = nullptr);

inline int Width(EsRectangle rectangle) {
	return rectangle.r - rectangle.l;
}

inline int Height(EsRectangle rectangle) {
	return rectangle.b - rectangle.t;
}

bool EsRectangleClip(EsRectangle parent, EsRectangle rectangle, EsRectangle *output) {
	EsRectangle current = parent;
	EsRectangle intersection;

	if (!((current.l > rectangle.r && current.r > rectangle.l)
			|| (current.t > rectangle.b && current.b > rectangle.t))) {
		intersection.l = current.l > rectangle.l ? current.l : rectangle.l;
		intersection.t = current.t > rectangle.t ? current.t : rectangle.t;
		intersection.r = current.r < rectangle.r ? current.r : rectangle.r;
		intersection.b = current.b < rectangle.b ? current.b : rectangle.b;
	} else {
		intersection = {};
	}

	if (output) {
		*output = intersection;
	}

	return intersection.l < intersection.r && intersection.t < intersection.b;
}

#define EsPanic(x, ...) do { fprintf(stderr, "EsPanic: %s\n", x); assert(false); } while (0)
#define EsAssert(x) assert(x)
#include "bitmap_font.h"
#include "linked_list.cpp"
#include "array.cpp"
#include "unicode.cpp"
#include "hash_table.cpp"

struct EsBuffer {
	union { const uint8_t *in; uint8_t *out; };
	size_t position, bytes;
	void *context;
	EsFileStore *fileStore;
	bool error, canGrow;
};

struct { float scale; } theming = { .scale = 1.0f };

#define ES_MEMORY_MAP_OBJECT_READ_ONLY (1 << 1)
void *EsFileStoreMap(EsFileStore *file, size_t *fileSize, uint32_t flags);
void FileStoreCloseHandle(EsFileStore *fileStore);

void EsDrawBlock(EsPainter *painter, EsRectangle bounds, uint32_t mainColor);
void EsDrawInvert(EsPainter *painter, EsRectangle bounds);

inline int AbsoluteInteger(int a) {
	return a > 0 ? a : -a;
}

#define ES_MACRO_SEARCH(_count, _compar, _result, _found) \
	do { \
		if (_count) { \
			intptr_t low = 0; \
			intptr_t high = _count - 1; \
			\
			while (low <= high) { \
				uintptr_t index = ((high - low) >> 1) + low; \
				int result; \
				_compar \
				\
				if (result < 0) { \
					high = index - 1; \
				} else if (result > 0) { \
					low = index + 1; \
				} else { \
					_result = index; \
					_found = true; \
					break; \
				} \
			} \
			\
			if (high < low) { \
				_result = low; \
				_found = false; \
			} \
		} else { \
			_result = 0; \
			_found = false; \
		} \
	} while (0)
#define EsLiteral(x) (char *) x, EsCStringLength((char *) x)

void EsDrawInvert(EsPainter *painter, EsRectangle bounds) {
	EsPaintTarget *target = painter->target;
	if (!EsRectangleClip(bounds, painter->clip, &bounds)) return;

	for (int y = bounds.t; y < bounds.b; y++) {
		for (int x = bounds.l; x < bounds.r; x++) {
			((uint32_t *) target->bits)[bounds.t * target->stride / 4 + bounds.l] ^= 0xFFFFFF;
		}
	}
}

void EsDrawBlock(EsPainter *painter, EsRectangle bounds, uint32_t color) {
	EsPaintTarget *target = painter->target;
	if (!EsRectangleClip(bounds, painter->clip, &bounds)) return;

	for (int y = bounds.t; y < bounds.b; y++) {
		for (int x = bounds.l; x < bounds.r; x++) {
			((uint32_t *) target->bits)[bounds.t * target->stride / 4 + bounds.l] = color;
		}
	}
}

void *EsFileStoreMap(EsFileStore *file, size_t *fileSize, uint32_t flags) {
	if (file->type == FILE_STORE_PATH) {
		int fd = open(file->cPath, O_RDONLY);
		struct stat statResult;
		fstat(fd, &statResult);
		*fileSize = statResult.st_size;
		return mmap(NULL, statResult.st_size, PROT_READ, MAP_SHARED, fd, 0);
	} else {
		EsAssert(false);
		return nullptr;
	}
}

void FileStoreCloseHandle(EsFileStore *fileStore) {
	EsAssert(fileStore->handles < 0x80000000);
	if (--fileStore->handles) return;
	EsHeapFree(fileStore);
}

/////////////////////////////////////

#ifdef USE_FREETYPE_AND_HARFBUZZ
#include <harfbuzz/hb.h>
#include <harfbuzz/hb-ft.h>
#ifndef FT_EXPORT
#define FT_EXPORT(x) extern "C" x
#endif
#include <ft2build.h>
#include FT_FREETYPE_H
#include <freetype/ftoutln.h>
#endif

#define CHARACTER_SUBPIXEL (1) // 24 bits per pixel; each byte specifies the alpha of each RGB channel.
#define CHARACTER_IMAGE    (2) // 32 bits per pixel, ARGB.
#define CHARACTER_RECOLOR  (3) // 32 bits per pixel, AXXX.

#define FREETYPE_UNIT_SCALE (64)

#define FALLBACK_SCRIPT_LANGUAGE ("en")
#define FALLBACK_SCRIPT (0x4C61746E) // "Latn"

struct Font {
#define FONT_TYPE_BITMAP (1)
#define FONT_TYPE_FREETYPE_AND_HARFBUZZ (2)
	uintptr_t type;

	union {
		const void *bitmapData; // All data has been validated to load the font.

#ifdef USE_FREETYPE_AND_HARFBUZZ
		struct {
			FT_Face ft;
			hb_font_t *hb;
		};
#endif
	};
};

struct GlyphCacheKey {
	uint32_t glyphIndex;
	uint16_t size;
	uint16_t fractionalPosition;
	Font font;
};

struct GlyphCacheEntry {
	uint8_t *data;
	size_t dataBytes;
	int width, height, xoff, yoff;
	int type;

	LinkedItem<GlyphCacheEntry> itemLRU;
	GlyphCacheKey key;
};

struct FontSubstitutionKey {
	EsFontFamily family;
	uint16_t _unused0;
	uint32_t script;
};

struct FontDatabaseEntry : EsFontInformation {
	EsFileStore *files[18];
	char *scripts;
	size_t scriptsBytes;
};

enum TextStyleDifference {
	TEXT_STYLE_NEW_FONT,   // A new font is selected.
	TEXT_STYLE_NEW_SHAPE,  // Shaping parameters have changed.
	TEXT_STYLE_NEW_RENDER, // Render-only properties have changed.
	TEXT_STYLE_IDENTICAL,  // The styles are the same.
};	

struct TextPiece {
	// Shaped glyphs, on the same line, and with constant style and script.
	int32_t ascent, descent, width;
	const EsTextStyle *style;
	uintptr_t glyphOffset;
	size_t glyphCount;
	uintptr_t start, end;
	bool isTabPiece, isEllipsisPiece;
};

struct TextLine {
	int32_t ascent, descent, width;
	bool hasEllipsis;
	uintptr_t ellipsisPieceIndex;
	uintptr_t pieceOffset;
	size_t pieceCount;
};

struct TextRun {
	EsTextStyle style;
	uint32_t offset;
	uint32_t script;
};

#ifdef USE_FREETYPE_AND_HARFBUZZ
typedef hb_glyph_info_t TextGlyphInfo;
typedef hb_glyph_position_t TextGlyphPosition;
typedef hb_segment_properties_t TextSegmentProperties;
typedef hb_buffer_t TextShapeBuffer;
typedef hb_feature_t TextFeature;
typedef hb_script_t TextScript;
#else
struct TextGlyphInfo {
	uint32_t codepoint;
	uint32_t cluster;
};

struct TextGlyphPosition {
	int32_t x_advance;
	int32_t y_advance;
	int32_t x_offset;
	int32_t y_offset;
};

struct TextSegmentProperties {
	uint32_t direction;
	uint32_t script;
	uint32_t language;
};

struct TextShapeBuffer {
	uint8_t _unused0;
};

struct TextFeature {
	uint8_t _unused0;
};

typedef uint32_t TextScript;
#endif

struct EsTextPlan {
	TextShapeBuffer *buffer;
	TextSegmentProperties segmentProperties;

	const char *string; 

	Array<TextRun> textRuns; 
	uintptr_t textRunPosition;

	const EsTextStyle *currentTextStyle;
	Font font;

	BreakState breaker;

	Array<TextGlyphInfo> glyphInfos;
	Array<TextGlyphPosition> glyphPositions;

	Array<TextPiece> pieces;
	Array<TextLine> lines;

	int32_t totalHeight, totalWidth;

	bool singleUse;

	EsTextPlanProperties properties;
};

struct {
	// Database.
	HashStore<FontSubstitutionKey, EsFontFamily> substitutions;
	Array<FontDatabaseEntry> database;
	uintptr_t sans, serif, monospaced, fallback;
	char *sansName, *serifName, *monospacedName, *fallbackName;

	// Rendering.
#ifdef USE_FREETYPE_AND_HARFBUZZ
	FT_Library freetypeLibrary;
#endif

	// Caching.
	HashStore<EsFont, Font> loaded; // TODO How many fonts to keep loaded? Reference counting?
#define GLYPH_CACHE_MAX_SIZE (4194304)
	HashStore<GlyphCacheKey, GlyphCacheEntry *> glyphCache;
	LinkedList<GlyphCacheEntry> glyphCacheLRU;
	size_t glyphCacheBytes;
} fontManagement;

struct {
	EsBuffer pack;
	const uint8_t *standardPack;
	size_t standardPackSize;
	char *buffer;
	size_t bufferPosition, bufferAllocated;
} iconManagement;

// --------------------------------- Glyph cache.

void GlyphCacheFreeEntry() {
	GlyphCacheEntry *entry = fontManagement.glyphCacheLRU.lastItem->thisItem;
	fontManagement.glyphCacheLRU.Remove(&entry->itemLRU);
	fontManagement.glyphCache.Delete(&entry->key);
	EsAssert(fontManagement.glyphCacheBytes >= entry->dataBytes);
	fontManagement.glyphCacheBytes -= entry->dataBytes;
	EsHeapFree(entry->data);
	EsHeapFree(entry);
}

void RegisterGlyphCacheEntry(GlyphCacheKey key, GlyphCacheEntry *entry) {
	// Free space in the glyph cache.
	// Do this before adding the new glyph to the cache,
	// in case the new glyph doesn't fit in the cache at all.
	while (fontManagement.glyphCacheBytes > GLYPH_CACHE_MAX_SIZE) {
		GlyphCacheFreeEntry();
	}

	entry->itemLRU.thisItem = entry;
	entry->key = key;
	*fontManagement.glyphCache.Put(&key) = entry;
	fontManagement.glyphCacheLRU.InsertStart(&entry->itemLRU);
	fontManagement.glyphCacheBytes += entry->dataBytes;
}

GlyphCacheEntry *LookupGlyphCacheEntry(GlyphCacheKey key) {
	GlyphCacheEntry *entry = fontManagement.glyphCache.Get1(&key);

	if (!entry) {
		return (GlyphCacheEntry *) EsHeapAllocate(sizeof(GlyphCacheEntry), true);
	} else {
		fontManagement.glyphCacheLRU.Remove(&entry->itemLRU);
		fontManagement.glyphCacheLRU.InsertStart(&entry->itemLRU);
		return entry;
	}
}

// --------------------------------- Font backend abstraction layer.

bool FontLoad(Font *font, const void *data, size_t dataBytes) {
	if (dataBytes > sizeof(BitmapFontHeader)) {
		const BitmapFontHeader *header = (const BitmapFontHeader *) data;

		if (header->signature == BITMAP_FONT_SIGNATURE && (size_t) header->headerBytes >= sizeof(BitmapFontHeader) 
				&& (size_t) header->headerBytes < 0x80 && (size_t) header->glyphBytes < 0x80 && (size_t) header->glyphCount < 0x8000 
				&& dataBytes > (size_t) header->headerBytes + (size_t) header->glyphCount * (size_t) header->glyphBytes
				&& header->glyphCount >= 1 /* index 0 is used as a fallback glyph, which must be present */) {
			for (uintptr_t i = 0; i < header->glyphCount; i++) {
				const BitmapFontGlyph *glyph = ((const BitmapFontGlyph *) ((const uint8_t *) data + header->headerBytes + i * header->glyphBytes));

				size_t bytesPerRow = (glyph->bitsWidth + 7) / 8;
				size_t bitsStorage = (size_t) glyph->bitsHeight * bytesPerRow; 
				size_t kerningStorage = (size_t) glyph->kerningEntryCount * sizeof(BitmapFontKerningEntry);

				if (glyph->bitsWidth > 0x4000 || glyph->bitsHeight > 0x4000 || glyph->kerningEntryCount > 0x4000
						|| glyph->bitsOffset >= dataBytes
						|| bitsStorage > dataBytes - glyph->bitsOffset
						|| kerningStorage > dataBytes - glyph->bitsOffset - bitsStorage) {
					return false;
				}
			}

			font->bitmapData = data;
			font->type = FONT_TYPE_BITMAP;
			return true;
		}
	}

#ifdef USE_FREETYPE_AND_HARFBUZZ
	if (!fontManagement.freetypeLibrary) {
		FT_Init_FreeType(&fontManagement.freetypeLibrary);
	}

	if (!FT_New_Memory_Face(fontManagement.freetypeLibrary, (uint8_t *) data, dataBytes, 0, &font->ft)) {
		font->hb = hb_ft_font_create(font->ft, nullptr);

		if (font->hb) {
			font->type = FONT_TYPE_FREETYPE_AND_HARFBUZZ;
			return true;
		} else {
			FT_Done_Face(font->ft);
		}
	}
#endif

	return false;
}

void FontSetSize(Font *font, uint32_t size) {
	if (font->type == FONT_TYPE_BITMAP) {
		(void) size;
		// Ignored.
	}

#ifdef USE_FREETYPE_AND_HARFBUZZ
	if (font->type == FONT_TYPE_FREETYPE_AND_HARFBUZZ) {
		FT_Set_Char_Size(font->ft, 0, size, 72, 72);
		hb_ft_font_changed(font->hb);
	}
#endif
}

int32_t FontGetAscent(Font *font) {
	if (font->type == FONT_TYPE_BITMAP) {
		const BitmapFontHeader *header = (const BitmapFontHeader *) font->bitmapData;
		return header->yAscent * FREETYPE_UNIT_SCALE;
	}

#ifdef USE_FREETYPE_AND_HARFBUZZ
	if (font->type == FONT_TYPE_FREETYPE_AND_HARFBUZZ) {
		return font->ft->size->metrics.ascender;
	}
#endif

	return 0;
}

int32_t FontGetDescent(Font *font) {
	if (font->type == FONT_TYPE_BITMAP) {
		const BitmapFontHeader *header = (const BitmapFontHeader *) font->bitmapData;
		return header->yDescent * -FREETYPE_UNIT_SCALE;
	}

#ifdef USE_FREETYPE_AND_HARFBUZZ
	if (font->type == FONT_TYPE_FREETYPE_AND_HARFBUZZ) {
		return font->ft->size->metrics.descender;
	}
#endif

	return 0;
}

int32_t FontGetEmWidth(Font *font) {
	if (font->type == FONT_TYPE_BITMAP) {
		const BitmapFontHeader *header = (const BitmapFontHeader *) font->bitmapData;
		return header->xEmWidth;
	}

#ifdef USE_FREETYPE_AND_HARFBUZZ
	if (font->type == FONT_TYPE_FREETYPE_AND_HARFBUZZ) {
		return font->ft->size->metrics.x_ppem;
	}
#endif

	return 0;
}

bool FontRenderGlyph(GlyphCacheKey key, GlyphCacheEntry *entry) {
	if (key.font.type == FONT_TYPE_BITMAP) {
		const BitmapFontHeader *header = (const BitmapFontHeader *) key.font.bitmapData;
		const BitmapFontGlyph *glyph = ((const BitmapFontGlyph *) ((const uint8_t *) key.font.bitmapData + header->headerBytes + key.glyphIndex * header->glyphBytes));

		entry->width = glyph->bitsWidth;
		entry->height = glyph->bitsHeight;
		entry->xoff = -glyph->xOrigin;
		entry->yoff = -glyph->yOrigin;
		entry->dataBytes = sizeof(uint8_t) * 4 * glyph->bitsWidth * glyph->bitsHeight;
		entry->data = (uint8_t *) EsHeapAllocate(entry->dataBytes, false);

		if (!entry->data) {
			return false;
		}

		size_t bytesPerRow = (glyph->bitsWidth + 7) / 8;

		for (uintptr_t i = 0; i < glyph->bitsHeight; i++) {
			const uint8_t *row = (const uint8_t *) key.font.bitmapData + glyph->bitsOffset + bytesPerRow * i;

			for (uintptr_t j = 0; j < glyph->bitsWidth; j++) {
				// TODO More efficient storage.
				uint8_t byte = (row[j / 8] & (1 << (j % 8))) ? 0xFF : 0x00;
				uint32_t copy = (byte << 24) | (byte << 16) | (byte << 8) | byte;
				((uint32_t *) entry->data)[i * entry->width + j] = copy;
			}
		}

		return true;
	}

#ifdef USE_FREETYPE_AND_HARFBUZZ
	if (key.font.type == FONT_TYPE_FREETYPE_AND_HARFBUZZ) {
		FT_Load_Glyph(key.font.ft, key.glyphIndex, FT_LOAD_DEFAULT);
		FT_Outline_Translate(&key.font.ft->glyph->outline, key.fractionalPosition, 0);

		int width; 
		int height; 
		int xoff; 
		int yoff; 
		uint8_t *output;

#ifdef USE_LCD_RENDERING
		FT_Render_Glyph(key.font.ft->glyph, FT_RENDER_MODE_LCD);
#else
		FT_Render_Glyph(key.font.ft->glyph, FT_RENDER_MODE_NORMAL);
#endif

		FT_Bitmap *bitmap = &key.font.ft->glyph->bitmap;
#ifdef USE_LCD_RENDERING
		width = bitmap->width / 3;
#else
		width = bitmap->width;
#endif
		height = bitmap->rows;
		xoff = key.font.ft->glyph->bitmap_left;
		yoff = -key.font.ft->glyph->bitmap_top;

		entry->dataBytes = 1 /*stupid hack for whitespace*/ + width * height * 4;
		output = (uint8_t *) EsHeapAllocate(entry->dataBytes, false);

		if (!output) {
			return false;
		}

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
#ifdef USE_LCD_RENDERING
				int32_t r = (int32_t) ((uint8_t *) bitmap->buffer)[x * 3 + y * bitmap->pitch + 0];
				int32_t g = (int32_t) ((uint8_t *) bitmap->buffer)[x * 3 + y * bitmap->pitch + 1];
				int32_t b = (int32_t) ((uint8_t *) bitmap->buffer)[x * 3 + y * bitmap->pitch + 2];

				// Reduce how noticible the colour fringes are.
				// TODO Make this adjustable?
				int32_t average = (r + g + b) / 3;
				r -= (r - average) / 3;
				g -= (g - average) / 3;
				b -= (b - average) / 3;

				output[(x + y * width) * 4 + 0] = (uint8_t) r;
				output[(x + y * width) * 4 + 1] = (uint8_t) g;
				output[(x + y * width) * 4 + 2] = (uint8_t) b;
				output[(x + y * width) * 4 + 3] = 0xFF;
#else
				int32_t a = (int32_t) ((uint8_t *) bitmap->buffer)[x + y * bitmap->pitch];
				output[(x + y * width) * 4 + 0] = (uint8_t) a;
				output[(x + y * width) * 4 + 1] = (uint8_t) a;
				output[(x + y * width) * 4 + 2] = (uint8_t) a;
				output[(x + y * width) * 4 + 3] = 0xFF;
#endif
			}
		}

		if (output) {
			entry->data = output;
			entry->width = width;
			entry->height = height;
			entry->xoff = xoff;
			entry->yoff = yoff;
			return true;
		}

		return false;
	}
#endif

	return false;
}

void FontShapeTextDone(EsTextPlan *plan, uint32_t glyphCount, TextGlyphInfo *_glyphInfos, TextGlyphPosition *_glyphPositions) {
	if (plan->font.type == FONT_TYPE_BITMAP) {
		(void) glyphCount;
		Array<TextGlyphInfo> glyphInfos = { .array = _glyphInfos };
		Array<TextGlyphPosition> glyphPositions = { .array = _glyphPositions };
		glyphInfos.Free();
		glyphPositions.Free();
	}
}

void FontShapeText(EsTextPlan *plan, const char *string, size_t stringBytes, 
		uintptr_t sectionOffsetBytes, size_t sectionCountBytes,
		TextFeature *features, size_t featureCount,
		uint32_t *glyphCount, TextGlyphInfo **_glyphInfos, TextGlyphPosition **_glyphPositions) {
	if (plan->font.type == FONT_TYPE_BITMAP) {
		(void) features;
		(void) featureCount;
		(void) stringBytes;

		Array<TextGlyphInfo> glyphInfos = {};
		Array<TextGlyphPosition> glyphPositions = {};

		Font *font = &plan->font;
		const BitmapFontHeader *header = (const BitmapFontHeader *) font->bitmapData;

		const char *text = string + sectionOffsetBytes;

		while (text < string + sectionOffsetBytes + sectionCountBytes) {
			TextGlyphInfo info = {};
			TextGlyphPosition position = {};
			info.cluster = text - (string + sectionOffsetBytes);
			uint32_t codepoint = utf8_value(text, string + sectionOffsetBytes + sectionCountBytes - text, nullptr);
			if (!codepoint) break;
			text = utf8_advance(text);

			bool glyphFound = false;
			uint32_t glyphIndex = 0;
			ES_MACRO_SEARCH(header->glyphCount, result = codepoint - ((const BitmapFontGlyph *) ((const uint8_t *) 
							font->bitmapData + header->headerBytes + index * header->glyphBytes))->codepoint;, glyphIndex, glyphFound);
			info.codepoint = glyphFound ? glyphIndex : 0;

			const BitmapFontGlyph *glyph = ((const BitmapFontGlyph *) ((const uint8_t *) 
						font->bitmapData + header->headerBytes + info.codepoint * header->glyphBytes));
			size_t bytesPerRow = (glyph->bitsWidth + 7) / 8;
			position.x_advance = glyph->xAdvance * FREETYPE_UNIT_SCALE;

			if (glyph->kerningEntryCount && text < string + sectionOffsetBytes + sectionCountBytes) {
				uint32_t nextCodepoint = utf8_value(text, string + sectionOffsetBytes + sectionCountBytes - text, nullptr);

				for (uintptr_t i = 0; i < glyph->kerningEntryCount; i++) {
					BitmapFontKerningEntry entry;
					EsMemoryCopy(&entry, (const uint8_t *) font->bitmapData + glyph->bitsOffset 
							+ bytesPerRow * glyph->bitsHeight + sizeof(BitmapFontKerningEntry) * i, sizeof(BitmapFontKerningEntry));

					if (entry.rightCodepoint == nextCodepoint) {
						position.x_advance += entry.xOffset * FREETYPE_UNIT_SCALE;
						break;
					}
				}
			}

			glyphInfos.Add(info);
			glyphPositions.Add(position);
		}

		*glyphCount = glyphInfos.Length();
		*_glyphInfos = &glyphInfos[0];
		*_glyphPositions = &glyphPositions[0];

		return;
	}

#ifdef USE_FREETYPE_AND_HARFBUZZ
	if (plan->font.type == FONT_TYPE_FREETYPE_AND_HARFBUZZ) {
		hb_buffer_clear_contents(plan->buffer);
		hb_buffer_set_segment_properties(plan->buffer, &plan->segmentProperties);
		hb_buffer_add_utf8(plan->buffer, string, stringBytes, sectionOffsetBytes, sectionCountBytes);
		hb_shape(plan->font.hb, plan->buffer, features, featureCount);
		*_glyphInfos = hb_buffer_get_glyph_infos(plan->buffer, glyphCount);
		*_glyphPositions = hb_buffer_get_glyph_positions(plan->buffer, glyphCount);
		return;
	}
#endif
}

uint32_t FontGetScriptFromCodepoint(uint32_t codepoint, bool *inheritingScript) {
#ifdef USE_FREETYPE_AND_HARFBUZZ
	static hb_unicode_funcs_t *unicodeFunctions = nullptr;

	if (!unicodeFunctions) {
		// Multiple threads could call this at the same time, but it doesn't matter,
		// since they should always return the same thing anyway...
		unicodeFunctions = hb_unicode_funcs_get_default();
	}

	uint32_t script = hb_unicode_script(unicodeFunctions, codepoint);
	*inheritingScript = script == HB_SCRIPT_COMMON || script == HB_SCRIPT_INHERITED;
	return script;
#else
	(void) codepoint;
	*inheritingScript = false;
	return FALLBACK_SCRIPT;
#endif
}

void FontInitialiseShaping(EsTextPlan *plan) {
#ifdef USE_FREETYPE_AND_HARFBUZZ
	plan->buffer = hb_buffer_create();
	hb_buffer_set_cluster_level(plan->buffer, HB_BUFFER_CLUSTER_LEVEL_MONOTONE_CHARACTERS);
	plan->segmentProperties.direction = (plan->properties.flags & ES_TEXT_PLAN_RTL) ? HB_DIRECTION_RTL : HB_DIRECTION_LTR;
	plan->segmentProperties.script = (TextScript) FALLBACK_SCRIPT;
	plan->segmentProperties.language = hb_language_from_string(plan->properties.cLanguage ?: FALLBACK_SCRIPT_LANGUAGE, -1);
#else
	(void) plan;
#endif
}

void FontDestroyShaping(EsTextPlan *plan) {
#ifdef USE_FREETYPE_AND_HARFBUZZ
	hb_buffer_destroy(plan->buffer);
	plan->buffer = nullptr;
#else
	(void) plan;
#endif
}

// --------------------------------- Font management.

void FontInitialise() {
	if (fontManagement.database.Length()) {
		return;
	}
	
#if 0
	fontManagement.sansName       = EsSystemConfigurationReadString(EsLiteral("ui_fonts"), EsLiteral("sans"));
	fontManagement.serifName      = EsSystemConfigurationReadString(EsLiteral("ui_fonts"), EsLiteral("serif"));
	fontManagement.monospacedName = EsSystemConfigurationReadString(EsLiteral("ui_fonts"), EsLiteral("mono"));
	fontManagement.fallbackName   = EsSystemConfigurationReadString(EsLiteral("ui_fonts"), EsLiteral("fallback")); 

	FontDatabaseEntry nullFont = {};
	fontManagement.database.Add(nullFont);

	EsMutexAcquire(&api.systemConfigurationMutex);

	for (uintptr_t i = 0; i < api.systemConfigurationGroups.Length(); i++) {
		SystemConfigurationGroup *g = &api.systemConfigurationGroups[i];

		if (0 == EsStringCompareRaw(g->section, g->sectionBytes, EsLiteral("font"))) {
			FontDatabaseEntry entry = {};
			const char *name = nullptr;
			size_t nameBytes = 0;

			for (uintptr_t i = 0; i < g->itemCount; i++) {
				SystemConfigurationItem *item = g->items + i;

				if (0 == EsStringCompareRaw(item->key, item->keyBytes, EsLiteral("name"))) {
					name = item->value, nameBytes = item->valueBytes;
				} else if (0 == EsStringCompareRaw(item->key, item->keyBytes, EsLiteral("category"))) {
					entry.categoryBytes = MinimumInteger(item->valueBytes, sizeof(entry.category));
					EsMemoryCopy(entry.category, item->value, entry.categoryBytes);
				} else if (0 == EsStringCompareRaw(item->key, item->keyBytes, EsLiteral("scripts"))) {
					entry.scripts = item->value;
					entry.scriptsBytes = item->valueBytes;
				} else if ((item->keyBytes == 2 && item->key[0] == '.' && EsCRTisdigit(item->key[1]))
						|| (item->keyBytes == 3 && item->key[0] == '.' && EsCRTisdigit(item->key[1]) && item->key[2] == 'i') ) {
					int weight = item->key[1] - '0';
					bool italic = item->keyBytes == 3;

					if (italic) {
						entry.availableWeightsItalic |= 1 << weight;
					} else {
						entry.availableWeightsNormal |= 1 << weight;
					}

					size_t fileIndex = weight - 1 + italic * 9;

					if (item->valueBytes && item->value[0] == ':') {
						entry.files[fileIndex] = FileStoreCreateFromEmbeddedFile(&bundleDesktop, item->value + 1, item->valueBytes - 1);
					} else {
						entry.files[fileIndex] = FileStoreCreateFromPath(item->value, item->valueBytes);
					}
				}
			}

			if (0 == EsStringCompareRaw(name, nameBytes, EsLiteral(fontManagement.sansName))) {
				fontManagement.sans = fontManagement.database.Length();
			}

			if (0 == EsStringCompareRaw(name, nameBytes, EsLiteral(fontManagement.serifName))) {
				fontManagement.serif = fontManagement.database.Length();
			}

			if (0 == EsStringCompareRaw(name, nameBytes, EsLiteral(fontManagement.monospacedName))) {
				fontManagement.monospaced = fontManagement.database.Length();
			}

			if (0 == EsStringCompareRaw(name, nameBytes, EsLiteral(fontManagement.fallbackName))) {
				fontManagement.fallback = fontManagement.database.Length();
			}

			entry.nameBytes = MinimumInteger(nameBytes, sizeof(entry.name));
			EsMemoryCopy(entry.name, name, entry.nameBytes);
			entry.id = fontManagement.database.Length();

			fontManagement.database.Add(entry);
		}
	}

	EsMutexRelease(&api.systemConfigurationMutex);
#else
	fontManagement.sansName       = strdup("sans");
	fontManagement.serifName      = strdup("serif");
	fontManagement.monospacedName = strdup("mono");
	fontManagement.fallbackName   = strdup("fallback"); 

	FontDatabaseEntry nullFont = {};
	fontManagement.database.Add(nullFont);
#endif
}

void TextAddFont(const char *cName, const char *cCategory, const char *cScripts, const char *const *const cPaths, const char *const *const cTypes, size_t fileCount) {
	FontInitialise();

	size_t nameBytes = EsCStringLength(cName);
	size_t categoryBytes = EsCStringLength(cCategory);
	size_t scriptsBytes = EsCStringLength(cScripts);

	FontDatabaseEntry entry = {};
	EsAssert(nameBytes < sizeof(entry.name));
	EsAssert(categoryBytes < sizeof(entry.category));
	entry.nameBytes = nameBytes;
	entry.categoryBytes = categoryBytes;
	EsMemoryCopy(entry.name, cName, entry.nameBytes);
	EsMemoryCopy(entry.category, cCategory, categoryBytes);
	entry.scripts = (char *) cScripts;
	entry.scriptsBytes = scriptsBytes;

	for (uintptr_t i = 0; i < fileCount; i++) {
		int weight = cTypes[i][0] - '0';
		bool italic = cTypes[i][EsCStringLength(cTypes[i]) - 1] == 'i';

		if (italic) {
			entry.availableWeightsItalic |= 1 << weight;
		} else {
			entry.availableWeightsNormal |= 1 << weight;
		}

		EsFileStore *store = (EsFileStore *) EsHeapAllocate(sizeof(EsFileStore) + EsCStringLength(cPaths[i]) + 1, true);
		store->type = FILE_STORE_PATH;
		store->handles = 1;
		store->cPath = (const char *) (store + 1);
		EsMemoryCopy((char *) (store + 1), cPaths[i], EsCStringLength(cPaths[i]) + 1);

		size_t fileIndex = weight - 1 + italic * 9;
		entry.files[fileIndex] = store;
	}

	if (0 == EsStringCompareRaw(cName, nameBytes, EsLiteral(fontManagement.sansName)))       fontManagement.sans       = fontManagement.database.Length();
	if (0 == EsStringCompareRaw(cName, nameBytes, EsLiteral(fontManagement.serifName)))      fontManagement.serif      = fontManagement.database.Length();
	if (0 == EsStringCompareRaw(cName, nameBytes, EsLiteral(fontManagement.monospacedName))) fontManagement.monospaced = fontManagement.database.Length();
	if (0 == EsStringCompareRaw(cName, nameBytes, EsLiteral(fontManagement.fallbackName)))   fontManagement.fallback   = fontManagement.database.Length();

	entry.id = fontManagement.database.Length();
	fontManagement.database.Add(entry);
}

EsFontFamily FontGetStandardFamily(EsFontFamily family) {
	FontInitialise();

	if (family == 0 || family == ES_FONT_SANS) {
		return fontManagement.sans ?: fontManagement.fallback;
	} else if (family == ES_FONT_SERIF) {
		return fontManagement.serif ?: fontManagement.fallback;
	} else if (family == ES_FONT_MONOSPACED) {
		return fontManagement.monospaced ?: fontManagement.fallback;
	} else {
		return family;
	}
}

bool EsFontDatabaseLookupByName(const char *name, ptrdiff_t nameBytes, EsFontInformation *information) {
	FontInitialise();
	EsMemoryZero(information, sizeof(EsFontInformation));

	for (uintptr_t i = 0; i < fontManagement.database.Length(); i++) {
		if (0 == EsStringCompareRaw(name, nameBytes, fontManagement.database[i].name, fontManagement.database[i].nameBytes)) {
			EsMemoryCopy(information, &fontManagement.database[i], sizeof(EsFontInformation));
			return true;
		}
	}

	return false;
}

bool EsFontDatabaseLookupByID(EsFontFamily id, EsFontInformation *information) {
	FontInitialise();
	EsMemoryZero(information, sizeof(EsFontInformation));

	id = FontGetStandardFamily(id);

	if (id >= fontManagement.database.Length()) {
		return false;
	}

	EsMemoryCopy(information, &fontManagement.database[id], sizeof(EsFontInformation));
	
	return true;
}

EsFontFamily EsFontDatabaseInsertFile(const EsFontInformation *information, EsFileStore *store) {
	// TODO Locking.

	FontInitialise();
	EsAssert(store->handles);
	FontDatabaseEntry *entry = nullptr;

	if (information->nameBytes) {
		for (uintptr_t i = 1; i < fontManagement.database.Length(); i++) {
			FontDatabaseEntry *entry = &fontManagement.database[i];
			EsAssert(entry->id == i);

			if (0 == EsStringCompareRaw(information->name, information->nameBytes, 
						fontManagement.database[i].name, fontManagement.database[i].nameBytes)) {
				if ((information->availableWeightsItalic & entry->availableWeightsItalic)
						|| (information->availableWeightsNormal & entry->availableWeightsNormal)) {
					// The variant is already in the database.
					return entry->id;
				}

				goto addFileToFamily;
			}
		}
	}

	{
		// The family is not yet in the database; add it.
		FontDatabaseEntry e = {};
		EsMemoryCopy(&e, information, sizeof(EsFontInformation));
		e.id = fontManagement.database.Length();

		if (fontManagement.database.Add(e)) {
			entry = &fontManagement.database.Last();
		} else {
			return 0;
		}
	}

	addFileToFamily:;

	store->handles++;

	entry->availableWeightsNormal |= information->availableWeightsNormal;
	entry->availableWeightsItalic |= information->availableWeightsItalic;

	for (uintptr_t i = 0; i < 18; i++) {
		if ((i < 9 && (information->availableWeightsNormal & (1 << i))) 
				|| (i >= 9 && (information->availableWeightsItalic & (1 << i)))) {
			store->handles++;
			entry->files[i] = store;
			return entry->id;
		}
	}

	EsAssert(false);
	return 0;
}

EsFontInformation *EsFontDatabaseEnumerate(size_t *count) {
	// TODO Locking.
	FontInitialise();
	*count = 0;
	EsFontInformation *result = (EsFontInformation *) EsHeapAllocate(sizeof(EsFontInformation) * (fontManagement.database.Length() - 1), true);
	if (!result) return nullptr;
	*count = fontManagement.database.Length() - 1;
	for (uintptr_t i = 1; i <= *count; i++) EsFontDatabaseLookupByID(i, &result[i - 1]);
	return result;
}

bool FontSupportsScript(FontDatabaseEntry *entry, uint32_t _script, bool first) {
	if (!entry->scriptsBytes) {
		return first;
	}

	char script[4];
	script[0] = (char) (_script >> 24);
	script[1] = (char) (_script >> 16);
	script[2] = (char) (_script >>  8);
	script[3] = (char) (_script >>  0);

	for (uintptr_t i = 0; i <= entry->scriptsBytes - 4; i += 5) {
		if (script[0] == entry->scripts[i + 0] 
				&& script[1] == entry->scripts[i + 1] 
				&& script[2] == entry->scripts[i + 2] 
				&& script[3] == entry->scripts[i + 3]) {
			return true;
		}
	}

	return false;
}

EsFontFamily FontApplySubstitution(EsTextPlanProperties *properties, EsFontFamily family, uint32_t script) {
	FontInitialise();

	if (properties->flags & ES_TEXT_PLAN_NO_FONT_SUBSTITUTION) {
		return family;
	}

	FontSubstitutionKey key = {};
	key.family = FontGetStandardFamily(family);
	key.script = script;
	EsFontFamily result = fontManagement.substitutions.Get1(&key);
	if (result) return result;

	EsAssert(key.family < fontManagement.database.Length());
	FontDatabaseEntry *entry = &fontManagement.database[key.family];

	if (FontSupportsScript(entry, script, true)) {
		*fontManagement.substitutions.Put(&key) = key.family;
		return key.family;
	}

	EsFontFamily firstMatch = (EsFontFamily) -1;

	for (uintptr_t i = 1; i < fontManagement.database.Length(); i++) {
		if (&fontManagement.database[i] == entry) continue;
		if (!FontSupportsScript(&fontManagement.database[i], script, false)) continue;

		if (firstMatch == (EsFontFamily) -1) {
			firstMatch = i;
		}

		if (0 == EsStringCompareRaw(fontManagement.database[i].category, fontManagement.database[i].categoryBytes,
					entry->category, entry->categoryBytes)) {
			*fontManagement.substitutions.Put(&key) = i;
			return i;
		}
	}

	if (firstMatch != (EsFontFamily) -1) {
		*fontManagement.substitutions.Put(&key) = firstMatch;
		return firstMatch;
	} else {
		// No installed font supports the script.
		*fontManagement.substitutions.Put(&key) = key.family;
		return result;
	}
}

Font FontGet(EsFont key) {
	FontInitialise();

	if (key.weight == 0) {
		key.weight = ES_FONT_REGULAR;
	}

	key.family = FontGetStandardFamily(key.family);

	Font *_font = fontManagement.loaded.Get(&key);
	if (_font) return *_font;

	EsFileStore *file = nullptr;
	int matchDistance = 1000;

	EsAssert(key.family < fontManagement.database.Length());
	FontDatabaseEntry *entry = &fontManagement.database[key.family];

	for (uintptr_t i = 0; i < 18; i++) {
		if (entry->files[i]) {
			int weight = (i % 9) + 1;
			bool italic = i >= 9;
			int distance = ((italic != ((key.flags & ES_FONT_ITALIC) ? true : false)) ? 10 : 0) + AbsoluteInteger(weight - key.weight);

			if (distance < matchDistance) {
				matchDistance = distance;
				file = entry->files[i];
			}
		}
	}

	if (!file) {
		// EsPrint("Could not load font (f%d/w%d/%X).\n", key.family, key.weight, key.flags);
		key.family = fontManagement.fallback;
		return FontGet(key);
	}

	// EsPrint("Loading font from '%z' (f%d/w%d/i%d).\n", file, key.family, key.weight, key.italic);

	size_t size;
	void *data = EsFileStoreMap(file, &size, ES_MEMORY_MAP_OBJECT_READ_ONLY);

	if (!data) {
		// EsPrint("Could not load font (f%d/w%d/%X).\n", key.family, key.weight, key.flags);
		key.family = fontManagement.fallback;
		return FontGet(key);
	}

	Font font = {};

	if (!FontLoad(&font, data, size)) {
		// EsPrint("Could not load font (f%d/w%d/%X).\n", key.family, key.weight, key.flags);
		key.family = fontManagement.fallback;
		return FontGet(key);
	}

	*fontManagement.loaded.Put(&key) = font;
	return font;
}

void FontDatabaseFree() {
	while (fontManagement.glyphCacheLRU.count) {
		GlyphCacheFreeEntry();
	}

	for (uintptr_t i = 0; i < fontManagement.loaded.Count(); i++) {
		// TODO Unmap file store data.
		Font font = fontManagement.loaded[i];

		if (font.type == FONT_TYPE_BITMAP) {
			// Nothing to be done.
		}

#ifdef USE_FREETYPE_AND_HARFBUZZ
		if (font.type == FONT_TYPE_FREETYPE_AND_HARFBUZZ) {
			hb_font_destroy(font.hb);
			FT_Done_Face(font.ft);
		}
#endif
	}

	for (uintptr_t i = 0; i < fontManagement.database.Length(); i++) {
		FontDatabaseEntry *entry = &fontManagement.database[i];

		for (uintptr_t j = 0; j < sizeof(entry->files) / sizeof(entry->files[0]); j++) {
			if (entry->files[j]) {
				FileStoreCloseHandle(entry->files[j]);
			}
		}
	}

	EsAssert(fontManagement.glyphCache.Count() == 0);
	EsAssert(fontManagement.glyphCacheBytes == 0);

	EsHeapFree(fontManagement.sansName);
	EsHeapFree(fontManagement.serifName);
	EsHeapFree(fontManagement.monospacedName);
	EsHeapFree(fontManagement.fallbackName);

	fontManagement.glyphCache.Free();
	fontManagement.substitutions.Free();
	fontManagement.database.Free();
	fontManagement.loaded.Free();

#ifdef USE_FREETYPE_AND_HARFBUZZ
	FT_Done_FreeType(fontManagement.freetypeLibrary);
#endif
}

// --------------------------------- Blitting rendered glyphs.

__attribute__((no_instrument_function))
inline static void DrawStringPixel(int oX, int oY, void *bitmap, size_t stride, uint32_t textColor, 
		uint32_t selectionColor, int32_t backgroundColor, uint32_t pixel, bool selected, bool fullAlpha) {
	uint32_t *destination = (uint32_t *) ((uint8_t *) bitmap + (oX) * 4 + (oY) * stride);
	uint8_t alpha = (textColor & 0xFF000000) >> 24;

	if (pixel == 0xFFFFFF && alpha == 0xFF) {
		*destination = 0xFF000000 | textColor;
	} else if (pixel && fullAlpha) {
		uint32_t original;

		if (selected) {
			original = selectionColor;
		} else if (backgroundColor < 0) {
			original = *destination;
		} else {
			original = backgroundColor;
		}

		uint32_t ga = (((pixel & 0x0000FF00) >> 8) * alpha) >> 8;
		uint32_t alphaD2 = (255 - ga) * ((original & 0xFF000000) >> 24);
		uint32_t alphaOut = ga + (alphaD2 >> 8);

		if (alphaOut) {
			uint32_t m2 = alphaD2 / alphaOut;
			uint32_t m1 = (ga << 8) / alphaOut;
			if (m2 == 0x100) m2--;
			if (m1 == 0x100) m1--;

			uint32_t r2 = m2 * ((original & 0x000000FF) >> 0);
			uint32_t g2 = m2 * ((original & 0x0000FF00) >> 8);
			uint32_t b2 = m2 * ((original & 0x00FF0000) >> 16);
			uint32_t r1 = m1 * ((textColor & 0x000000FF) >> 0);
			uint32_t g1 = m1 * ((textColor & 0x0000FF00) >> 8);
			uint32_t b1 = m1 * ((textColor & 0x00FF0000) >> 16);

			uint32_t result = 
				(0x00FF0000 & ((b1 + b2) << 8)) 
				| (0x0000FF00 & ((g1 + g2) << 0)) 
				| (0x000000FF & ((r1 + r2) >> 8))
				| (alphaOut << 24);

			*destination = result;
		}
	} else if (pixel) {
		uint32_t original;

		if (selected) {
			original = selectionColor;
		} else if (backgroundColor < 0) {
			original = *destination;
		} else {
			original = backgroundColor;
		}

		uint32_t ra = (((pixel & 0x000000FF) >> 0) * alpha) >> 8;
		uint32_t ga = (((pixel & 0x0000FF00) >> 8) * alpha) >> 8;
		uint32_t ba = (((pixel & 0x00FF0000) >> 16) * alpha) >> 8;
		uint32_t r2 = (255 - ra) * ((original & 0x000000FF) >> 0);
		uint32_t g2 = (255 - ga) * ((original & 0x0000FF00) >> 8);
		uint32_t b2 = (255 - ba) * ((original & 0x00FF0000) >> 16);
		uint32_t r1 = ra * ((textColor & 0x000000FF) >> 0);
		uint32_t g1 = ga * ((textColor & 0x0000FF00) >> 8);
		uint32_t b1 = ba * ((textColor & 0x00FF0000) >> 16);

		uint32_t result = 0xFF000000 | (0x00FF0000 & ((b1 + b2) << 8)) 
			| (0x0000FF00 & ((g1 + g2) << 0)) 
			| (0x000000FF & ((r1 + r2) >> 8));

		*destination = result;
	}
}

void DrawSingleCharacter(int width, int height, int xoff, int yoff, 
		EsPoint outputPosition, EsRectangle region, EsPaintTarget *target,  
		int blur, int type, bool selected, uint8_t *output,
		uint32_t color, uint32_t selectionColor, int32_t backgroundColor, bool fullAlpha) {
	// TODO Rewrite.

	if (type != CHARACTER_SUBPIXEL) {
		blur = 0;
	}

	uint8_t alpha = color >> 24;

	int xOut = outputPosition.x + xoff;
	int yOut = outputPosition.y + yoff;
	int xFrom = xOut, xTo = xOut + width;
	int yFrom = yOut, yTo = yOut + height;

	if (blur) {
		xFrom -= blur;
		yFrom -= blur;
		xTo   += blur;
		yTo   += blur;
	}

	if (xFrom < region.l) xFrom = region.l; else if (xFrom >= region.r) xFrom = region.r;
	if (xFrom < 0) xFrom = 0; else if (xFrom >= (int) target->width) xFrom = target->width;
	if (xTo < region.l) xTo = region.l; else if (xTo >= region.r) xTo = region.r;
	if (xTo < 0) xTo = 0; else if (xTo >= (int) target->width) xTo = target->width;

	if (yFrom < region.t) yFrom = region.t; else if (yFrom >= region.b) yFrom = region.b;
	if (yFrom < 0) yFrom = 0; else if (yFrom >= (int) target->height) yFrom = target->height;
	if (yTo < region.t) yTo = region.t; else if (yTo >= region.b) yTo = region.b;
	if (yTo < 0) yTo = 0; else if (yTo >= (int) target->height) yTo = target->height;

	float blurExponentDenominator = -1.0f / (2.0f * (blur / 3.0f) * (blur / 3.0f));

	for (int oY = yFrom; oY < yTo; oY++) {
		int y = oY - yOut;

		for (int oX = xFrom; oX < xTo; oX++) {
			int x = oX - xOut;

			if (blur) {
				float c = 0, d = 0;

				for (int i = y - blur; i <= y + blur; i++) {
					for (int j = x - blur; j <= x + blur; j++) {
						float weight = EsCRTexpf(blurExponentDenominator * ((i - y) * (i - y) + (j - x) * (j - x)));
						d += weight;

						if (i >= 0 && j >= 0 && i < height && j < width) {
							uint32_t pixel = *((uint32_t *) (output + (j * 4 + i * width * 4)));
							c += (pixel & 0xFF00) * weight;
						}
					}
				}

				uint32_t a = c / (d * 256.0f);
				DrawStringPixel(oX, oY, target->bits, target->stride, color, selectionColor, backgroundColor, 
						a | (a << 8) | (a << 16), selected, fullAlpha);
			} else if (type == CHARACTER_IMAGE || type == CHARACTER_RECOLOR) {
				uint32_t pixel = *((uint32_t *) (output + (x * 4 + y * width * 4)));
				uint32_t *destination = (uint32_t *) ((uint8_t *) target->bits + (oX) * 4 + (oY) * target->stride);

				if (type == CHARACTER_RECOLOR) {
					pixel = (pixel & 0xFF000000) | (color & 0x00FFFFFF);
				}

				if ((pixel >> 24) == 0xFF && alpha == 0xFF) {
					*destination = pixel;
				} else if (pixel && fullAlpha) {
					uint32_t original = *destination;
					uint32_t alphaSource = ((pixel >> 24) * alpha) >> 8;
					uint32_t alphaDestination = ((original & 0xFF000000) >> 24) * (255 - alphaSource);
					uint32_t alphaOut = alphaSource + (alphaDestination >> 8);

					if (alphaOut) {
						uint32_t m2 = alphaDestination / alphaOut;
						uint32_t m1 = (alphaSource << 8) / alphaOut;
						if (m2 == 0x100) m2--;
						if (m1 == 0x100) m1--;
						uint32_t r2 = m2 * ((original & 0x000000FF) >> 0);
						uint32_t g2 = m2 * ((original & 0x0000FF00) >> 8);
						uint32_t b2 = m2 * ((original & 0x00FF0000) >> 16);
						uint32_t r1 = m1 * ((pixel & 0x000000FF) >> 0);
						uint32_t g1 = m1 * ((pixel & 0x0000FF00) >> 8);
						uint32_t b1 = m1 * ((pixel & 0x00FF0000) >> 16);
						uint32_t result = (alphaOut << 24) | (0x00FF0000 & ((b1 + b2) << 8)) 
							| (0x0000FF00 & ((g1 + g2) << 0)) 
							| (0x000000FF & ((r1 + r2) >> 8));
						*destination = result;
					}
				} else if (pixel) {
					uint32_t original = *destination;
					uint32_t a = ((pixel >> 24) * alpha) >> 8;
					uint32_t r2 = (255 - a) * ((original & 0x000000FF) >> 0);
					uint32_t g2 = (255 - a) * ((original & 0x0000FF00) >> 8);
					uint32_t b2 = (255 - a) * ((original & 0x00FF0000) >> 16);
					uint32_t r1 = a * ((pixel & 0x000000FF) >> 0);
					uint32_t g1 = a * ((pixel & 0x0000FF00) >> 8);
					uint32_t b1 = a * ((pixel & 0x00FF0000) >> 16);
					uint32_t result = 0xFF000000 
						| (0x00FF0000 & ((b1 + b2) << 8)) 
						| (0x0000FF00 & ((g1 + g2) << 0)) 
						| (0x000000FF & ((r1 + r2) >> 8));
					*destination = result;
				}
			} else if (type == CHARACTER_SUBPIXEL) {
				uint32_t pixel = *((uint32_t *) (output + (x * 4 + y * width * 4)));
				DrawStringPixel(oX, oY, target->bits, target->stride, color, selectionColor, backgroundColor, pixel, selected, fullAlpha);
			}
		}
	}
}

// --------------------------------- Text shaping.

TextStyleDifference CompareTextStyles(const EsTextStyle *style1, const EsTextStyle *style2) {
	if (!style1) return TEXT_STYLE_NEW_FONT;
	if (style1->font.family 	!= style2->font.family) 	return TEXT_STYLE_NEW_FONT;
	if (style1->font.weight 	!= style2->font.weight) 	return TEXT_STYLE_NEW_FONT;
	if (style1->font.flags 		!= style2->font.flags) 		return TEXT_STYLE_NEW_FONT;
	if (style1->size 		!= style2->size) 		return TEXT_STYLE_NEW_FONT;
	if (style1->baselineOffset 	!= style2->baselineOffset) 	return TEXT_STYLE_NEW_SHAPE;
	if (style1->tracking 		!= style2->tracking) 		return TEXT_STYLE_NEW_SHAPE;
	if (style1->figures 		!= style2->figures) 		return TEXT_STYLE_NEW_SHAPE;
	if (style1->alternateDirection	!= style2->alternateDirection)	return TEXT_STYLE_NEW_SHAPE;
	if (style1->color		!= style2->color)		return TEXT_STYLE_NEW_RENDER;
	if (style1->blur		!= style2->blur)		return TEXT_STYLE_NEW_RENDER;
	if (style1->decorations		!= style2->decorations)		return TEXT_STYLE_NEW_RENDER;
	if (style1->decorationsColor	!= style2->decorationsColor)	return TEXT_STYLE_NEW_RENDER;
	return TEXT_STYLE_IDENTICAL;
}

#define TEXT_GET_CHARACTER_AT_POINT_MIDDLE  (1 << 0)
ptrdiff_t TextGetCharacterAtPoint(EsElement *element, const EsTextStyle *textStyle, const char *string, size_t stringBytes, int *_pointX, uint32_t flags) {
	// TODO Better integration with the EsTextPlan API.

	EsTextPlanProperties properties = {};
	EsTextRun textRuns[2] = {};
	textRuns[0].style = *textStyle;
	textRuns[1].offset = stringBytes;
	EsTextPlan *plan = EsTextPlanCreate(element, &properties, {}, string, textRuns, 1); 
	if (!plan) return 0;

	EsAssert(plan->lines.Length() == 1);
	bool useMiddle = flags & TEXT_GET_CHARACTER_AT_POINT_MIDDLE;
	int pointX = *_pointX;
	pointX *= FREETYPE_UNIT_SCALE;
	int currentX = 0, priorMiddle = 0;
	ptrdiff_t result = -1;

	for (uintptr_t j = 0; j < plan->lines[0].pieceCount; j++) {
		TextPiece *piece = &plan->pieces[plan->lines[0].pieceOffset + j];
		TextGlyphInfo *glyphs = &plan->glyphInfos[piece->glyphOffset];
		TextGlyphPosition *glyphPositions = &plan->glyphPositions[piece->glyphOffset];

		for (uintptr_t i = 0; i < piece->glyphCount; i++) {
			int left = useMiddle ? priorMiddle : currentX;
			int right = currentX + glyphPositions[i].x_advance / (useMiddle ? 2 : 1);

			priorMiddle = right;

			if (pointX >= left && pointX < right) {
				result = glyphs[i].cluster;
				goto done;
			}

			currentX += glyphPositions[i].x_advance;
		}
	}

	done:;
	*_pointX = currentX / FREETYPE_UNIT_SCALE;
	EsTextPlanDestroy(plan);
	return result;
}

int TextGetPartialStringWidth(EsElement *element, const EsTextStyle *textStyle, const char *fullString, size_t fullStringBytes, size_t measureBytes) {
	// TODO Better integration with the EsTextPlan API.

	EsTextPlanProperties properties = {};
	EsTextRun textRuns[3] = {};
	textRuns[0].style = *textStyle;
	textRuns[1].style = *textStyle;
	textRuns[1].offset = measureBytes;
	textRuns[2].offset = fullStringBytes;
	EsTextPlan *plan = EsTextPlanCreate(element, &properties, {}, fullString, textRuns, 2); 
	if (!plan) return 0;

	int width = 0;
	EsAssert(plan->lines.Length() == 1);

	for (uintptr_t i = 0; i < plan->lines[0].pieceCount; i++) {
		TextPiece *piece = &plan->pieces[plan->lines[0].pieceOffset + i];

		if (piece->start < measureBytes) {
			width += piece->width;
		}
	}

	EsTextPlanDestroy(plan);
	return width / FREETYPE_UNIT_SCALE;
}

int TextGetStringWidth(EsElement *element, const EsTextStyle *textStyle, const char *string, size_t stringBytes) {
	return TextGetPartialStringWidth(element, textStyle, string, stringBytes, stringBytes);
}

void TextTrimSpaces(EsTextPlan *plan) {
	if (~plan->properties.flags & ES_TEXT_PLAN_TRIM_SPACES) {
		return;
	}

	for (uintptr_t i = 0; i < plan->lines.Length(); i++) {
		TextLine *line = &plan->lines[i];

		if (!line->pieceCount) {
			continue;
		}

		TextPiece *firstPiece = &plan->pieces[line->pieceOffset];
		TextPiece *lastPiece = &plan->pieces[line->pieceOffset + line->pieceCount - 1];

		while (firstPiece->glyphCount && firstPiece->start != firstPiece->end
				&& plan->glyphInfos[firstPiece->glyphOffset].cluster == firstPiece->start 
				&& plan->string[firstPiece->start] == ' ') {
			line->width -= plan->glyphPositions[firstPiece->glyphOffset].x_advance;
			firstPiece->width -= plan->glyphPositions[firstPiece->glyphOffset].x_advance;
			firstPiece->glyphOffset++;
			firstPiece->glyphCount--;
			firstPiece->start++;
		}

		while (lastPiece->glyphCount && lastPiece->start != lastPiece->end
				&& plan->glyphInfos[lastPiece->glyphOffset + lastPiece->glyphCount - 1].cluster == lastPiece->end - 1
				&& plan->string[lastPiece->end - 1] == ' ') {
			line->width -= plan->glyphPositions[lastPiece->glyphOffset + lastPiece->glyphCount - 1].x_advance;
			lastPiece->width -= plan->glyphPositions[lastPiece->glyphOffset + lastPiece->glyphCount - 1].x_advance;
			lastPiece->glyphCount--;
			lastPiece->end--;
		}
	}
}

void TextPlaceEmergencyBreaks(EsTextPlan *plan, int32_t maximumLineWidth) {
	if ((plan->properties.flags & ES_TEXT_PLAN_CLIP_UNBREAKABLE_LINES) || maximumLineWidth == -1) {
		return;
	}

	repeat:;
	TextLine *line = &plan->lines.Last();
	if (line->width <= maximumLineWidth) return;
	EsAssert(line->pieceCount >= 1);

	int32_t x = 0, x0 = 0;
	uintptr_t j, piece;

	for (piece = 0; piece < line->pieceCount; piece++) {
		TextPiece *p = &plan->pieces[line->pieceOffset + piece];
		x0 = x;

		for (j = 0; j < p->glyphCount; j++) {
			int32_t width = plan->glyphPositions[p->glyphOffset + j].x_advance;

			if (x + width > maximumLineWidth && (j || piece)) {
				goto foundBreakPoint;
			}

			x += width;
		}
	}

	return; // One glyph on the line; we can't do anything.
	foundBreakPoint:;

	// Split the line.

	TextPiece *piece0 = &plan->pieces[line->pieceOffset + piece];
	TextPiece  piece1 = *piece0;
	piece1.width = piece0->width - (x - x0);
	piece1.glyphOffset += j;
	piece1.glyphCount = piece0->glyphCount - j;
	piece1.start = plan->glyphInfos[piece0->glyphOffset + j].cluster;
	piece0->end = piece1.start;
	piece0->width = x - x0;
	piece0->glyphCount = j;
	plan->pieces.Insert(piece1, line->pieceOffset + piece + 1);

	TextLine *line0 = line;
	TextLine  line1 = *line;
	line1.width -= x;
	line0->width = x;
	line1.pieceOffset += piece + 1;
	line1.pieceCount = line0->pieceCount - piece;
	line0->pieceCount = piece + 1;
	plan->lines.Add(line1);

	goto repeat;
}

void TextAddEllipsis(EsTextPlan *plan, int32_t maximumLineWidth, bool needFinalEllipsis, int32_t boundsWidth) {
	if (!boundsWidth || (~plan->properties.flags & ES_TEXT_ELLIPSIS)) {
		return;
	}

	bool needEllipsis = false;

	if (maximumLineWidth == -1) {
		for (uintptr_t i = 0; i < plan->lines.Length(); i++) {
			if (plan->lines[i].width > boundsWidth * FREETYPE_UNIT_SCALE) {
				maximumLineWidth = boundsWidth * FREETYPE_UNIT_SCALE;
				needEllipsis = true;
			}
		}
	} else {
		// Word-wrapping was enabled so lines won't exceed the boundary width.
	}

	if (!needEllipsis && !needFinalEllipsis) {
		return;
	}

	uint8_t ellipsisUTF8[3] = { 0xE2, 0x80, 0xA6 };

	// Shape and measure the ellipsis character.

	unsigned int glyphCount;
	TextGlyphInfo *glyphInfos;
	TextGlyphPosition *glyphPositions;
	FontShapeText(plan, (const char *) ellipsisUTF8, sizeof(ellipsisUTF8), 0, sizeof(ellipsisUTF8), nullptr, 0, &glyphCount, &glyphInfos, &glyphPositions);
	
	int32_t ellipsisWidth = 0;

	for (uintptr_t i = 0; i < glyphCount; i++) {
		ellipsisWidth += glyphPositions[i].x_advance;
	}

	for (uintptr_t i = needEllipsis ? 0 : plan->lines.Length() - 1; i < plan->lines.Length(); i++) {
		TextLine *line = &plan->lines[i];

		if (i == plan->lines.Length() - 1 && needFinalEllipsis) {
			// The maximum number of lines was exceeded, and this is the last permitted line, so add an ellipsis.
		} else if (line->width > boundsWidth * FREETYPE_UNIT_SCALE) {
			// This line exceeds the width boundary (and hence word-wrapping was disabled), so add an ellipsis.
		} else {
			continue;
		}
			
		// Make space for the ellipsis.

		int32_t spaceNeeded = ellipsisWidth - (maximumLineWidth - line->width);

		while (line->pieceCount && spaceNeeded > 0) {
			TextPiece *piece = &plan->pieces[line->pieceOffset + line->pieceCount - 1];

			if (piece->isTabPiece) {
				spaceNeeded -= piece->width;
				line->pieceCount--;
			} else if (piece->start == piece->end || !piece->glyphCount) {
				line->pieceCount--;
			} else {
				piece->end = plan->glyphInfos[piece->glyphOffset + piece->glyphCount - 1].cluster;
				int32_t width = plan->glyphPositions[piece->glyphOffset + piece->glyphCount - 1].x_advance;
				spaceNeeded -= width, line->width -= width, piece->width -= width;
				piece->glyphCount--;

				while (piece->glyphCount) {
					if (plan->glyphInfos[piece->glyphOffset + piece->glyphCount - 1].cluster == piece->end) {
						// TODO Test this branch!
						int32_t width = plan->glyphPositions[piece->glyphOffset + piece->glyphCount - 1].x_advance;
						spaceNeeded -= width, line->width -= width, piece->width -= width;
						piece->glyphCount--;
					} else {
						break;
					}
				}
			}
		}

		// Add the ellipsis.

		TextPiece piece = {};
		piece.style = plan->currentTextStyle;
		piece.glyphOffset = plan->glyphInfos.Length();
		piece.ascent  =  FontGetAscent (&plan->font) + plan->currentTextStyle->baselineOffset, 
		piece.descent = -FontGetDescent(&plan->font) - plan->currentTextStyle->baselineOffset;
		piece.isEllipsisPiece = true;

		for (uintptr_t i = 0; i < glyphCount; i++) {
			if (!plan->glyphInfos.Add(glyphInfos[i])) break;
			if (!plan->glyphPositions.Add(glyphPositions[i])) break;
			piece.glyphCount++;
			int32_t width = glyphPositions[i].x_advance;
			piece.width += width, line->width += width;
		}

		line->hasEllipsis = true;
		line->ellipsisPieceIndex = plan->pieces.Length();
		plan->pieces.Add(piece);
	}

	FontShapeTextDone(plan, glyphCount, glyphInfos, glyphPositions);
}

void TextItemizeByScript(EsTextPlan *plan, const EsTextRun *runs, size_t runCount, float sizeScaleFactor) {
	uint32_t lastAssignedScript = FALLBACK_SCRIPT;

	for (uintptr_t i = 0; i < runCount; i++) {
		uintptr_t offset = runs[i].offset;

		for (uintptr_t j = offset; j < runs[i + 1].offset;) {
			uint32_t codepoint = utf8_value(plan->string + j);
			uint32_t script;
			bool inheritingScript = false;

			if (codepoint == '\t') {
				// Tab characters should go in their own section.
				script = '\t';
			} else {
				script = FontGetScriptFromCodepoint(codepoint, &inheritingScript);
			}

			if (inheritingScript) {
				// TODO If this is a closing character, restore the last assigned script before the most recent opening character.
				script = lastAssignedScript == '\t' ? FALLBACK_SCRIPT : lastAssignedScript;
			}

			if (lastAssignedScript != script && j != runs[i].offset) {
				TextRun run = {};
				run.style = runs[i].style;
				run.offset = offset;
				run.script = lastAssignedScript;
				run.style.font.family = FontApplySubstitution(&plan->properties, run.style.font.family, run.script);
				run.style.size *= sizeScaleFactor;
				plan->textRuns.Add(run);
				offset = j;
			}

			lastAssignedScript = script;
			j = utf8_advance(plan->string + j) - plan->string;
		}

		TextRun run = {};
		run.style = runs[i].style;
		run.offset = offset;
		run.script = lastAssignedScript;
		run.style.font.family = FontApplySubstitution(&plan->properties, run.style.font.family, run.script);
		run.style.size *= sizeScaleFactor;
		plan->textRuns.Add(run);
	}

	TextRun run = {};
	run.offset = runs[runCount].offset;
	plan->textRuns.Add(run);
}

void TextUpdateFont(EsTextPlan *plan, const EsTextStyle *style) {
	if (TEXT_STYLE_NEW_FONT == CompareTextStyles(plan->currentTextStyle, style)) {
		plan->font = FontGet(style->font);
		FontSetSize(&plan->font, style->size * FREETYPE_UNIT_SCALE);
	}

	plan->currentTextStyle = style;
}

int32_t TextExpandTabs(EsTextPlan *plan, uintptr_t pieceOffset, int32_t width) {
	int32_t addedWidth = 0;

	for (uintptr_t i = pieceOffset; i < plan->pieces.Length(); i++) {
		TextPiece *piece = &plan->pieces[i];

		if (piece->isTabPiece) {
			TextUpdateFont(plan, piece->style);
			int32_t emWidth = FontGetEmWidth(&plan->font) * FREETYPE_UNIT_SCALE;
			int32_t tabWidth = emWidth * 4; // TODO Make spaces-per-tab customizable.
			int32_t firstWidth = emWidth + tabWidth - (width + emWidth) % tabWidth;
			piece->width = firstWidth + tabWidth * (piece->end - piece->start - 1);
			addedWidth += piece->width;
			piece->glyphOffset = plan->glyphInfos.Length();
			piece->glyphCount = 0;
			piece->ascent = FontGetAscent(&plan->font);
			piece->descent = -FontGetDescent(&plan->font);

			for (uintptr_t i = 0; i < piece->end - piece->start; i++) {
				TextGlyphInfo info = {};
				info.cluster = piece->start + i;
				info.codepoint = 0xFFFFFFFF;
				TextGlyphPosition position = {};
				position.x_advance = i ? tabWidth : firstWidth;
				if (!plan->glyphInfos.Add(info)) break;
				if (!plan->glyphPositions.Add(position)) break;
				piece->glyphCount++;
			}
		}

		width += piece->width;
	}

	return addedWidth;
}

int32_t TextBuildTextPieces(EsTextPlan *plan, uintptr_t sectionStart, uintptr_t sectionEnd) {
	// Find the first run that contains the section.

	for (; plan->textRunPosition < plan->textRuns.Length() - 1; plan->textRunPosition++) {
		if (plan->textRuns[plan->textRunPosition].offset <= sectionStart && plan->textRuns[plan->textRunPosition + 1].offset > sectionStart) {
			break;
		}
	}

	EsAssert(plan->textRunPosition != plan->textRuns.Length() - 1);

	// Iterate through each run in the section.

	int32_t width = 0;

	while (plan->textRunPosition != plan->textRuns.Length() - 1) {
		TextRun *run = &plan->textRuns[plan->textRunPosition];

		uintptr_t start = sectionStart > run[0].offset ? sectionStart : run[0].offset;
		uintptr_t end   = sectionEnd   < run[1].offset ? sectionEnd   : run[1].offset;

		if (end <= start) {
			break;
		}

		// Update the font to match the run.

		TextUpdateFont(plan, &run->style);

		// Don't shape newline characters.

		while (start < end && plan->string[start] == '\n') start++;
		while (end - 1 > start && plan->string[end - 1] == '\n') end--;

		if (end == start) {
			plan->textRunPosition++;
			continue;
		}

		EsAssert(end > start);

		// Handle tab characters specially.

		if (plan->string[start] == '\t') {
			TextPiece _piece = {};

			if (plan->pieces.Add(_piece)) {
				TextPiece *piece = &plan->pieces.Last();
				piece->style = plan->currentTextStyle;
				piece->glyphOffset = 0;
				piece->glyphCount = 0;
				piece->start = start;
				piece->end = end;
				piece->isTabPiece = true;
			}

			plan->textRunPosition++;
			continue;
		}

		// Shape the run.

		TextFeature features[4] = {};
		size_t featureCount = 0;

#ifdef USE_FREETYPE_AND_HARFBUZZ
		if (plan->currentTextStyle->figures == ES_TEXT_FIGURE_OLD) hb_feature_from_string("onum", -1, features + (featureCount++));
		if (plan->currentTextStyle->figures == ES_TEXT_FIGURE_TABULAR) hb_feature_from_string("tnum", -1, features + (featureCount++));
#endif
		plan->segmentProperties.script = (TextScript) run->script;

		unsigned int glyphCount;
		TextGlyphInfo *glyphInfos;
		TextGlyphPosition *glyphPositions;
		FontShapeText(plan, plan->string, plan->breaker.bytes, start, end - start, features, featureCount, &glyphCount, &glyphInfos, &glyphPositions);

		// Create the text piece.

		TextPiece _piece = {};

		if (plan->pieces.Add(_piece)) {
			TextPiece *piece = &plan->pieces.Last();
			piece->style = plan->currentTextStyle;
			piece->glyphOffset = plan->glyphInfos.Length();
			piece->glyphCount = 0;
			piece->ascent  =  FontGetAscent (&plan->font) + plan->currentTextStyle->baselineOffset;
			piece->descent = -FontGetDescent(&plan->font) - plan->currentTextStyle->baselineOffset;
			piece->start = start;
			piece->end = end;

			for (uintptr_t i = 0; i < glyphCount; i++) {
				if (!plan->glyphInfos.Add(glyphInfos[i])) break;
				if (!plan->glyphPositions.Add(glyphPositions[i])) break;
				piece->glyphCount++;

				piece->width += glyphPositions[i].x_advance;

				if (i == glyphCount - 1 || glyphInfos[i].cluster != glyphInfos[i + 1].cluster) {
					piece->width += plan->currentTextStyle->tracking * FREETYPE_UNIT_SCALE;
				}

				// EsPrint("\t%d\n", glyphInfos[i].codepoint);
			}

			width += piece->width;
		}

		// Go to the next run.

		plan->textRunPosition++;

		FontShapeTextDone(plan, glyphCount, glyphInfos, glyphPositions);
	}

	plan->textRunPosition--;

	return width;
}

void TextPlanDestroy(EsTextPlan *plan) {
	plan->glyphInfos.Free();
	plan->glyphPositions.Free();
	plan->pieces.Free();
	plan->lines.Free();
	plan->textRuns.Free();
}

EsTextPlan *EsTextPlanCreate(EsElement *element, EsTextPlanProperties *properties, EsRectangle bounds, const char *string, const EsTextRun *formatRuns, size_t formatRunCount) {
	// TODO Bidirectional text (UAX9). 
	// TODO Vertical text layout (UAX50).
	// TODO Supporting arbitrary OpenType features.
	// TODO Reshaping lines once word wrapping is applied.
	// TODO Don't break lines in the middle of emails/URLs.

	// EsPrint("EsTextPlanCreate... width %d\n", Width(bounds) * FREETYPE_UNIT_SCALE);

	EsMessageMutexCheck();
#if 0
	EsAssert(element);
#endif
	float scale = theming.scale; // TODO Get the scale factor from the element's window.

	EsTextPlan plan = {};

	// Initialise the line breaker.

	plan.breaker.string = string;
	plan.breaker.bytes = formatRuns[formatRunCount].offset;
	EsAssert(plan.breaker.bytes < 0x80000000);

	if (!plan.breaker.bytes) {
		return nullptr; // Empty input.
	}

	// Initialise the plan.

	plan.string = string;
	plan.singleUse = properties->flags & ES_TEXT_PLAN_SINGLE_USE;
	plan.properties = *properties;

	TextLine blankLine = {};

	if (!plan.lines.Add(blankLine)) {
		return nullptr;
	}

	// Setup the HarfBuzz buffer.

	FontInitialiseShaping(&plan);

	// Subdivide the runs by character script.
	// This is also responsible for scaling the text sizes.
	
	TextItemizeByScript(&plan, formatRuns, formatRunCount, scale);

	// Layout the paragraph.

	int32_t maximumLineWidth = Width(bounds) && (properties->flags & ES_TEXT_WRAP) ? Width(bounds) * FREETYPE_UNIT_SCALE : -1;
	Break previousBreak = {};
	bool needEllipsis = false;

	while (previousBreak.position != plan.breaker.bytes) {
		// Find the next break opportunity.

		Break nextBreak = plan.breaker.Next();

		while (!plan.breaker.error && !nextBreak.forced && (~properties->flags & ES_TEXT_WRAP)) {
			nextBreak = plan.breaker.Next();
		}

		if (plan.breaker.error) {
			break;
		}

		// Build the text pieces for this section.

		uintptr_t pieceOffset = plan.pieces.Length();
		int32_t width = TextBuildTextPieces(&plan, previousBreak.position, nextBreak.position);
		width += TextExpandTabs(&plan, pieceOffset, plan.lines.Last().width);

		// Should we start a new line?

		if (previousBreak.forced || (maximumLineWidth != -1 && plan.lines.Last().width + width > maximumLineWidth)) {
			if (properties->maxLines == (int32_t) plan.lines.Length()) {
				needEllipsis = true;
				break;
			}

			if (plan.lines.Add(blankLine)) {
				plan.lines.Last().pieceOffset = pieceOffset;
			}
		}

#if 0
		EsPrint("\tadded section '%s' to line %d (%d pieces) at x=%d\n", 
				nextBreak.position - previousBreak.position, string + previousBreak.position, 
				ArrayLength(plan.lines) - 1, ArrayLength(plan.pieces) - pieceOffset, 
				plan.lines.Last().width);
#endif

		// Add the pieces to the line.

		TextLine *line = &plan.lines.Last();
		TextExpandTabs(&plan, pieceOffset, line->width);

		for (uintptr_t i = pieceOffset; i < plan.pieces.Length(); i++) {
			line->width += plan.pieces[i].width;
			line->pieceCount++;
		}

		TextPlaceEmergencyBreaks(&plan, maximumLineWidth);

		// Go to the next section.

		previousBreak = nextBreak;
	}

	// Calculate the ascent/descent of each line.

	for (uintptr_t i = 0; i < plan.lines.Length(); i++) {
		TextLine *line = &plan.lines[i];

		if (!line->pieceCount && i) {
			// If the line doesn't have any pieces, it must be from a double newline.
			// Inherit the ascent/descent of the previous line.

			line->ascent = line[-1].ascent;
			line->descent = line[-1].descent;
		}

		for (uintptr_t i = line->pieceOffset; i < line->pieceOffset + line->pieceCount; i++) {
			if (line->ascent < plan.pieces[i].ascent) line->ascent = plan.pieces[i].ascent;
			if (line->descent < plan.pieces[i].descent) line->descent = plan.pieces[i].descent;
		}
	}

	// Trim leading and trailing spaces.
	
	TextTrimSpaces(&plan);

	// Add a terminating ellipsis.

	TextAddEllipsis(&plan, maximumLineWidth, needEllipsis, Width(bounds));

	// Calculate the total width and height.

	for (uintptr_t i = 0; i < plan.lines.Length(); i++) {
		plan.totalHeight += plan.lines[i].ascent + plan.lines[i].descent;

		if (plan.lines[i].width > plan.totalWidth) {
			plan.totalWidth = plan.lines[i].width;
		}
	}

	// Destroy the HarfBuzz buffer.

	FontDestroyShaping(&plan);
	
	// Return the plan.

	EsTextPlan *copy = (EsTextPlan *) EsHeapAllocate(sizeof(EsTextPlan), true);
	
	if (copy) {
		*copy = plan;
		return copy;
	} else {
		TextPlanDestroy(&plan);
		return nullptr;
	}
}

void EsTextPlanDestroy(EsTextPlan *plan) {
	EsMessageMutexCheck();
	EsAssert(!plan->singleUse);
	TextPlanDestroy(plan);
	EsHeapFree(plan);
}

void EsTextPlanReplaceStyleRenderProperties(EsTextPlan *plan, const EsTextStyle *style) {
	for (uintptr_t i = 0; i < plan->textRuns.Length() - 1; i++) {
		plan->textRuns[i].style.color = style->color;
		plan->textRuns[i].style.blur = style->blur;
		plan->textRuns[i].style.decorations = style->decorations;
		plan->textRuns[i].style.decorationsColor = style->decorationsColor;
	}
}

int EsTextPlanGetWidth(EsTextPlan *plan) {
	return (plan->totalWidth + FREETYPE_UNIT_SCALE - 1) / FREETYPE_UNIT_SCALE;
}

int EsTextPlanGetHeight(EsTextPlan *plan) {
	return (plan->totalHeight + FREETYPE_UNIT_SCALE - 1) / FREETYPE_UNIT_SCALE;
}

size_t EsTextPlanGetLineCount(EsTextPlan *plan) {
	return plan->lines.Length();
}

EsTextStyle TextPlanGetPrimaryStyle(EsTextPlan *plan) {
	return plan->textRuns[0].style;
}

void DrawTextPiece(EsPainter *painter, EsTextPlan *plan, TextPiece *piece, TextLine *line,
		int32_t cursorX, int32_t cursorY, 
		const EsTextSelection *selection, uintptr_t caret, int32_t selectionBackgroundBottom) {
	if (cursorX / FREETYPE_UNIT_SCALE > painter->clip.r 
			|| (cursorX + piece->width) / FREETYPE_UNIT_SCALE < painter->clip.l
			|| cursorY > painter->clip.b
			|| (cursorY + (piece->ascent + piece->descent) / FREETYPE_UNIT_SCALE) < painter->clip.t) {
		return;
	}

#if 0
	EsPrint("\tdrawing piece; '%s' on line %d glyphOffset %d and glyphCount %d at %i, %i with caret %d\n", 
			piece->end - piece->start, plan->string + piece->start,
			line - plan->lines, piece->glyphOffset, piece->glyphCount, 
			cursorX / FREETYPE_UNIT_SCALE, cursorY, caret);
#endif

	// Prevent issues with negative numbers getting rounded differently...
	int32_t cursorXIntegerOffset = -(0x40000000 / FREETYPE_UNIT_SCALE);
	cursorX += 0x40000000;
	int32_t cursorXStart = cursorX;

	TextGlyphInfo *glyphs = &plan->glyphInfos[piece->glyphOffset];
	TextGlyphPosition *glyphPositions = &plan->glyphPositions[piece->glyphOffset];

	// Update the font to match the piece.

	TextUpdateFont(plan, piece->style);

	// Draw the selection background.

	if (selection->caret0 != selection->caret1 && !selection->hideCaret && !piece->isEllipsisPiece) {
		int sCursorX = cursorX, selectionStartX = -1, selectionEndX = -1;

		for (uintptr_t i = 0; i < piece->glyphCount; i++) {
			if (selectionStartX == -1 && (int32_t) glyphs[i].cluster >= selection->caret0) {
				selectionStartX = sCursorX;
			}

			if (selectionEndX == -1 && (int32_t) glyphs[i].cluster >= selection->caret1) {
				selectionEndX = sCursorX;
			}

			sCursorX += glyphPositions[i].x_advance;

			if (i == piece->glyphCount - 1 || glyphs[i].cluster != glyphs[i + 1].cluster) {
				sCursorX += plan->currentTextStyle->tracking;
			}
		}

		if (selectionStartX == -1 && selection->caret0 >= 0) {
			selectionStartX = sCursorX;
		}

		if (selectionEndX == -1) {
			selectionEndX = sCursorX;
		}

		EsRectangle s;
		s.l = (selectionStartX + FREETYPE_UNIT_SCALE / 2) / FREETYPE_UNIT_SCALE + cursorXIntegerOffset;
		s.t = cursorY;
		s.r = (selectionEndX + FREETYPE_UNIT_SCALE / 2) / FREETYPE_UNIT_SCALE + cursorXIntegerOffset;
		s.b = selectionBackgroundBottom;
		EsDrawBlock(painter, s, selection->background);
	}

	// Draw each glyph in the piece.

	int32_t caretX = -1, caretY = cursorY;

	for (uintptr_t i = 0; i < piece->glyphCount; i++) {
		uint32_t codepoint = glyphs[i].codepoint;

		int positionX = (glyphPositions[i].x_offset + cursorX) / FREETYPE_UNIT_SCALE + cursorXIntegerOffset, 
		    positionY = ((glyphPositions[i].y_offset + FREETYPE_UNIT_SCALE / 2) / FREETYPE_UNIT_SCALE + cursorY);
		uint32_t color = plan->currentTextStyle->color;

		GlyphCacheKey key = {};
		key.glyphIndex = codepoint;
		key.size = plan->currentTextStyle->size * FREETYPE_UNIT_SCALE;
		key.font = plan->font;
		GlyphCacheEntry *entry = nullptr;

		if (codepoint == 0xFFFFFFFF || key.size > 2000) {
			goto nextCharacter;
		}

		if (key.size > 25) {
			key.fractionalPosition = 0;
		} else if (key.size > 15) {
			key.fractionalPosition = ((glyphPositions[i].x_offset + cursorX) & 0x3F) & 0x20;
		} else {
			key.fractionalPosition = ((glyphPositions[i].x_offset + cursorX) & 0x3F) & 0x30;
		}

		entry = LookupGlyphCacheEntry(key);

		if (!entry) {
			goto nextCharacter;
		}

		if (!entry->data) {
			// EsPrint("Rendering '%c' in size %d\n", plan->string[glyphs[i].cluster], key.size);
			FontSetSize(&key.font, key.size);

			if (!FontRenderGlyph(key, entry)) {
				EsHeapFree(entry);
				goto nextCharacter;
			} else {
				RegisterGlyphCacheEntry(key, entry);
			}
		}

		EsAssert(0 == EsMemoryCompare(&entry->key, &key, sizeof(GlyphCacheKey)));

		if (selection->caret0 != selection->caret1 && !selection->hideCaret 
				&& (int32_t) glyphs[i].cluster >= selection->caret0 && (int32_t) glyphs[i].cluster < selection->caret1
				&& selection->foreground) {
			color = selection->foreground;
		}

		DrawSingleCharacter(entry->width, entry->height, entry->xoff, entry->yoff, 
				ES_POINT(positionX, positionY + line->ascent / FREETYPE_UNIT_SCALE), 
				painter->clip, painter->target, 
				plan->currentTextStyle->blur, 
				CHARACTER_SUBPIXEL, 
				false, entry->data, 
				color, 0, -1, painter->target->fullAlpha);

		nextCharacter:;

		if (caretX == -1 && glyphs[i].cluster >= caret) {
			caretX = (cursorX + FREETYPE_UNIT_SCALE / 2) / FREETYPE_UNIT_SCALE;
		}

		cursorX += glyphPositions[i].x_advance;
		cursorY += (glyphPositions[i].y_advance + FREETYPE_UNIT_SCALE / 2) / FREETYPE_UNIT_SCALE;

		if (i == piece->glyphCount - 1 || glyphs[i].cluster != glyphs[i + 1].cluster) {
			cursorX += plan->currentTextStyle->tracking * FREETYPE_UNIT_SCALE;
		}
	}

	if (caretX == -1) {
		caretX = (cursorX + FREETYPE_UNIT_SCALE / 2) / FREETYPE_UNIT_SCALE;
	}

	// Draw the caret.

	if (!selection->hideCaret && caret >= piece->start 
			&& (caret < piece->end || (caret == piece->end && piece == &plan->pieces[line->pieceOffset + line->pieceCount - 1]))) {
		caretX += cursorXIntegerOffset;

#if 0
		if (selection->snapCaretToInsets && selection->caret0 == selection->caret1) {
			EsRectangle insets = EsPainterBoundsInset(painter); 
			// EsPrint("%d, %d, %d\n", caretX + bounds.l, insets.l, insets.r);

			if (caretX >= insets.l - 1 && caretX <= insets.l + 1) {
				caretX = insets.l;
			} else if (caretX >= insets.r - 2 && caretX <= insets.r) {
				caretX = insets.r - 1;
			}
		}
#endif

		int caretWidth = theming.scale; // TODO Make this a system constant.
		EsDrawInvert(painter, ES_RECT_4(caretX, caretX + caretWidth, caretY, selectionBackgroundBottom));
	}

	// Draw decorations.

	{
		int32_t thickness = piece->style->size / 15 + 1;
		uint32_t color = piece->style->decorationsColor ?: piece->style->color;

		EsRectangle bounds;
		bounds.l = cursorXStart / FREETYPE_UNIT_SCALE + cursorXIntegerOffset;
		bounds.r = cursorX / FREETYPE_UNIT_SCALE + cursorXIntegerOffset;

		if (piece->style->decorations & ES_TEXT_DECORATION_STRIKE_THROUGH) {
			int32_t center = cursorY + (line->ascent + line->descent) / FREETYPE_UNIT_SCALE / 2;
			bounds.t = center - thickness / 2 + 1;
			bounds.b = center + (thickness + 1) / 2 + 1;
			EsDrawBlock(painter, bounds, color);
		}

		if (piece->style->decorations & ES_TEXT_DECORATION_UNDERLINE) {
			int32_t baseline = cursorY + line->ascent / FREETYPE_UNIT_SCALE;
			bounds.t = baseline + thickness;
			bounds.b = baseline + thickness * 2;
			EsDrawBlock(painter, bounds, color);
		}
	}
}

void EsDrawText(EsPainter *painter, EsTextPlan *plan, EsRectangle bounds, const EsRectangle *_clip, const EsTextSelection *_selection) {
	EsMessageMutexCheck();

	if (!plan) return;

	// EsPrint("EsDrawText... '%s' in %R\n", plan->textRuns[plan->textRunCount].offset, plan->string, bounds);

	// TODO Underlined text.
	// TODO Inline images and icons.

	// Work out the selection we should display.

	EsTextSelection selection = {};
	if (_selection) selection = *_selection;
	uintptr_t caret = selection.caret1;

	if (selection.caret0 > selection.caret1) {
		int swap = selection.caret1;
		selection.caret1 = selection.caret0;
		selection.caret0 = swap;
	} else if (!_selection) {
		selection.hideCaret = true;
	}

	// Calculate the area we're drawing into.

	int32_t maximumLineWidth = Width(bounds), maximumHeight = Height(bounds);
	EsRectangle oldClip = painter->clip;
	if (_clip) EsRectangleClip(*_clip, painter->clip, &painter->clip);
	int32_t cursorY = (plan->properties.flags & ES_TEXT_V_CENTER) ? (maximumHeight - EsTextPlanGetHeight(plan)) / 2 
		: (plan->properties.flags & ES_TEXT_V_BOTTOM) ? maximumHeight - EsTextPlanGetHeight(plan) : 0;

	// Iterate through each line.

	for (uintptr_t i = 0; i < plan->lines.Length(); i++) {
		TextLine *line = &plan->lines[i];

		int32_t cursorX = (plan->properties.flags & ES_TEXT_H_CENTER) ? ((maximumLineWidth * FREETYPE_UNIT_SCALE - line->width) / 2)
			: (plan->properties.flags & ES_TEXT_H_RIGHT) ? (maximumLineWidth * FREETYPE_UNIT_SCALE - line->width) : 0;

		int32_t selectionBackgroundBottom;
		
		if (plan->lines.Length() == 1 && (plan->properties.flags & ES_TEXT_V_CENTER)) {
			// If this is a single, centered line, make sure that the selection background bottom edge
			// is the same distance from the destination bounds as it is for the top edge.
			selectionBackgroundBottom = bounds.b - cursorY;
		} else {
			selectionBackgroundBottom = cursorY + bounds.t + (line->ascent + line->descent + FREETYPE_UNIT_SCALE / 2) / FREETYPE_UNIT_SCALE;
		}

		// Draw each text piece on the line.

		for (uintptr_t j = 0; j < line->pieceCount; j++) {
			TextPiece *piece = &plan->pieces[line->pieceOffset + j];
			DrawTextPiece(painter, plan, piece, line, cursorX + bounds.l * FREETYPE_UNIT_SCALE, cursorY + bounds.t, &selection, caret, selectionBackgroundBottom);
			cursorX += piece->width;
		}

		if (line->hasEllipsis) {
			TextPiece *piece = &plan->pieces[line->ellipsisPieceIndex];
			DrawTextPiece(painter, plan, piece, line, cursorX + bounds.l * FREETYPE_UNIT_SCALE, cursorY + bounds.t, &selection, caret, selectionBackgroundBottom);
			cursorX += piece->width;
		}

		cursorY += (line->ascent + line->descent + FREETYPE_UNIT_SCALE / 2) / FREETYPE_UNIT_SCALE;
	}

	// Destroy the plan if it is single use.

	if (plan->singleUse) {
		plan->singleUse = false;
		EsTextPlanDestroy(plan);
	}

	painter->clip = oldClip;
}

void EsDrawTextSimple(EsPainter *painter, EsElement *element, EsRectangle bounds, const char *string, ptrdiff_t stringBytes, EsTextStyle style, uint32_t flags) {
	EsTextPlanProperties properties = {};
	properties.flags = ES_TEXT_PLAN_SINGLE_USE | flags;
	EsTextRun textRuns[2] = {};
	textRuns[0].style = style;
	textRuns[1].offset = stringBytes == -1 ? EsCStringLength(string) : stringBytes;
	if (!textRuns[1].offset) return;
	EsDrawText(painter, EsTextPlanCreate(element, &properties, bounds, string, textRuns, 1), bounds); 
}
