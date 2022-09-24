#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

#define ES_FONT_ITALIC (1 << 0)
#define ES_TEXT_FIGURE_DEFAULT (0)
#define ES_TEXT_FIGURE_OLD (1)
#define ES_TEXT_FIGURE_TABULAR (2)
#define ES_TEXT_DECORATION_UNDERLINE (1 << 0)
#define ES_TEXT_DECORATION_STRIKE_THROUGH (1 << 1)
#define ES_TEXT_H_LEFT (1 << 0)
#define ES_TEXT_H_CENTER (1 << 1)
#define ES_TEXT_H_RIGHT (1 << 2)
#define ES_TEXT_V_TOP (1 << 3)
#define ES_TEXT_V_CENTER (1 << 4)
#define ES_TEXT_V_BOTTOM (1 << 5)
#define ES_TEXT_ELLIPSIS (1 << 6)
#define ES_TEXT_WRAP (1 << 7)
#define ES_TEXT_PLAN_SINGLE_USE (1 << 8)
#define ES_TEXT_PLAN_TRIM_SPACES (1 << 9)
#define ES_TEXT_PLAN_RTL (1 << 10)
#define ES_TEXT_PLAN_CLIP_UNBREAKABLE_LINES (1 << 11)
#define ES_TEXT_PLAN_NO_FONT_SUBSTITUTION (1 << 12)
#define ES_FONT_SANS       (0xFFFF)
#define ES_FONT_SERIF      (0xFFFE)
#define ES_FONT_MONOSPACED (0xFFFD)
#define ES_FONT_REGULAR  (4)
#define ES_FONT_SEMIBOLD (6)
#define ES_FONT_BOLD     (7)

typedef uint16_t EsFontFamily;
typedef const char *EsCString;

struct EsElement;
struct EsTextPlan;

typedef struct EsRectangle {
	int32_t l, r, t, b;
} EsRectangle;

typedef struct EsPaintTarget {
	void *bits;
	uint32_t width, height, stride;
	bool fullAlpha, readOnly, fromBitmap, forWindowManager;
} EsPaintTarget;

typedef struct EsPainter {
	EsRectangle clip;
	int32_t offsetX, offsetY, width, height;
	void *style;
	EsPaintTarget *target;
} EsPainter;

typedef struct EsFont {
	EsFontFamily family;
	uint8_t weight;
	uint8_t flags;
} EsFont;

typedef struct EsFontInformation {
	char name[96];
	size_t nameBytes;
	char category[32];
	size_t categoryBytes;
	EsFontFamily id;
	uint16_t availableWeightsNormal;
	uint16_t availableWeightsItalic;
} EsFontInformation;

typedef struct EsTextStyle {
	// TODO Indicating RTL/vertical writing modes.
	// TODO Support for features.
	// TODO Support for variable fonts.

	EsFont font;
	float size;

	uint16_t baselineOffset;
	int8_t tracking;
	uint8_t figures;
	bool alternateDirection;

	// Render properties:
	uint8_t blur;
	uint8_t decorations;
	uint32_t color;
	uint32_t decorationsColor;
} EsTextStyle;

typedef struct EsTextPlanProperties {
	EsCString cLanguage;
	uint32_t flags; 
	int maxLines; // Set to 0 for no limit.
} EsTextPlanProperties;

typedef struct EsTextRun {
	EsTextStyle style;
	uint32_t offset;
} EsTextRun;

typedef struct EsTextSelection {
	ptrdiff_t caret0, caret1;
	bool hideCaret;
	uint32_t foreground, background;
} EsTextSelection;

#ifdef __cplusplus
#define ES_EXTERN_C extern "C"
#else
#define ES_EXTERN_C
#endif

ES_EXTERN_C struct EsTextPlan *EsTextPlanCreate(struct EsElement *element, EsTextPlanProperties *properties, EsRectangle bounds, const char *string, const EsTextRun *textRuns, size_t textRunCount);
ES_EXTERN_C int EsTextPlanGetWidth(struct EsTextPlan *plan);
ES_EXTERN_C int EsTextPlanGetHeight(struct EsTextPlan *plan);
ES_EXTERN_C size_t EsTextPlanGetLineCount(struct EsTextPlan *plan);
ES_EXTERN_C void EsTextPlanDestroy(struct EsTextPlan *plan);
ES_EXTERN_C void EsTextPlanReplaceStyleRenderProperties(struct EsTextPlan *plan, const EsTextStyle *style);

#ifdef __cplusplus
ES_EXTERN_C void EsDrawText(EsPainter *painter, struct EsTextPlan *plan, EsRectangle bounds, const EsRectangle *clip = nullptr, const EsTextSelection *selectionProperties = nullptr);
#else
ES_EXTERN_C void EsDrawText(EsPainter *painter, struct EsTextPlan *plan, EsRectangle bounds, const EsRectangle *clip, const EsTextSelection *selectionProperties);
#endif
ES_EXTERN_C void EsDrawTextSimple(EsPainter *painter, struct EsElement *element, EsRectangle bounds, const char *string, ptrdiff_t stringBytes, EsTextStyle style, uint32_t flags);

ES_EXTERN_C bool EsFontDatabaseLookupByName(const char *name, ptrdiff_t nameBytes, EsFontInformation *information);
ES_EXTERN_C bool EsFontDatabaseLookupByID(EsFontFamily id, EsFontInformation *information);

ES_EXTERN_C void TextAddFont(const char *cName, const char *cCategory, const char *cScripts, const char *const *const cPaths, const char *const *const cTypes, size_t fileCount);

#define ES_RECT_1(x) ((EsRectangle) { (int32_t) (x), (int32_t) (x), (int32_t) (x), (int32_t) (x) })
#define ES_RECT_1I(x) ((EsRectangle) { (int32_t) (x), (int32_t) -(x), (int32_t) (x), (int32_t) -(x) })
#define ES_RECT_1S(x) ((EsRectangle) { 0, (int32_t) (x), 0, (int32_t) (x) })
#define ES_RECT_2(x, y) ((EsRectangle) { (int32_t) (x), (int32_t) (x), (int32_t) (y), (int32_t) (y) })
#define ES_RECT_2I(x, y) ((EsRectangle) { (int32_t) (x), (int32_t) -(x), (int32_t) (y), (int32_t) -(y) })
#define ES_RECT_2S(x, y) ((EsRectangle) { 0, (int32_t) (x), 0, (int32_t) (y) })
#define ES_RECT_4(x, y, z, w) ((EsRectangle) { (int32_t) (x), (int32_t) (y), (int32_t) (z), (int32_t) (w) })
#define ES_RECT_4PD(x, y, w, h) ((EsRectangle) { (int32_t) (x), (int32_t) ((x) + (w)), (int32_t) (y), (int32_t) ((y) + (h)) })
#define ES_RECT_WIDTH(_r) ((_r).r - (_r).l)
#define ES_RECT_HEIGHT(_r) ((_r).b - (_r).t)
#define ES_RECT_TOTAL_H(_r) ((_r).r + (_r).l)
#define ES_RECT_TOTAL_V(_r) ((_r).b + (_r).t)
#define ES_RECT_SIZE(_r) ES_RECT_WIDTH(_r), ES_RECT_HEIGHT(_r)
#define ES_RECT_TOP_LEFT(_r) (_r).l, (_r).t
#define ES_RECT_BOTTOM_LEFT(_r) (_r).l, (_r).b
#define ES_RECT_BOTTOM_RIGHT(_r) (_r).r, (_r).b
#define ES_RECT_ALL(_r) (_r).l, (_r).r, (_r).t, (_r).b
#define ES_RECT_VALID(_r) (ES_RECT_WIDTH(_r) > 0 && ES_RECT_HEIGHT(_r) > 0)
#define ES_POINT(x, y) ((EsPoint) { (int32_t) (x), (int32_t) (y) })
