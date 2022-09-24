// TODO [bin/render_vector res/FigmaTestDump.ini "Group 5" bin/test.png] is generating slightly transparent pixels around the text!?
// TODO Text sizing is inaccurate.
// TODO Shadow sizing/opacity is inaccurate.

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>

#define MyMalloc malloc
#define MyFree free
#define MyCalloc calloc
#define MyRealloc realloc

struct Matrix {
	float m00, m01, m02, m10, m11, m12;
};

struct Image {
	uint32_t width, height;
	uint32_t *bits;
};

void BlendPixel(uint32_t *destinationPixel, uint32_t modified);

#define STB_DS_IMPLEMENTATION
#include "stblib/stb_ds.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stblib/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stblib/stb_image_write.h"
#include "software_rasterizer.cpp"
#include "textlib/text.h"

struct INIState {
	char *buffer, *section, *key, *value;
	size_t bytes, sectionBytes, keyBytes, valueBytes;
};

struct Property {
	char *key;
	char *value;
};

struct Node {
	char *key;
	Property *value;
};

Node *nodes;

char emptyString;

bool INIParse(INIState *s) {
#define INI_READ(destination, counter, c1, c2) \
	s->destination = s->buffer, s->counter = 0; \
	while (s->bytes && *s->buffer != c1 && *s->buffer != c2) s->counter++, s->buffer++, s->bytes--; \
	if (s->bytes && *s->buffer == c1) s->buffer++, s->bytes--;

	while (s->bytes) {
		char c = *s->buffer;

		if (c == ' ' || c == '\n' || c == '\r') { 
			s->buffer++, s->bytes--; 
			continue;
		} else if (c == ';') {
			s->valueBytes = 0;
			INI_READ(key, keyBytes, '\n', 0);
		} else if (c == '[') {
			s->keyBytes = s->valueBytes = 0;
			s->buffer++, s->bytes--;
			INI_READ(section, sectionBytes, ']', 0);
		} else {
			INI_READ(key, keyBytes, '=', '\n');
			INI_READ(value, valueBytes, '\n', 0);
		}

		if (s->sectionBytes) s->section[s->sectionBytes] = 0; else s->section = &emptyString;
		if (s->keyBytes) s->key[s->keyBytes] = 0; else s->key = &emptyString;
		if (s->valueBytes) s->value[s->valueBytes] = 0; else s->value = &emptyString;

		return true;
	}

	return false;
}

void *LoadFile(const char *path, size_t *length) {
	FILE *file = fopen(path, "rb");
	if (!file) return NULL;
	fseek(file, 0, SEEK_END);
	size_t fileSize = ftell(file);
	fseek(file, 0, SEEK_SET);
	char *buffer = (char *) malloc(fileSize + 1);
	buffer[fileSize] = 0;
	fread(buffer, 1, fileSize, file);
	fclose(file);
	if (length) *length = fileSize;
	return buffer;
}

void BlendPixel(uint32_t *destinationPixel, uint32_t modified) {
	// TODO Is Figma doing gamma-correct blending?
	// 	rast-full.cpp and textlib will also need to be updated.

	if ((modified & 0xFF000000) == 0xFF000000) {
		*destinationPixel = modified;
		return;
	} else if ((modified & 0xFF000000) == 0x00000000) {
		return;
	}

	uint32_t m1, m2, a;
	uint32_t original = *destinationPixel;

	if ((*destinationPixel & 0xFF000000) != 0xFF000000) {
		uint32_t alpha1 = (modified & 0xFF000000) >> 24;
		uint32_t alpha2 = 255 - alpha1;
		uint32_t alphaD = (original & 0xFF000000) >> 24;
		uint32_t alphaD2 = alphaD * alpha2;
		uint32_t alphaOut = alpha1 + (alphaD2 >> 8);

		if (!alphaOut) {
			return;
		}

		m2 = alphaD2 / alphaOut;
		m1 = (alpha1 << 8) / alphaOut;
		if (m2 == 0x100) m2--;
		if (m1 == 0x100) m1--;
		a = alphaOut << 24;
	} else {
		m1 = (modified & 0xFF000000) >> 24;
		m2 = 255 - m1;
		a = 0xFF000000;
	}

	uint32_t r2 = m2 * (original & 0x00FF00FF);
	uint32_t g2 = m2 * (original & 0x0000FF00);
	uint32_t r1 = m1 * (modified & 0x00FF00FF);
	uint32_t g1 = m1 * (modified & 0x0000FF00);
	uint32_t result = a | (0x0000FF00 & ((g1 + g2) >> 8)) | (0x00FF00FF & ((r1 + r2) >> 8));
	*destinationPixel = result;
}

void Blit(Image *destination, Image *source, int32_t x0, int32_t y0, bool blend, float opacity) {
	for (uint32_t y = 0; y < source->height; y++) {
		for (uint32_t x = 0; x < source->width; x++) {
			int32_t x1 = x + x0, y1 = y + y0;

			if (x1 >= 0 && x1 < (int32_t) destination->width && y1 >= 0 && y1 < (int32_t) destination->height) {
				if (blend) {
					uint32_t sample = source->bits[x + y * source->width];
					sample = (sample & 0xFFFFFF) | ((uint32_t) ((float) (sample >> 24) * opacity) << 24);
					BlendPixel(&destination->bits[x1 + y1 * destination->width], sample);
				} else {
					destination->bits[x1 + y1 * destination->width] = source->bits[x + y * source->width];
				}
			}
		}
	}
}

RastVertex Apply(Matrix matrix, RastVertex vertex) {
	return { matrix.m00 * vertex.x + matrix.m01 * vertex.y + matrix.m02,
		matrix.m10 * vertex.x + matrix.m11 * vertex.y + matrix.m12 };
}

Matrix Multiply(Matrix a, Matrix b) {
	return (Matrix) {
		a.m00 * b.m00 + a.m01 * b.m10,
		a.m00 * b.m01 + a.m01 * b.m11,
		a.m00 * b.m02 + a.m01 * b.m12 + a.m02,
		a.m10 * b.m00 + a.m11 * b.m10,
		a.m10 * b.m01 + a.m11 * b.m11,
		a.m10 * b.m02 + a.m11 * b.m12 + a.m12,
	};
}

Matrix Inverse(Matrix a) {
	float dm1 = a.m00 * a.m11 - a.m01 * a.m10;

	return (Matrix) {
		a.m11 * dm1,
		a.m01 * dm1 * -1.0f,
		(a.m01 * a.m12 - a.m11 * a.m02) * dm1,
		a.m10 * dm1 * -1.0f,
		a.m00 * dm1,
		(a.m00 * a.m12 - a.m10 * a.m02) * dm1 * -1.0f,
	};
}

RastPaint LoadPaint(Property *cp, Matrix absoluteTransform, float width, float height, uint32_t **outImage) {
	RastPaint paint = {};
	const char *type = shget(cp, "type");

	if (0 == strcmp(type, "PAINT.SOLID")) {
		char *color = shget(cp, "color");
		float red = strtod(color, &color);
		color++;
		float green = strtod(color, &color);
		color++;
		float blue = strtod(color, &color);
		paint.type = RAST_PAINT_SOLID;
		paint.solid.color = ((uint32_t) (red * 255.0f) << 16) 
			| ((uint32_t) (green * 255.0f) << 8) 
			| ((uint32_t) (blue * 255.0f) << 0);
		paint.solid.alpha = atof(shget(cp, "opacity"));
	} else if (0 == strcmp(type, "PAINT.GRADIENT_LINEAR") || 0 == strcmp(type, "PAINT.GRADIENT_RADIAL") || 0 == strcmp(type, "PAINT.GRADIENT_ANGULAR")) {
		paint.type = 0 == strcmp(type, "PAINT.GRADIENT_RADIAL") ? RAST_PAINT_RADIAL_GRADIENT : 
			0 == strcmp(type, "PAINT.GRADIENT_ANGULAR") ? RAST_PAINT_ANGULAR_GRADIENT_2 : RAST_PAINT_LINEAR_GRADIENT;

		float opacity0 = atof(shget(cp, "opacity"));

		Matrix gradientTransform = {
			(float) atof(shget(cp, "gradientTransform00")) / width,
			(float) atof(shget(cp, "gradientTransform01")) / height,
			(float) atof(shget(cp, "gradientTransform02")),
			(float) atof(shget(cp, "gradientTransform10")) / width,
			(float) atof(shget(cp, "gradientTransform11")) / height,
			(float) atof(shget(cp, "gradientTransform12")),
		};

		Matrix result = Multiply(gradientTransform, Inverse(absoluteTransform));
		memcpy(&paint.gradient.transform[0], &result, sizeof(Matrix));

		if (paint.type != RAST_PAINT_LINEAR_GRADIENT) {
			paint.gradient.transform[0] *= 2.0f;
			paint.gradient.transform[1] *= 2.0f;
			paint.gradient.transform[2] *= 2.0f;
			paint.gradient.transform[3] *= 2.0f;
			paint.gradient.transform[4] *= 2.0f;
			paint.gradient.transform[5] *= 2.0f;
			paint.gradient.transform[2] -= 1.0f;
			paint.gradient.transform[5] -= 1.0f;
		}

		RastGradientStop *stops = {};
		char *string = shget(cp, "gradientStops");

		while (*string) {
			float red = strtod(string, &string);
			string++;
			float green = strtod(string, &string);
			string++;
			float blue = strtod(string, &string);
			string++;
			float opacity = opacity0 * strtod(string, &string);
			string++;
			float position = strtod(string, &string);
			string++;
			RastGradientStop stop = {};
			stop.color = ((uint32_t) (opacity * 255.0f) << 24) 
				| ((uint32_t) (red   * 255.0f) << 16) 
				| ((uint32_t) (green * 255.0f) <<  8) 
				| ((uint32_t) (blue  * 255.0f) <<  0);
			stop.position = position;
			arrput(stops, stop);
		}

		RastGradientInitialise(&paint, stops, arrlenu(stops), false);
		arrfree(stops);
	} else if (0 == strcmp(type, "PAINT.IMAGE")) {
		// TODO Scale modes: fill, fit, crop, tile.
		// TODO Rotation.
		// TODO Filters.

		paint.type = RAST_PAINT_IMAGE;

		char *data = shget(cp, "imageData");
		size_t dataBytes = strlen(data) / 2;
		uint8_t *buffer = (uint8_t *) malloc(dataBytes);

		for (uintptr_t i = 0; i < dataBytes; i++) {
			buffer[i] = (data[i * 2 + 0] >= 'A' ? data[i * 2 + 0] - 'A' + 10 : data[i * 2 + 0] - '0') 
				 | ((data[i * 2 + 1] >= 'A' ? data[i * 2 + 1] - 'A' + 10 : data[i * 2 + 1] - '0') << 4);
		}

		int imageWidth, imageHeight, channels;
		*outImage = (uint32_t *) stbi_load_from_memory(buffer, dataBytes, &imageWidth, &imageHeight, &channels, 4);

		for (intptr_t i = 0; i < imageWidth * imageHeight; i++) {
			(*outImage)[i] = ((*outImage)[i] & 0xFF00FF00) | (((*outImage)[i] & 0xFF0000) >> 16) | (((*outImage)[i] & 0xFF) << 16);
		}

		paint.image.bits = *outImage;
		paint.image.width = imageWidth;
		paint.image.height = imageHeight;
		free(buffer);
		assert(*outImage && imageWidth && imageHeight);

		bool mode = (float) imageWidth / (float) imageHeight < (float) width / (float) height;
		float  scaleX = mode == false ? (float) height / (float) width : 1.0f;
		float  scaleY = mode ==  true ? (float) width / (float) height : 1.0f;
		float offsetX = mode == false ? (((imageWidth  * height / imageHeight) * 0.5f - width  * 0.5f) * (imageWidth  / width))  : 0.0f;
		float offsetY = mode ==  true ? (((imageHeight * width  / imageWidth)  * 0.5f - height * 0.5f) * (imageHeight / height)) : 0.0f;

		// TODO "Nearest neighbour" offset.
		Matrix imageTransform = { scaleX / width, 0.0f, offsetX / width, 0.0f, scaleY / height, offsetY / height };
		Matrix result = Multiply(imageTransform, Inverse(absoluteTransform));
		memcpy(&paint.image.transform[0], &result, sizeof(Matrix));
	}

	return paint;
}

void PathAppendCorner(RastPath *path, RastVertex v1, RastVertex v2, RastVertex center) {
	// Because Figma's transforms have matrix columns of unit length,
	// the radius of the arc is constant, thus we may use RastPathAppendArc.
	_RastPathAddVertex(path, v1);
	float radius = sqrtf((v1.x - center.x) * (v1.x - center.x) + (v1.y - center.y) * (v1.y - center.y));
	float startAngle = atan2f(center.y - v1.y, center.x - v1.x) + M_PI;
	float endAngle = atan2f(center.y - v2.y, center.x - v2.x) + M_PI;
	if (endAngle > startAngle) endAngle -= M_PI * 2.0f;
	if (radius > 0.1f) RastPathAppendArc(path, center, radius, startAngle, endAngle);
	_RastPathAddVertex(path, v2);
}

void Blur(float *buffer, uint32_t width, uint32_t height, float radius, bool clamp /* if false, use void color */) {
	float *voidColor = &buffer[width * height * 4];
	float *buffer2 = (float *) malloc(width * height * sizeof(float) * 4);
	float standardDeviation = radius * 0.5f; // Assuming the radius property is interpreted the same as in the CSS specification.
						 
	// Calculate the kernel.

	int32_t kernelSize = (int32_t) ceilf(standardDeviation * 3.0f);
	float *kernel = (float *) malloc((kernelSize * 2 + 1) * sizeof(float));
	float kernelSum = 0.0f;

	for (int32_t x = -kernelSize; x <= kernelSize; x++) {
		kernel[x + kernelSize] = expf(x * x * -0.5f / (standardDeviation * standardDeviation)) / (2.506628275f * standardDeviation);
		kernelSum += kernel[x + kernelSize];
	}

	float kernelSumReciprocal = 1.0f / kernelSum;

	// Blur vertically.

	for (uint32_t column = 0; column < width; column++) {
		for (uint32_t row = 0; row < height; row++) {
			float sums[4] = {};

			for (int32_t x0 = -kernelSize; x0 <= kernelSize; x0++) {
				int32_t x = row + x0;

				if (clamp) {
					if (x < 0) x = 0;
					if (x >= (int32_t) height) x = height - 1;
				}
				
				if (x >= 0 && x < (int32_t) height) {
					sums[0] += buffer[(column + x * width) * 4 + 0] * kernel[x0 + kernelSize];
					sums[1] += buffer[(column + x * width) * 4 + 1] * kernel[x0 + kernelSize];
					sums[2] += buffer[(column + x * width) * 4 + 2] * kernel[x0 + kernelSize];
					sums[3] += buffer[(column + x * width) * 4 + 3] * kernel[x0 + kernelSize];
				} else {
					sums[0] += voidColor[0] * kernel[x0 + kernelSize];
					sums[1] += voidColor[1] * kernel[x0 + kernelSize];
					sums[2] += voidColor[2] * kernel[x0 + kernelSize];
					sums[3] += voidColor[3] * kernel[x0 + kernelSize];
				}
			}

			buffer2[(column + row * width) * 4 + 0] = sums[0] * kernelSumReciprocal;
			buffer2[(column + row * width) * 4 + 1] = sums[1] * kernelSumReciprocal;
			buffer2[(column + row * width) * 4 + 2] = sums[2] * kernelSumReciprocal;
			buffer2[(column + row * width) * 4 + 3] = sums[3] * kernelSumReciprocal;
		}
	}

	// Blur horizontally.

	for (uint32_t row = 0; row < height; row++) {
		for (uint32_t column = 0; column < width; column++) {
			float sums[4] = {};

			for (int32_t x0 = -kernelSize; x0 <= kernelSize; x0++) {
				int32_t x = column + x0;

				if (clamp) {
					if (x < 0) x = 0;
					if (x >= (int32_t) width) x = width - 1;
				}
				
				if (x >= 0 && x < (int32_t) width) {
					sums[0] += buffer2[(x + row * width) * 4 + 0] * kernel[x0 + kernelSize];
					sums[1] += buffer2[(x + row * width) * 4 + 1] * kernel[x0 + kernelSize];
					sums[2] += buffer2[(x + row * width) * 4 + 2] * kernel[x0 + kernelSize];
					sums[3] += buffer2[(x + row * width) * 4 + 3] * kernel[x0 + kernelSize];
				} else {
					sums[0] += voidColor[0] * kernel[x0 + kernelSize];
					sums[1] += voidColor[1] * kernel[x0 + kernelSize];
					sums[2] += voidColor[2] * kernel[x0 + kernelSize];
					sums[3] += voidColor[3] * kernel[x0 + kernelSize];
				}
			}

			buffer[(column + row * width) * 4 + 0] = sums[0] * kernelSumReciprocal;
			buffer[(column + row * width) * 4 + 1] = sums[1] * kernelSumReciprocal;
			buffer[(column + row * width) * 4 + 2] = sums[2] * kernelSumReciprocal;
			buffer[(column + row * width) * 4 + 3] = sums[3] * kernelSumReciprocal;
		}
	}

	// Free intermediate buffers.

	free(buffer2);
	free(kernel);
}

void RenderNode(const char *id, float drawXOffset, float drawYOffset, Image *destinationBitmap) {
	// TODO Clipping (including frames).
	// TODO Gradient dithering.

	Property *p = shget(nodes, id);

	// Check the node is visible.

	const char *visible = shget(p, "visible");
	if (visible && 0 == strcmp(visible, "false")) return;

	// Create the bitmap.

	// TODO The bitmap may need to be larger than this, otherwise offset drop shadows will not be rendered correctly.
	// 	In this case, the entire node needs to be visible, not just the portion that will be blitted onto the parent's buffer.
	Image bitmap = {};
	bitmap.width = destinationBitmap->width;
	bitmap.height = destinationBitmap->height, 
	bitmap.bits = (uint32_t *) calloc(4, destinationBitmap->width * destinationBitmap->height);
	RastSurface surface = {};
	surface.buffer = (uint32_t *) bitmap.bits;
	surface.stride = bitmap.width * 4;
	RastSurfaceInitialise(&surface, bitmap.width, bitmap.height, true);
	RastSurface mask = {};
	RastSurfaceInitialise(&mask, bitmap.width, bitmap.height, false);

	// Get the transform.

	Matrix matrix = {};
	float width = 0.0f, height = 0.0f;

	if (shget(p, "absoluteTransform00")) {
		matrix.m00 = atof(shget(p, "absoluteTransform00"));
		matrix.m01 = atof(shget(p, "absoluteTransform01"));
		matrix.m02 = atof(shget(p, "absoluteTransform02")) + drawXOffset;
		matrix.m10 = atof(shget(p, "absoluteTransform10"));
		matrix.m11 = atof(shget(p, "absoluteTransform11"));
		matrix.m12 = atof(shget(p, "absoluteTransform12")) + drawYOffset;
		width = atof(shget(p, "width"));
		height = atof(shget(p, "height"));
	}

	// Get the basic properties and arrays.

	const char *type = shget(p, "type");
	const char *fills = shget(p, "fills");
	const char *children = shget(p, "children");
	const char *effects = shget(p, "effects");
	const char *strokes = shget(p, "strokes");

	RastPath path = {};
	RastContourStyle contourStyle = {};
	bool evenOdd = false;

	// Set the path.

	// TODO Corner smoothing, see https://www.figma.com/blog/desperately-seeking-squircles/.

	if (0 == strcmp(type, "RECTANGLE")) {
		// TODO Different contour sizes on each edge.

		float maximumCornerRadius = width * 0.5f > height * 0.5f ? height * 0.5f : width * 0.5f;
		float topLeftRadius = atof(shget(p, "topLeftRadius"));
		float topRightRadius = atof(shget(p, "topRightRadius"));
		float bottomLeftRadius = atof(shget(p, "bottomLeftRadius"));
		float bottomRightRadius = atof(shget(p, "bottomRightRadius"));
		if (topLeftRadius > maximumCornerRadius) topLeftRadius = maximumCornerRadius;
		if (topRightRadius > maximumCornerRadius) topRightRadius = maximumCornerRadius;
		if (bottomLeftRadius > maximumCornerRadius) bottomLeftRadius = maximumCornerRadius;
		if (bottomRightRadius > maximumCornerRadius) bottomRightRadius = maximumCornerRadius;

		PathAppendCorner(&path, Apply(matrix, { 0, height - bottomLeftRadius }), Apply(matrix, { bottomLeftRadius, height }), 
				Apply(matrix, { bottomLeftRadius, height - bottomRightRadius }));
		PathAppendCorner(&path, Apply(matrix, { width - bottomRightRadius, height }), Apply(matrix, { width, height - bottomRightRadius }), 
				Apply(matrix, { width - bottomRightRadius, height - bottomRightRadius }));
		PathAppendCorner(&path, Apply(matrix, { width, topRightRadius }), Apply(matrix, { width - topRightRadius, 0 }), 
				Apply(matrix, { width - topRightRadius, topRightRadius }));
		PathAppendCorner(&path, Apply(matrix, { topLeftRadius, 0 }), Apply(matrix, { 0, topLeftRadius }), 
				Apply(matrix, { topLeftRadius, topLeftRadius }));
	} else if (0 == strcmp(type, "POLYGON")) {
		size_t pointCount = atoi(shget(p, "pointCount"));
		float cornerRadius = atof(shget(p, "cornerRadius"));

		if (cornerRadius >= 1.0f) {
			for (uintptr_t i = 0; i < pointCount; i++) {
				RastVertex vertex0 = Apply(matrix, { (sinf(M_PI * 2.0f * (i - 1.0f) / pointCount) * -0.5f + 0.5f) * width, 
					(cosf(M_PI * 2.0f * (i - 1.0f) / pointCount) * -0.5f + 0.5f) * height });
				RastVertex vertex2 = Apply(matrix, { (sinf(M_PI * 2.0f * (i + 1.0f) / pointCount) * -0.5f + 0.5f) * width, 
					(cosf(M_PI * 2.0f * (i + 1.0f) / pointCount) * -0.5f + 0.5f) * height });
				RastVertex vertex1 = Apply(matrix, { (sinf(M_PI * 2.0f * (i + 0.0f) / pointCount) * -0.5f + 0.5f) * width, 
					(cosf(M_PI * 2.0f * (i + 0.0f) / pointCount) * -0.5f + 0.5f) * height });

				RastVertex direction0 = { vertex0.x - vertex1.x, vertex0.y - vertex1.y };
				float length0 = sqrtf(direction0.x * direction0.x + direction0.y * direction0.y);
				direction0.x /= length0; direction0.y /= length0;

				RastVertex direction2 = { vertex2.x - vertex1.x, vertex2.y - vertex1.y };
				float length2 = sqrtf(direction2.x * direction2.x + direction2.y * direction2.y);
				direction2.x /= length2; direction2.y /= length2;

				RastVertex direction1 = { direction0.x + direction2.x, direction0.y + direction2.y };
				float length1 = sqrtf(direction1.x * direction1.x + direction1.y * direction1.y);
				direction1.x /= length1; direction1.y /= length1;

				float theta = atan2f(direction2.y, direction2.x) - atan2f(direction1.y, direction1.x);
				float tanTheta = tanf(theta);
				float dt = cornerRadius / tanTheta;
				// TODO Clamp the cornerRadius (see Group 6).
				float dc = dt * tanTheta / sinf(theta);

				PathAppendCorner(&path, { vertex1.x + direction0.x * dt, vertex1.y + direction0.y * dt },
						{ vertex1.x + direction2.x * dt, vertex1.y + direction2.y * dt },
						{ vertex1.x + direction1.x * dc, vertex1.y + direction1.y * dc });
			}
		} else {
			for (uintptr_t i = 0; i < pointCount; i++) {
				RastVertex vertex1 = Apply(matrix, { (sinf(M_PI * 2.0f * (i + 0) / pointCount) * -0.5f + 0.5f) * width, 
					(cosf(M_PI * 2.0f * (i + 0) / pointCount) * -0.5f + 0.5f) * height });
				_RastPathAddVertex(&path, vertex1);
			}
		}
	} else if (0 == strcmp(type, "ELLIPSE")) {
		// TODO Arc starting angle, ending angle and inner radius.
		// TODO Corner radii.

		float radius0 = width * 0.5f;
		float radius1 = height * 0.5f;
		float radius0P = powf(radius0, -2.0f);
		float radius1P = powf(radius1, -2.0f);
		float radiusMaximum = radius0 > radius1 ? radius0 : radius1;
		float deltaAngle = acosf(1 - 0.5f * RAST_FLATTEN_TOLERANCE * RAST_FLATTEN_TOLERANCE / radiusMaximum / radiusMaximum); // From cosine rule.
		size_t steps = 2.0f * M_PI / deltaAngle;

		for (uintptr_t i = 0; i <= steps; i++) {
			float angle = 2.0f * M_PI * ((float) i / steps);
			float angleC = cosf(angle);
			float angleS = sinf(angle);
			float radius = powf(radius0P * angleC * angleC + radius1P * angleS * angleS, -0.5f);
			_RastPathAddVertex(&path, Apply(matrix, { radius0 + radius * angleC, radius1 + radius * angleS }));
		}
	} else if (0 == strcmp(type, "BOOLEAN_OPERATION") || 0 == strcmp(type, "VECTOR")) {
		children = nullptr;
		char *data = shget(p, "vectorPath");
		float px = 0.0f, py = 0.0f;

		while (*data) {
			char c = *data++;

			if (c == 'M' || c == 'L') {
				px = strtof(data, &data); data++;
				py = strtof(data, &data);
				_RastPathAddVertex(&path, Apply(matrix, { px, py }));
			} else if (c == 'C') {
				float x0 = strtof(data, &data); data++;
				float y0 = strtof(data, &data); data++;
				float x1 = strtof(data, &data); data++;
				float y1 = strtof(data, &data); data++;
				float x  = strtof(data, &data); data++;
				float y  = strtof(data, &data);
				_RastFlattenBezierRecursive(&path, Apply(matrix, { px, py }), Apply(matrix, { x0, y0 }), 
						Apply(matrix, { x1, y1 }), Apply(matrix, { x, y }), 0);
				px = x, py = y;
			} else if (c == 'Z') {
				RastPathCloseSegment(&path);
			}
		}

		evenOdd = true;
	} else if (0 == strcmp(type, "TEXT")) {
		// TODO This needs to be converted to a vector path and rendered as other nodes!
		// 	At the moment, things like gradients, strokes, rotations, etc. don't work.
		// TODO Something's off with the kerning/hinting...
		// TODO Justified text, bidi, etc.

		const char *alignH = shget(p, "textAlignHorizontal");
		const char *alignV = shget(p, "textAlignVertical");

		EsTextPlanProperties properties = {};
		properties.flags = ES_TEXT_PLAN_SINGLE_USE | ES_TEXT_PLAN_NO_FONT_SUBSTITUTION | ES_TEXT_WRAP
			| (0 == strcmp(alignH, "LEFT") ? ES_TEXT_H_LEFT : 0 == strcmp(alignH, "RIGHT") ? ES_TEXT_H_RIGHT : ES_TEXT_H_CENTER)
			| (0 == strcmp(alignV, "TOP") ? ES_TEXT_V_TOP : 0 == strcmp(alignV, "BOTTOM") ? ES_TEXT_V_BOTTOM : ES_TEXT_V_CENTER);

		EsTextRun *textRuns = nullptr;
		char *string = nullptr;

		char *copyBase = strdup(shget(p, "textRuns"));
		char *copy = copyBase;

		while (*copy) {
			const char *n = copy;
			strtol(copy, &copy, 10);
			if (*copy == ',') { *copy = 0; copy++; }
			Property *cp = shget(nodes, n);

			char *fillsCopy = strdup(shget(cp, "fills"));
			char *fillsCopyEnd;
			strtol(fillsCopy, &fillsCopyEnd, 10);
			*fillsCopyEnd = 0;
			Property *fp = shget(nodes, fillsCopy);
			free(fillsCopy);

			const char *fontStyle = shget(cp, "fontStyle");

			char *color = shget(fp, "color");
			float red = strtod(color, &color);
			color++;
			float green = strtod(color, &color);
			color++;
			float blue = strtod(color, &color);
			float opacity = atof(shget(fp, "opacity"));

			EsFontInformation fontInformation;
			bool foundFont = EsFontDatabaseLookupByName(shget(cp, "fontFamily"), -1, &fontInformation);

			EsTextRun run = {};
			run.offset = arrlenu(string);
			run.style.font.family = foundFont ? fontInformation.id : 0;
			run.style.font.weight = strstr(fontStyle, "Thin") ? 1 : strstr(fontStyle, "Extra Light") ? 2 : strstr(fontStyle, "Light") ? 3
				: strstr(fontStyle, "Black") ? 9 : strstr(fontStyle, "Extra Bold") ? 8 : strstr(fontStyle, "Semi Bold") ? 6
				: strstr(fontStyle, "SemiBold") ? 6 : strstr(fontStyle, "Bold") ? 7 : 4;
			run.style.font.flags = strstr(fontStyle, "Italic") ? ES_FONT_ITALIC : 0;
			run.style.size = atoi(shget(cp, "fontSize")) - 5.0f / 64.0f /* Figma seems to favour slightly smaller font sizes */;
			run.style.decorations = 0 == strcmp(shget(cp, "textDecoration"), "UNDERLINE") ? ES_TEXT_DECORATION_UNDERLINE
				: 0 == strcmp(shget(cp, "textDecoration"), "STRIKETHROUGH") ? ES_TEXT_DECORATION_STRIKE_THROUGH : 0;
			run.style.color = run.style.decorationsColor = ((uint32_t) (opacity * 255.0f) << 24) | ((uint32_t) (red * 255.0f) << 16) 
				| ((uint32_t) (green * 255.0f) << 8) | ((uint32_t) (blue * 255.0f) << 0);
			// TODO Increase tracking slightly to match Figma's rendering better?

			const char *characters = shget(cp, "characters");
			assert(*characters == '"');
			characters++;

			while (*characters != '"') {
				if (*characters == '\\') {
					// TODO Unicode.
					arrput(string, characters[1] == 'n' ? '\n' : characters[1]);
					characters += 2;
				} else {
					arrput(string, *characters++);
				}
			}

			arrput(textRuns, run);
		}

		free(copyBase);

		EsTextRun terminatingRun = {};
		terminatingRun.offset = arrlenu(string);
		arrput(textRuns, terminatingRun);

		EsRectangle bounds = ES_RECT_4(Apply(matrix, { 0, 0 }).x, Apply(matrix, { width, height }).x, 
				Apply(matrix, { 0, 0 }).y, Apply(matrix, { width, height }).y);
		EsTextPlan *plan = EsTextPlanCreate(nullptr, &properties, bounds, string, textRuns, arrlenu(textRuns) - 1);
		EsPaintTarget target = { 0 };
		target.bits = bitmap.bits;
		target.width = bitmap.width;
		target.height = bitmap.height;
		target.stride = bitmap.width * 4;
		target.fullAlpha = true;
		EsPainter painter = { 0 };
		painter.clip = ES_RECT_4(0, bitmap.width, 0, bitmap.height);
		painter.target = &target;
		EsDrawText(&painter, plan, bounds);

		// TODO This doesn't work for translucent text.
		for (uintptr_t i = 0; i < bitmap.height; i++) {
			for (uintptr_t j = 0; j < bitmap.width; j++) {
				((float *) mask.buffer)[i * bitmap.width + j] = (float) (bitmap.bits[i * bitmap.width + j] >> 24) / 255.0f;
			}
		}

		arrfree(textRuns);
		arrfree(string);
	}

	// Render fills.

	if (fills) {
		char *copy = strdup(fills);
		char *copyBase = copy;

		while (*copy) {
			const char *n = copy;
			strtol(copy, &copy, 10);
			if (*copy == ',') { *copy = 0; copy++; }

			uint32_t *image = nullptr;
			RastPaint paint = LoadPaint(shget(nodes, n), matrix, width, height, &image);
			RastSurfaceFill(surface, RastShapeCreateSolid(&path), paint, evenOdd);
			RastSurfaceFill(mask, RastShapeCreateSolid(&path), { .type = RAST_PAINT_MAXIMUM_A }, evenOdd);
			RastGradientDestroy(&paint);
			free(image);
		}

		free(copyBase);
	}

	// Render strokes.

	if (strokes) {
		char *copy = strdup(strokes);
		char *copyBase = copy;

		float strokeWeight = atof(shget(p, "strokeWeight"));
		const char *strokeAlign = shget(p, "strokeAlign");
		const char *strokeJoin = shget(p, "strokeJoin");
		const char *strokeCap = shget(p, "strokeCap");

		contourStyle.internalWidth = 0 == strcmp(strokeAlign, "INSIDE") ? strokeWeight 
			: 0 == strcmp(strokeAlign, "CENTER") ? strokeWeight * 0.5f : 0.0f;
		contourStyle.externalWidth = 0 == strcmp(strokeAlign, "OUTSIDE") ? strokeWeight 
			: 0 == strcmp(strokeAlign, "CENTER") ? strokeWeight * 0.5f : 0.0f;
		contourStyle.miterLimit = 0 == strcmp(strokeJoin, "BEVEL") ? 0.0f : atof(shget(p, "strokeMiterLimit"));
		contourStyle.joinMode = 0 == strcmp(strokeJoin, "ROUND") ? RAST_LINE_JOIN_ROUND : RAST_LINE_JOIN_MITER;
		contourStyle.capMode = 0 == strcmp(strokeCap, "ROUND") ? RAST_LINE_CAP_ROUND 
			: 0 == strcmp(strokeCap, "SQUARE") ? RAST_LINE_CAP_SQUARE : RAST_LINE_CAP_FLAT;
		// TODO Dashed; arrows.
		// TODO For internal opaque strokes, prevent the fill's colour leaking outside into the antialiasing of the stroke.

		while (*copy) {
			const char *n = copy;
			strtol(copy, &copy, 10);
			if (*copy == ',') { *copy = 0; copy++; }

			uint32_t *image = nullptr;
			RastPaint paint = LoadPaint(shget(nodes, n), matrix, width, height, &image);
			RastSurfaceFill(surface, RastShapeCreateContour(&path, contourStyle, false), paint, false);
			RastSurfaceFill(mask, RastShapeCreateContour(&path, contourStyle, false), { .type = RAST_PAINT_MAXIMUM_A }, false);
			RastGradientDestroy(&paint);
			free(image);
		}

		free(copyBase);
	}

	// Render the children.

	if (children) {
		char *copy = strdup(children);
		char *copyBase = copy;

		while (*copy) {
			const char *n = copy;
			strtol(copy, &copy, 10);
			if (*copy == ',') { *copy = 0; copy++; }
			RenderNode(n, drawXOffset, drawYOffset, &bitmap);
		}

		free(copyBase);
	}

	// Apply shadow effects.

	if (effects) {
		char *copy = strdup(effects);
		char *copyBase = copy;

		while (*copy) {
			const char *n = copy;
			strtol(copy, &copy, 10);
			if (*copy == ',') { *copy = 0; copy++; }

			Property *cp = shget(nodes, n);
			const char *type = shget(cp, "type");

			if (0 == strcmp(shget(cp, "visible"), "false")) {
				// The effect is not visible.
			} else if (0 == strcmp(type, "EFFECT.DROP_SHADOW") || 0 == strcmp(type, "EFFECT.INNER_SHADOW")) {
				bool inner = 0 == strcmp(type, "EFFECT.INNER_SHADOW");
				bool showShadowBehindNode = !inner && 0 == strcmp(shget(cp, "showShadowBehindNode"), "true");
				char *color = shget(cp, "color");
				float red = strtod(color, &color); color++;
				float green = strtod(color, &color); color++;
				float blue = strtod(color, &color); color++;
				float opacity = strtod(color, &color); color++;
				int offsetX = atoi(shget(cp, "offsetX")); // TODO This could be fractional.
				int offsetY = atoi(shget(cp, "offsetY")); // TODO This could be fractional.
				float radius = atof(shget(cp, "radius"));
				float spread = atof(shget(cp, "spread")); // TODO. (Rectangular shapes only.)
				(void) spread;

				int left = offsetX > 0 ? 0 : offsetX;
				int right = offsetX + bitmap.width > bitmap.width ? offsetX + bitmap.width : bitmap.width;
				int top = offsetY > 0 ? 0 : offsetY;
				int bottom = offsetY + bitmap.height > bitmap.height ? offsetY + bitmap.height : bitmap.height;
									  
				uint32_t bufferWidth = right - left;
				uint32_t bufferHeight = bottom - top;
				float *buffer = (float *) malloc((bufferWidth * bufferHeight + 1) * sizeof(float) * 4);

				for (uint32_t y = 0; y < bufferHeight; y++) {
					for (uint32_t x = 0; x < bufferWidth; x++) {
						float sample = 0.0f;
						int x0 = x - offsetX + left;
						int y0 = y - offsetY + top;

						if (x0 >= 0 && x0 < (int32_t) bitmap.width && y0 >= 0 && y0 < (int32_t) bitmap.height) {
							sample = ((float *) mask.buffer)[x0 + y0 * bitmap.width];
						}

						buffer[(x + y * bufferWidth) * 4 + 0] = blue;
						buffer[(x + y * bufferWidth) * 4 + 1] = green;
						buffer[(x + y * bufferWidth) * 4 + 2] = red;
						buffer[(x + y * bufferWidth) * 4 + 3] = opacity * (inner ? 1.0f - sample : sample);
					}
				}

				// Void color.
				buffer[bufferWidth * bufferHeight * 4 + 0] = blue;
				buffer[bufferWidth * bufferHeight * 4 + 1] = green;
				buffer[bufferWidth * bufferHeight * 4 + 2] = red;
				buffer[bufferWidth * bufferHeight * 4 + 3] = opacity * (inner ? 1.0f : 0.0f);

				if (radius) {
					Blur(buffer, bufferWidth, bufferHeight, radius, false /* use the void color */);
				}

				if (!showShadowBehindNode && !inner) {
					for (uint32_t y = 0; y < bitmap.height; y++) {
						for (uint32_t x = 0; x < bitmap.width; x++) {
							float sample = 1.0f - ((float *) mask.buffer)[y * bitmap.width + x];
							buffer[((y - top) * bufferWidth + (x - left)) * 4 + 3] *= sample;
						}
					}
				}

				if (inner) {
					for (uint32_t y = 0; y < bitmap.height; y++) {
						for (uint32_t x = 0; x < bitmap.width; x++) {
							float sample = ((float *) mask.buffer)[y * bitmap.width + x];
							buffer[((y - top) * bufferWidth + (x - left)) * 4 + 3] *= sample;
						}
					}
				}

				Image asImage = {};
				asImage.width = bufferWidth;
				asImage.height = bufferHeight;
				asImage.bits = (uint32_t *) malloc(bufferWidth * bufferHeight * sizeof(uint32_t));

				for (uintptr_t i = 0; i < bufferWidth * bufferHeight; i++) {
					asImage.bits[i] = ((uint32_t) (buffer[i * 4 + 0] * 255.0f) << 0)
						| ((uint32_t) (buffer[i * 4 + 1] * 255.0f) << 8)
						| ((uint32_t) (buffer[i * 4 + 2] * 255.0f) << 16)
						| ((uint32_t) (buffer[i * 4 + 3] * 255.0f) << 24);
				}

				// TODO Blend modes.

				if (inner) {
					Blit(&bitmap, &asImage, left, top, true, 1.0f);
				} else {
					Blit(&asImage, &bitmap, -left, -top, true, 1.0f);
					Blit(&bitmap, &asImage, left, top, false, 1.0f);
				}

				free(asImage.bits);
				free(buffer);
			}
		}

		free(copyBase);
	}

	// Apply the layer blur effect (if it exists).

	if (effects) {
		char *copy = strdup(effects);
		char *copyBase = copy;

		while (*copy) {
			const char *n = copy;
			strtol(copy, &copy, 10);
			if (*copy == ',') { *copy = 0; copy++; }

			Property *cp = shget(nodes, n);
			const char *type = shget(cp, "type");

			if (0 == strcmp(shget(cp, "visible"), "false")) {
				// The effect is not visible.
			} else if (0 == strcmp(type, "EFFECT.LAYER_BLUR")) {
				float radius = atof(shget(cp, "radius"));
				float *buffer = (float *) malloc((bitmap.width * bitmap.height + 1) * sizeof(float) * 4);

				for (uint32_t y = 0; y < bitmap.height; y++) {
					for (uint32_t x = 0; x < bitmap.width; x++) {
						buffer[(x + y * bitmap.width) * 4 + 0] = (float) (0xFF & (bitmap.bits[x + y * bitmap.width] >>  0)) / 255.0f;
						buffer[(x + y * bitmap.width) * 4 + 1] = (float) (0xFF & (bitmap.bits[x + y * bitmap.width] >>  8)) / 255.0f;
						buffer[(x + y * bitmap.width) * 4 + 2] = (float) (0xFF & (bitmap.bits[x + y * bitmap.width] >> 16)) / 255.0f;
						buffer[(x + y * bitmap.width) * 4 + 3] = (float) (0xFF & (bitmap.bits[x + y * bitmap.width] >> 24)) / 255.0f;
					}
				}

				if (radius) {
					Blur(buffer, bitmap.width, bitmap.height, radius, true /* clamp to edges */);
				}

				for (uintptr_t i = 0; i < bitmap.width * bitmap.height; i++) {
					bitmap.bits[i] = ((uint32_t) (buffer[i * 4 + 0] * 255.0f) << 0)
						| ((uint32_t) (buffer[i * 4 + 1] * 255.0f) << 8)
						| ((uint32_t) (buffer[i * 4 + 2] * 255.0f) << 16)
						| ((uint32_t) (buffer[i * 4 + 3] * 255.0f) << 24);
				}

				free(buffer);
			}
		}

		free(copyBase);
	}

	// TODO Apply the background blur effect (if it exists).

	// Blit on the destination image.

	// TODO Blend modes.
	float opacity = atof(shget(p, "opacity"));
	Blit(destinationBitmap, &bitmap, 0, 0, true, opacity);

	// Free temporary data.

	RastPathDestroy(&path);
	free(bitmap.bits);
	RastSurfaceDestroy(&surface);
	RastSurfaceDestroy(&mask);
}

int main(int argc, char **argv) {
	if (argc != 4) {
		fprintf(stderr, "Usage: %s <input file> <node name> <output file>\n", argv[0]);
		return 1;
	}

	size_t inputLength;
	void *input = LoadFile(argv[1], &inputLength);

	if (!input) {
		fprintf(stderr, "Error: Couldn't read input file.\n");
		return 1;
	}

	INIState s = { .buffer = (char *) input, .bytes = inputLength };
	bool keepGoing = true;
	INIParse(&s);

	while (keepGoing) {
		Property *properties = nullptr;
		const char *section = s.section;

		while (true) {
			keepGoing = INIParse(&s);
			if (!keepGoing || !s.keyBytes) break;
			shput(properties, s.key, s.value);
		}

		shput(nodes, section, properties);
	}

	{
		const char *fontFiles[] = { 
			"res/Fonts/Barlow-Thin.ttf", "res/Fonts/Barlow-ThinItalic.ttf",
			"res/Fonts/Barlow-ExtraLight.ttf", "res/Fonts/Barlow-ExtraLightItalic.ttf",
			"res/Fonts/Barlow-Light.ttf", "res/Fonts/Barlow-LightItalic.ttf",
			"res/Fonts/Barlow-Regular.ttf", "res/Fonts/Barlow-Italic.ttf",
			"res/Fonts/Barlow-Medium.ttf", "res/Fonts/Barlow-MediumItalic.ttf",
			"res/Fonts/Barlow-SemiBold.ttf", "res/Fonts/Barlow-SemiBoldItalic.ttf",
			"res/Fonts/Barlow-Bold.ttf", "res/Fonts/Barlow-BoldItalic.ttf",
			"res/Fonts/Barlow-ExtraBold.ttf", "res/Fonts/Barlow-ExtraBoldItalic.ttf",
			"res/Fonts/Barlow-Black.ttf", "res/Fonts/Barlow-BlackItalic.ttf",
		};

		const char *fontTypes[] = { 
			"1", "1i",
			"2", "2i",
			"3", "3i",
			"4", "4i",
			"5", "5i",
			"6", "6i",
			"7", "7i",
			"8", "8i",
			"9", "9i",
		};

		TextAddFont("Barlow", "Sans", "Latn,Grek,Cyrl", fontFiles, fontTypes, sizeof(fontFiles) / sizeof(fontFiles[0]));
	}

	{
		const char *fontFiles[] = { 
			"res/Fonts/Inter Thin.otf", "res/Fonts/Inter Thin Italic.otf",
			"res/Fonts/Inter Extra Light.otf", "res/Fonts/Inter Extra Light Italic.otf",
			"res/Fonts/Inter Light.otf", "res/Fonts/Inter Light Italic.otf",
			"res/Fonts/Inter Regular.otf", "res/Fonts/Inter Italic.otf",
			"res/Fonts/Inter Medium.otf", "res/Fonts/Inter Medium Italic.otf",
			"res/Fonts/Inter Semi Bold.otf", "res/Fonts/Inter Semi Bold Italic.otf",
			"res/Fonts/Inter Bold.otf", "res/Fonts/Inter Bold Italic.otf",
			"res/Fonts/Inter Extra Bold.otf", "res/Fonts/Inter Extra Bold Italic.otf",
			"res/Fonts/Inter Black.otf", "res/Fonts/Inter Black Italic.otf",
		};

		const char *fontTypes[] = { 
			"1", "1i",
			"2", "2i",
			"3", "3i",
			"4", "4i",
			"5", "5i",
			"6", "6i",
			"7", "7i",
			"8", "8i",
			"9", "9i",
		};

		TextAddFont("Inter", "Sans", "Latn,Grek,Cyrl", fontFiles, fontTypes, sizeof(fontFiles) / sizeof(fontFiles[0]));
	}

	for (uintptr_t i = 0; i < shlenu(nodes); i++) {
		const char *name = shget(nodes[i].value, "name");

		if (name && 0 == strcmp(argv[2], name)) {
			Image bitmap = {};
			bitmap.width = atoi(shget(nodes[i].value, "width"));
			bitmap.height = atoi(shget(nodes[i].value, "height"));
			bitmap.bits = (uint32_t *) calloc(sizeof(uint32_t), bitmap.width * bitmap.height);
			RenderNode(nodes[i].key, -atof(shget(nodes[i].value, "absoluteTransform02")), 
					-atof(shget(nodes[i].value, "absoluteTransform12")), &bitmap);

			for (uint32_t i = 0; i < bitmap.width * bitmap.height; i++) {
				bitmap.bits[i] = (bitmap.bits[i] & 0xFF00FF00) | ((bitmap.bits[i] & 0xFF) << 16) | ((bitmap.bits[i] >> 16) & 0xFF);
			}

			stbi_write_png(argv[3], bitmap.width, bitmap.height, 4, bitmap.bits, bitmap.width * 4);
			free(bitmap.bits);
			break;
		}
	}

	for (uintptr_t i = 0; i < shlenu(nodes); i++) {
		shfree(nodes[i].value);
	}

	shfree(nodes);
	free(input);

	return 0;
}
