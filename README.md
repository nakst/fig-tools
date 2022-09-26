# fig-tools

Please note that this is very work in progress! If you have time to spare, please help out contributing features and resolving issues. Many are marked with `TODO` in the code.

## Exporting a document

1. Open your document in Figma, and open your web browser's developer tools JavaScript command prompt.
2. Copy the contents of `exporter.js` in and run it.
3. A file called `fig-export.ini` will be downloaded.

## Rendering a node from a document

1. Run `build_renderer.sh` in a UNIX terminal emulator.
2. Run something like `bin/renderer ~/Downloads/fig-export.ini <node name> output.png` where `<node name>` is the name of the node in Figma you want to export.

## Rough TODO list

- Text sizing is inaccurate.
- Shadow sizing/opacity is inaccurate.
- Allow the user to load custom fonts.
- Error reporting when a font is missing.
- Image scaling modes: fit, crop and tile.
- Image rotation.
- Image filters.
- Image nearest neighbour offset.
- Cache the results of `ComputeNodeDrawBounds`.
- Support clipping children in frames.
- Gradient dithering.
- Corner smoothing, see https://www.figma.com/blog/desperately-seeking-squircles/.
- Different stroke sizes on each edge of rectangles/frames.
- Clamping corner radii for polygons.
- Compute the boolean operation vector paths manually.
- Ellipse arc start/end and inner radius.
- Ellipse corner radii.
- Dashed strokes.
- Arrows.
- Open paths.
- For internal opaque strokes, prevent the fill's colour leaking outside into the antialiasing of the stroke.
- Supporting gradients, strokes, rotations, masking etc on text. The text will need to be converted to a vector path.
- Fix the kerning/tracking/hinting(?) of the text. Currently it looks off.
- Justified text paragraphs.
- Bidirectional text.
- Support Unicode escape sequences in text content.
- Shadows with fractional offsets.
- Shadow spread for rectangles/frames.
- Blend modes.
- Background blur effect.
- Support the isMaskOutline property.
- Bugs related to gradient positioning.
- When are there multiple vector paths?
- Don't duplicate image data in the exported file.
- Preserve instance/component information in the exporter, and update the renderer to manually handle it.
- Component set nodes.
- Code block nodes.
- Exporting auto-layout information.
- Line nodes.
- Shape with text nodes.
- Star nodes.
