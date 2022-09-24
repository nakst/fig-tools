# fig-tools

Please note that this is very work in progress! If you have time to spare, please help out contributing features and resolving issues. Many are marked with `TODO` in the code.

## Exporting a document

1. Open your document in Figma, and open your web browser's developer tools JavaScript command prompt.
2. Copy the contents of `exporter.js` in and run it.
3. Copy the output message (you might need to expand the message first, depending on your browser) and save it to a file.

## Rendering a node from a document

1. Run `build_renderer.sh` in a UNIX terminal emulator.
2. Run `bin/renderer <input file> <node name> <output file>` where `<input file>` is the path to the file you exported, `<node name>` is the name of the node in Figma, and `<output file>` is the PNG output file.
