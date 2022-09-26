# fig-tools

Please note that this is very work in progress! If you have time to spare, please help out contributing features and resolving issues. Many are marked with `TODO` in the code.

## Exporting a document

1. Open your document in Figma, and open your web browser's developer tools JavaScript command prompt.
2. Copy the contents of `exporter.js` in and run it.
3. A file called `fig-export.ini` will be downloaded.

## Rendering a node from a document

1. Run `build_renderer.sh` in a UNIX terminal emulator.
2. Run something like `bin/renderer ~/Downloads/fig-export.ini <node name> output.png` where `<node name>` is the name of the node in Figma you want to export.
