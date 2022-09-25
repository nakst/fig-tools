// TODO Download the result as a file, rather than making the user copy and paste the message out.

var outText = "";
var objectIDAllocator = 1;

async function DumpPaint(paint) {
	var id = objectIDAllocator++;
	outText += "[" + id.toString() + "]\n";
	outText += "type=PAINT." + paint.type + "\n";
	outText += "visible=" + paint.visible + "\n";
	outText += "opacity=" + paint.opacity + "\n";
	outText += "blendMode=" + paint.blendMode + "\n";

	if (paint.type == "SOLID") {
		outText += "color=" + paint.color.r + "," + paint.color.g + "," + paint.color.b + "\n";
	} else if (paint.type == "IMAGE") {
		outText += "scaleMode=" + paint.scaleMode + "\n";
		outText += "preCropTransform00=" + paint.imageTransform[0][0] + "\n";
		outText += "preCropTransform01=" + paint.imageTransform[0][1] + "\n";
		outText += "preCropTransform02=" + paint.imageTransform[0][2] + "\n";
		outText += "preCropTransform10=" + paint.imageTransform[1][0] + "\n";
		outText += "preCropTransform11=" + paint.imageTransform[1][1] + "\n";
		outText += "preCropTransform12=" + paint.imageTransform[1][2] + "\n";
		outText += "tileScaling=" + paint.scalingFactor + "\n";
		outText += "rotationDegrees=" + paint.rotation + "\n"; // Ignored for crop.
		outText += "filterExposure=" + paint.filters.exposure + "\n";
		outText += "filterContrast=" + paint.filters.contrast + "\n";
		outText += "filterSaturation=" + paint.filters.saturation + "\n";
		outText += "filterTemperature=" + paint.filters.temperature + "\n";
		outText += "filterTint=" + paint.filters.tint + "\n";
		outText += "filterHighlights=" + paint.filters.highlights + "\n";
		outText += "filterShadows=" + paint.filters.shadows + "\n";

		// TODO Store image data separately, by hash.

		const image = figma.getImageByHash(paint.imageHash);
		const bytes = await image.getBytesAsync();
		outText += "imageData=";

		for (var i = 0; i < bytes.byteLength; i += 1) {
			var b = bytes[i];
			var b0 = (b >> 4) & 0xF;
			var b1 = b & 0xF;
			var hexChars = "0123456789ABCDEF";
			outText += hexChars[b1] + hexChars[b0];
		}

		outText += "\n";
	} else {
		outText += "gradientTransform00=" + paint.gradientTransform[0][0] + "\n";
		outText += "gradientTransform01=" + paint.gradientTransform[0][1] + "\n";
		outText += "gradientTransform02=" + paint.gradientTransform[0][2] + "\n";
		outText += "gradientTransform10=" + paint.gradientTransform[1][0] + "\n";
		outText += "gradientTransform11=" + paint.gradientTransform[1][1] + "\n";
		outText += "gradientTransform12=" + paint.gradientTransform[1][2] + "\n";

		var gradientStops = "";

		for (var i in paint.gradientStops) {
			gradientStops += paint.gradientStops[i].color.r + "," + paint.gradientStops[i].color.g + "," 
				+ paint.gradientStops[i].color.b + "," + paint.gradientStops[i].color.a + "," + paint.gradientStops[i].position + ";";
		}

		outText += "gradientStops=" + gradientStops + "\n";
	}

	return id;
}

async function DumpEffect(effect) {
	var id = objectIDAllocator++;
	outText += "[" + id.toString() + "]\n";
	outText += "type=EFFECT." + effect.type + "\n";
	outText += "visible=" + effect.visible + "\n";

	if (effect.type == "DROP_SHADOW" || effect.type == "INNER_SHADOW") {
		outText += "color=" + effect.color.r + "," + effect.color.g + "," + effect.color.b + "," + effect.color.a + "\n";
		outText += "offsetX=" + effect.offset.x + "\n";
		outText += "offsetY=" + effect.offset.y + "\n";
		outText += "radius=" + effect.radius + "\n";
		outText += "spread=" + effect.spread + "\n";
		outText += "blendMode=" + effect.blendMode + "\n";

		if (effect.type == "DROP_SHADOW") {
			outText += "showShadowBehindNode=" + effect.showShadowBehindNode + "\n";
		}
	} else if (effect.type == "LAYER_BLUR" || effect.type == "BACKGROUND_BLUR") {
		outText += "radius=" + effect.radius + "\n";
	}

	return id;
}

async function DumpTextSegment(segment) {
	var fillIDsArray = [];

	for (const fill in segment.fills) {
		fillIDsArray.push(await DumpPaint(segment.fills[fill]));
	}

	var id = objectIDAllocator++;

	outText += "[" + id.toString() + "]\n";
	outText += "type=TEXT_SEGMENT" + "\n";
	outText += "characters=" + JSON.stringify(segment.characters) + "\n";
	outText += "fontSize=" + segment.fontSize + "\n";
	outText += "fontFamily=" + segment.fontName.family + "\n";
	outText += "fontStyle=" + segment.fontName.style + "\n";
	outText += "textDecoration=" + segment.textDecoration + "\n";
	outText += "textCase=" + segment.textCase + "\n";
	outText += "fills=" + fillIDsArray.toString() + "\n";

	return id;
}

async function DumpNode(node, parentID) {
	var childrenIDsArray = [];
	var fillIDsArray = [];
	var strokeIDsArray = [];
	var effectIDsArray = [];
	var textRunIDsArray = [];
	var id = objectIDAllocator++;

	if ("children" in node) {
		for (const child in node.children) {
			childrenIDsArray.push(await DumpNode(node.children[child], parentID));
		}
	}

	if ("fills" in node) {
		for (const fill in node.fills) {
			fillIDsArray.push(await DumpPaint(node.fills[fill]));
		}
	}

	if ("strokes" in node) {
		for (const stroke in node.strokes) {
			strokeIDsArray.push(await DumpPaint(node.strokes[stroke]));
		}
	}

	if ("effects" in node) {
		for (const effect in node.effects) {
			effectIDsArray.push(await DumpEffect(node.effects[effect]));
		}
	}

	if (node.type == "TEXT") {
		var segments = node.getStyledTextSegments(["fontSize", "fontName", "textDecoration", 
			"textCase", "lineHeight", "letterSpacing", "fills", "textStyleId", "fillStyleId", "listOptions", "indentation", "hyperlink"]);

		for (const segment in segments) {
			textRunIDsArray.push(await DumpTextSegment(segments[segment]));
		}
	}

	outText += "[" + id.toString() + "]\n";
	outText += "type=" + node.type + "\n";
	outText += "name=" + node.name + "\n";
	outText += "parent=" + parentID.toString() + "\n";
	outText += "children=" + childrenIDsArray.toString() + "\n";
	outText += "isMask=" + node.isMask + "\n";

	if (node.type == "BOOLEAN_OPERATION") {
		outText += "booleanOperation=" + node.booleanOperation + "\n";
		outText += "vectorPath=" + node.fillGeometry[0].data + "\n";
	} else if (node.type == "POLYGON") {
		outText += "pointCount=" + node.pointCount + "\n";
	} else if (node.type == "VECTOR") {
		// TODO When are there more than one?
		outText += "vectorPath=" + node.vectorPaths[0].data + "\n";
	} else if (node.type == "ELLIPSE") {
		outText += "arcStartingAngle=" + node.arcData.startingAngle + "\n";
		outText += "arcEndingAngle=" + node.arcData.endingAngle + "\n";
		outText += "arcInnerRadius=" + node.arcData.innerRadius + "\n";
	} else if (node.type == "TEXT") {
		outText += "textAlignHorizontal=" + node.textAlignHorizontal + "\n";
		outText += "textAlignVertical=" + node.textAlignVertical + "\n";
		outText += "paragraphIndent=" + node.paragraphIndent + "\n";
		outText += "paragraphSpacing=" + node.paragraphSpacing + "\n";
		outText += "textRuns=" + textRunIDsArray.toString() + "\n";
	} else if (node.type == "COMPONENT" || node.type == "FRAME" || node.type == "RECTANGLE" || node.type == "INSTANCE") {
		if (node.type == "INSTANCE") outText += "scaleFactor=" + node.scaleFactor + "\n";
		outText += "topLeftRadius=" + node.topLeftRadius + "\n";
		outText += "topRightRadius=" + node.topRightRadius + "\n";
		outText += "bottomLeftRadius=" + node.bottomLeftRadius + "\n";
		outText += "bottomRightRadius=" + node.bottomRightRadius + "\n";
		outText += "strokeLeftWeight=" + node.strokeLeftWeight + "\n";
		outText += "strokeRightWeight=" + node.strokeRightWeight + "\n";
		outText += "strokeTopWeight=" + node.strokeTopWeight + "\n";
		outText += "strokeBottomWeight=" + node.strokeBottomWeight + "\n";
	}

	if (node.type == "BOOLEAN_OPERATION" || node.type == "COMPONENT" || node.type == "ELLIPSE" || node.type == "FRAME" 
			|| node.type == "INSTANCE" || node.type == "RECTANGLE" || node.type == "POLYGON" || node.type == "VECTOR" 
			|| node.type == "TEXT" || node.type == "GROUP") {
		outText += "visible=" + node.visible + "\n";
		outText += "opacity=" + node.opacity + "\n";
		outText += "blendMode=" + node.blendMode + "\n";
		outText += "effects=" + effectIDsArray.toString() + "\n";
		outText += "absoluteTransform00=" + node.absoluteTransform[0][0] + "\n";
		outText += "absoluteTransform01=" + node.absoluteTransform[0][1] + "\n";
		outText += "absoluteTransform02=" + node.absoluteTransform[0][2] + "\n";
		outText += "absoluteTransform10=" + node.absoluteTransform[1][0] + "\n";
		outText += "absoluteTransform11=" + node.absoluteTransform[1][1] + "\n";
		outText += "absoluteTransform12=" + node.absoluteTransform[1][2] + "\n";
		outText += "width=" + node.width + "\n";
		outText += "height=" + node.height + "\n";

		if (node.type != "GROUP") {
			outText += "fills=" + fillIDsArray.toString() + "\n";
			outText += "strokes=" + strokeIDsArray.toString() + "\n";
			outText += "strokeWeight=" + node.strokeWeight + "\n";
			outText += "strokeJoin=" + node.strokeJoin + "\n";
			outText += "strokeAlign=" + node.strokeAlign + "\n";
			outText += "dashPattern=" + node.dashPattern.toString() + "\n";
			outText += "strokeCap=" + node.strokeCap + "\n";
			outText += "strokeMiterLimit=" + node.strokeMiterLimit + "\n";

			if (node.type != "TEXT") {
				outText += "cornerRadius=" + node.cornerRadius + "\n";
				outText += "cornerSmoothing=" + node.cornerSmoothing + "\n";
			}
		}
	}

	outText += "\n";
	return id;
}

await DumpNode(figma.root, 0);
console.log(outText);
