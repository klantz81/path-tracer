#version 330

in vec2 textureCoord;
// in vec4 colorOut;
uniform sampler2D textureSample;

out vec4 fragColor;

void main (void) {
	fragColor = texture(textureSample, textureCoord);
// 	fragColor = colorOut*texture(textureSample, textureCoord);
}