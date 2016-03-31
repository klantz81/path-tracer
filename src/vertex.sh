#version 330

in vec3 vertex;
in vec3 texCoord;
// uniform vec4 color;
uniform mat4 Projection;
uniform mat4 View;
uniform mat4 Model;

out vec2 textureCoord;
// out vec4 colorOut;

void main() {
	gl_Position = Projection * View * Model * vec4(vertex, 1.0);
	textureCoord = texCoord.xy;
// 	colorOut = color;
}