#version 330

layout(location = 0) in vec4 position;
layout(location = 1) in vec4 color;

smooth out vec4 theColor;

uniform vec3 offset;
uniform mat4 perspectiveMatrix;
uniform mat4 viewMatrix;
uniform mat4 objMatrix;

void main()
{
	vec4 cameraPos = position + vec4(offset.x, offset.y, offset.z, 0.0);

	gl_Position = perspectiveMatrix * viewMatrix * objMatrix * cameraPos;
	theColor = color;
}

