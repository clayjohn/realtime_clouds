#version 330 core
in vec3 position;
in vec3 color;
in vec3 texCoords;

uniform float size;
uniform float pos;

out vec3 ourColor;
out vec3 TexCoords;

uniform mat4 MVPM;

void main()
{

    gl_Position = vec4(position.xy, 0.0, 1.0);
    ourColor = color;
    TexCoords = texCoords;
}
