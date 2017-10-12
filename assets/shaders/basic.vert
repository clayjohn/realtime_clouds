#version 330 core
in vec3 position;
in vec3 color;
in vec3 texCoords;

out vec3 ourColor;
out vec3 TexCoords;

uniform mat4 MVPM;

void main()
{

    gl_Position = MVPM * vec4(position, 1.0f);
    ourColor = color;
    TexCoords = texCoords;
}
