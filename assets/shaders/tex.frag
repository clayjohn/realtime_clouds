#version 330 core
in vec3 ourColor;
in vec3 TexCoords;

uniform sampler2D fbo;

out vec4 color;

void main()
{
  color = texture(fbo, TexCoords.xy);  
}
