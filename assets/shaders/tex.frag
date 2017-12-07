#version 330 core
in vec2 TexCoords;

uniform sampler2D fbo;

out vec4 color;

void main()
{
  color = texture(fbo, TexCoords, 0.0);  
}
