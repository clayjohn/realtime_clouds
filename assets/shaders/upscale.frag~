#version 330 core
in vec3 ourColor;
in vec3 TexCoords;

uniform sampler2D buff;
uniform int check;

out vec4 color;


int check_pos(vec2 x, float size) {
	return int(mod(floor(x.x), size) + mod(floor(x.y), size)*size);
}

void main()
{
 	//if (check_pos(gl_FragCoord.xy/4.0, 4.0)!=check){
		//discard;
	//}
	//compare pixel color to mipped pixel color
	vec4 pix = texture(buff, TexCoords.xy);
	//vec3 mip = texture(buff, TexCoords.xy, 1.0).xyz;
	//float diff = abs(length(pix.xyz)-length(mip))*pix.w;
	color.a = 1.0;
	color.xyz = pix.xyz;
	//if (diff>1.0) {
		//color.xyz = mip;
	//}
}
