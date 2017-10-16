#version 330 core
in vec3 ourColor;
in vec3 TexCoords;
in vec4 fragPos;

uniform sampler2D buff;
uniform sampler2D pong;
uniform int check;
uniform mat4 MVPM;
uniform mat4 LFMVPM;

out vec4 color;


int check_pos(vec2 x, float size) {
	return int(mod(floor(x.x), size) + mod(floor(x.y), size)*size);
}

void main()
{
	vec2 shift = vec2(floor(float(check)/4.0), mod(float(check), 4.0));	

	vec2 uv = floor(gl_FragCoord.xy/4.0);
	uv = uv/vec2(128.0);

	vec4 col = vec4(0.0);
	if (check_pos(gl_FragCoord.xy/1.0, 4.0)!=check&&true==true){
		//reprojection from http://john-chapman-graphics.blogspot.ca/2013/01/what-is-motion-blur-motion-pictures-are.html
		//look into running all this on cpu
		//discard;
		uv = gl_FragCoord.xy/512.0;
		vec2 uvd = uv-vec2(0.5);
		uvd *= 2.0;
		vec4 uvdir = (vec4(uvd, 1.0, 1.0));
		mat4 invmat = inverse(MVPM);
		vec4 worldPos = (inverse((MVPM))*uvdir);
		vec4 current = worldPos;
    vec4 previous = LFMVPM * current;
    previous.xyz /= previous.w;
    previous.xy = previous.xy * 0.5 + 0.5;
    vec2 blurVec = previous.xy - uv.xy;
		vec2 lookup = uv.xy+blurVec;
		float mip = 0.0;
		if (lookup.x<0.0||lookup.x>1.0||lookup.y<0.0||lookup.y>1.0) {
			lookup = clamp(lookup, 0.0, 1.0);
			lookup = uv.xy;
			col = texture(buff, lookup);
		} else {
			uv = gl_FragCoord.xy/512.0;
			col = texture(pong, lookup);
		}
	} else {
	col = texture(buff, uv.xy);
	}
	color.xyz = col.xyz;
	color.a = 1.0;
}
