#version 330 core
in vec3 ourColor;
in vec3 TexCoords;
in vec4 fragPos;

//TODO update cloud shaping
//TODO fix up weather texture
//TODO add more flexibility to parameters
	//rain
	//

uniform sampler3D perlworl;
uniform sampler3D worl;
uniform sampler2D curl;
uniform sampler2D weather;
uniform sampler2D atmosphere;

uniform int check;
uniform mat4 MVPM; 
uniform mat4 LFMVPM;
uniform mat4 VM;
uniform float aspect;
uniform float time;
uniform int camera_dirty;

out vec4 color;

const float g_radius = 600000.0; //ground radius
const float sky_b_radius = 601000.0;//bottom of cloud layer
const float sky_t_radius = 602300.0;//top of cloud layer
const float c_radius = 6008400.0; //2d noise layer

const vec3 RANDOM_VECTORS[6] = vec3[6]
(
	vec3( 0.38051305f,  0.92453449f, -0.02111345f),
	vec3(-0.50625799f, -0.03590792f, -0.86163418f),
	vec3(-0.32509218f, -0.94557439f,  0.01428793f),
	vec3( 0.09026238f, -0.27376545f,  0.95755165f),
	vec3( 0.28128598f,  0.42443639f, -0.86065785f),
	vec3(-0.16852403f,  0.14748697f,  0.97460106f)
	);


vec3 U2Tone(vec3 x) {
	const float A = 0.15;
	const float B = 0.50;
	const float C = 0.10;
	const float D = 0.20;
	const float E = 0.02;
	const float F = 0.30;

   return ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}


mat3 rotate_around_x(float angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(1, 0, 0, 0, _cos, -_sin, 0, _sin, _cos);
}

vec3 getSunDirection() {
	vec3 sun_dir = vec3(0.0, 1.0, 0.0);

	mat3 rot = rotate_around_x(-abs(sin(time / 20.)) * 90.);
	sun_dir *= rot;
	return sun_dir;
}

vec3 getSunColor() {
	vec3 dir = getSunDirection();
	return mix(vec3(1.0, 0.3, 0.0), vec3(0.99), smoothstep(0.05, 0.2, dir.y));
}

vec3 getSkyColor() {
	return texture(atmosphere, vec2(0.5, 0.9), 7).xyz;
}

int check_pos(vec2 x, float size) {
	return int(mod(floor(x.x), size) + mod(floor(x.y), size)*size);
}

float hash(vec3 p) {
	return texture(worl, p*0.01).x;
}


// fractional value for sample position in the cloud layer
float GetHeightFractionForPoint(vec3 inPosition, vec2 inCloudMinMax)
{ // get global fractional position in cloud zone
	float height_fraction = (inPosition.y - inCloudMinMax.x ) / (inCloudMinMax.y - inCloudMinMax.x); 
	return clamp(height_fraction, 0.0, 1.0);
}

vec4 mixGradients( float cloudType){

	const vec4 STRATUS_GRADIENT = vec4(0.02f, 0.05f, 0.09f, 0.11f);
	const vec4 STRATOCUMULUS_GRADIENT = vec4(0.02f, 0.2f, 0.48f, 0.625f);
	const vec4 CUMULUS_GRADIENT = vec4(0.01f, 0.0625f, 0.78f, 1.0f); // these fractions would need to be altered if cumulonimbus are added to the same pass
	float stratus = 1.0f - clamp(cloudType * 2.0f, 0.0, 1.0);
	float stratocumulus = 1.0f - abs(cloudType - 0.5f) * 2.0f;
	float cumulus = clamp(cloudType - 0.5f, 0.0, 1.0) * 2.0f;
	return STRATUS_GRADIENT * stratus + STRATOCUMULUS_GRADIENT * stratocumulus + CUMULUS_GRADIENT * cumulus;
}

float densityHeightGradient(float heightFrac, float cloudType) {
	vec4 cloudGradient = mixGradients(cloudType);
	return smoothstep(cloudGradient.x, cloudGradient.y, heightFrac) - smoothstep(cloudGradient.z, cloudGradient.w, heightFrac);
}

float intersectSphere(vec3 pos, vec3 dir, float r) {
    float a = dot(dir, dir);
    float b = 2.0 * dot(dir, pos);
    float c = dot(pos, pos) - (r * r);
		float d = sqrt((b*b) - 4.0*a*c);
		float p = -b - d;
		float p2 = -b + d;
    return max(p, p2)/(2.0*a);
}

// Utility function that maps a value from one range to another. 

// the remap function used in the shaders as described in Gpu Pro 7. It must match when using pre packed textures
float remap(float originalValue, float originalMin, float originalMax, float newMin, float newMax)
{
	return newMin + (((originalValue - originalMin) / (originalMax - originalMin)) * (newMax - newMin));
}

//From Thomas Shander, retreived from: https://www.shadertoy.com/view/4sjBDG
float numericalMieFit(float costh)
{
	float s = sign(costh);
	costh = sqrt(abs(costh))*s;
    // This function was optimized to minimize (delta*delta)/reference in order to capture
    // the low intensity behavior.
    float bestParams[10];
    bestParams[0]=9.805233e-06;
    bestParams[1]=-6.500000e+01;
    bestParams[2]=-5.500000e+01;
    bestParams[3]=8.194068e-01;
    bestParams[4]=1.388198e-01;
    bestParams[5]=-8.370334e+01;
    bestParams[6]=7.810083e+00;
    bestParams[7]=2.054747e-03;
    bestParams[8]=2.600563e-02;
    bestParams[9]=-4.552125e-12;
    
    float p1 = costh + bestParams[3];
    vec4 expValues = exp(vec4(bestParams[1] *costh+bestParams[2], bestParams[5] *p1*p1, bestParams[6] *costh, bestParams[9] *costh));
    vec4 expValWeight= vec4(bestParams[0], bestParams[4], bestParams[7], bestParams[8]);
    return dot(expValues, expValWeight);
}



float HG(vec3 inv, vec3 outv, float g) {
	float costheta = dot(inv, outv);
	const float k = 0.0795774715459; 
	return k*(1.0-g*g)/(pow(1.0+g*g-2.0*g*costheta, 1.5));
}

float density(vec3 p,vec3 weather, bool hq, float LOD) {
	p.x += time*10.0;
	float height_fraction = GetHeightFractionForPoint(p, vec2(float(sky_b_radius), float(sky_t_radius)));
	vec4 n = textureLod(perlworl, p*0.0003, LOD);
	float fbm = n.g*0.625+n.b*0.25+n.a*0.125;
	weather.x = smoothstep(0.6, 1.2, weather.x);
	float g = densityHeightGradient(height_fraction, weather.z);
	float base_cloud = remap(n.r, -(1.0-fbm), 1.0, 0.0, 1.0);
	float cloud_coverage = weather.x;
	base_cloud = remap(base_cloud*g, 1.0-cloud_coverage, 1.0, 0.0, 1.0); 
	base_cloud *= cloud_coverage;
	if (hq) {
		vec2 whisp = texture(curl, p.xy*0.001).xy;
		p.xy += whisp*300.0*(1.0-height_fraction);
		vec3 hn = textureLod(worl, p*0.002, LOD).xyz;
		float hfbm = hn.r*0.625+hn.g*0.25+hn.b*0.125;
		hfbm = mix(hfbm, 1.0-hfbm, clamp(height_fraction*10.0, 0.0, 1.0));
		base_cloud = remap(base_cloud, hfbm*0.2, 1.0, 0.0, 1.0);
	}
	return clamp(base_cloud, 0.0, 1.0);
}

vec4 march(vec3 pos, vec3 end, vec3 dir, int depth) {
	float T = 1.0;
	float alpha = 0.0;
	vec3 p = pos;
	float ss = length(dir);
	const float t_dist = sky_t_radius-sky_b_radius;
	float lss = t_dist/float(depth);
	vec3 ldir = getSunDirection()*ss;
	vec3 L = vec3(0.0);//getSkyColor();
	int count=0;
	float t = 1.0;
	float phase = max(HG(normalize(ldir), normalize(dir), 0.6), HG(normalize(ldir), normalize(dir), (0.99-1.5*normalize(ldir).y)));
	for (int i=0;i<depth;i++) {
		p += dir;
		float height_fraction = GetHeightFractionForPoint(p, vec2(float(sky_b_radius), float(sky_t_radius)));
		float weather_scale = 0.0001;
		vec3 weather_sample = texture(weather, p.xz*weather_scale).xyz;
		t = density(p, weather_sample, true, 1.0);
		const float ldt = 0.5;
		float dt = exp(-ldt*t*ss);//, exp(-ldt*0.25*t*ss)*0.8);
		T *= dt;		
		vec3 lp = p;
		float lT = 1.0;
		const float ld = 1.0;
		float lt = 1.0;
		float ncd = 0.0;
		float cd = 0.0;
		if (t>0.0) {
			for (int j=0;j<6;j++) {
				lp += (ldir+(RANDOM_VECTORS[j]*float(j+1))*lss);
				vec3 lweather = texture(weather, lp.xz*weather_scale).xyz;
				lt = density(lp, lweather, false, float(j));
				cd += lt;
				ncd += (lt * (1.0-(cd*(1.0/(lss*6.0)))));
				}
				lp += ldir*12.0;
				vec3 lweather = texture(weather, lp.xz*weather_scale).xyz;
				lt = density(lp, lweather, false, 7.0);
				cd += lt;
				ncd += (lt * (1.0-(cd*(1.0/(lss*6.0)))));

		}
		float beers = max(exp(-ld*ncd*lss), exp(-ld*0.25*ncd*lss)*0.7);
		float powshug = 1.0-exp(-ld*ncd*lss*2.0);

		vec3 ambient = 5.0*getSkyColor()*mix(0.25, 1.0, height_fraction);
		vec3 sunC = getSunColor()*30.0;
		L += (ambient+sunC*beers*powshug*2.0*phase)*(t)*T*ss;		
		alpha += (1.0-dt)*(1.0-alpha);
	}
	L = U2Tone(L);
	L /= U2Tone(vec3(40.0));
	L = sqrt(L);

	return vec4(L, alpha);
}


void main()
{
	vec2 shift = vec2(floor(float(check)/4.0), mod(float(check), 4.0));	
	//shift = vec2(0.0);
	vec2 uv = (gl_FragCoord.xy*4.0+shift.yx)/vec2(512.0);
	uv = uv-vec2(0.5);
	uv *= 2.0;
	uv.x *= aspect;
	vec4 uvdir = (vec4(uv.xy, 1.0, 1.0));
	mat4 invmat = inverse(MVPM);
	vec4 worldPos = (inverse((MVPM))*uvdir);
	//worldPos.xyz /= worldPos.w;
	vec3 camPos = vec3(invmat[3]);
	vec3 dir = normalize(worldPos.xyz/worldPos.w);

	vec4 col = vec4(0.0);
	vec3 background = textureLod(atmosphere, TexCoords.xy, 0.0).xyz;
	//vec3 background = vec3(dot(normalize(dir), normalize(getSunDirection())));
	if (dir.y>0.0) {
		vec3 start = camPos+vec3(0.0, g_radius, 0.0)+dir*intersectSphere(camPos+vec3(0.0, g_radius, 0.0), dir, sky_b_radius);
		vec3 end = camPos+vec3(0.0, g_radius, 0.0)+dir*intersectSphere(camPos+vec3(0.0, g_radius, 0.0), dir, sky_t_radius);
		const float t_dist = sky_t_radius-sky_b_radius;
		float shelldist = (length(end-start));
		vec4 volume;
		int steps = int(mix(96.0, 54.0, dot(dir, vec3(0.0, 1.0, 0.0))));
		float dmod = smoothstep(0.0, 1.0, (shelldist/t_dist)/14.0);
		float s_dist = mix(t_dist, t_dist*4.0, dmod)/float(steps);
		vec3 raystep = dir*s_dist;//(shelldist/float(steps));
		volume = march(start, end, raystep, steps);
		col = vec4(background*(1.0-volume.a)+volume.xyz*volume.a, 1.0);
		if (volume.a>1.0) {col = vec4(1.0, 0.0, 0.0, 1.0);}
	} else {
		col = vec4(vec3(0.4), 1.0);
		//col.xyz = texture(curl, uv).yyy;
	}
	
	color = col;
}
