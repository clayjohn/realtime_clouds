#version 330 core
in vec3 ourColor;
in vec3 TexCoords;

//uniform vec2 resolution;
const vec2 resolution = vec2(512.0);
uniform float time;
uniform mat4 MVPM;
uniform int check;
out vec4 color;

// ----------------------------------------------------------------------------
// Rayleigh and Mie scattering atmosphere system
//
// implementation of the techniques described here:
// http://www.scratchapixel.com/old/lessons/3d-advanced-lessons/simulating-the-colors-of-the-sky/atmospheric-scattering/
// ----------------------------------------------------------------------------

#define _in(T) const in T
#define _inout(T) inout T
#define _out(T) out T
#define _begin(type) type (
#define _end )
#define mul(a, b) (a) * (b)

#define PI 3.14159265359

// Shadertoy specific uniforms
#define u_res resolution
#define u_time time

struct ray_t {
	vec3 origin;
	vec3 direction;
};
#define BIAS 1e-4 // small offset to avoid self-intersections

struct sphere_t {
	vec3 origin;
	float radius;
	int material;
};

struct plane_t {
	vec3 direction;
	float distance;
	int material;
};

int check_pos(vec2 x, float size) {
	return int(mod(floor(x.x), size) + mod(floor(x.y), size)*size);
}

mat3 rotate_around_x(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(1, 0, 0, 0, _cos, -_sin, 0, _sin, _cos);
}


ray_t get_primary_ray(
	_in(vec3) cam_local_point,
	_inout(vec3) cam_origin,
	_inout(vec3) cam_look_at
){
	vec3 fwd = normalize(cam_look_at - cam_origin);
	vec3 up = vec3(0, 1, 0);
	vec3 right = cross(up, fwd);
	up = cross(fwd, right);

	ray_t r = _begin(ray_t)
		cam_origin,
		normalize(fwd + up * cam_local_point.y + right * cam_local_point.x)
		_end;
	return r;
}

bool isect_sphere(_in(ray_t) ray, _in(sphere_t) sphere, _inout(float) t0, _inout(float) t1)
{
	vec3 rc = sphere.origin - ray.origin;
	float radius2 = sphere.radius * sphere.radius;
	float tca = dot(rc, ray.direction);
	float d2 = dot(rc, rc) - tca * tca;
	if (d2 > radius2) return false;
	float thc = sqrt(radius2 - d2);
	t0 = tca - thc;
	t1 = tca + thc;

	return true;
}

// scattering coefficients at sea level (m)
const vec3 betaR = vec3(5.5e-6, 13.0e-6, 22.4e-6); // Rayleigh 
const vec3 betaM = vec3(21e-6); // Mie

// scale height (m)
// thickness of the atmosphere if its density were uniform
const float hR = 7994.0; // Rayleigh
const float hM = 1200.0; // Mie

float rayleigh_phase_func(float mu)
{
	return
			3. * (1. + mu*mu)
	/ //------------------------
				(16. * PI);
}

// Henyey-Greenstein phase function factor [-1, 1]
// represents the average cosine of the scattered directions
// 0 is isotropic scattering
// > 1 is forward scattering, < 1 is backwards
const float g = 0.76;
float henyey_greenstein_phase_func(float mu)
{
	return
						(1. - g*g)
	/ //---------------------------------------------
		((4. + PI) * pow(1. + g*g - 2.*g*mu, 1.5));
}

// Schlick Phase Function factor
// Pharr and  Humphreys [2004] equivalence to g above
const float k = 1.55*g - 0.55 * (g*g*g);
float schlick_phase_func(float mu)
{
	return
					(1. - k*k)
	/ //-------------------------------------------
		(4. * PI * (1. + k*mu) * (1. + k*mu));
}

const float earth_radius = 6360e3; // (m)
const float atmosphere_radius = 6420e3; // (m)

vec3 sun_dir = vec3(0, 1, 0);
const float sun_power = 20.0;

const sphere_t atmosphere = _begin(sphere_t)
	vec3(0, 0, 0), atmosphere_radius, 0
_end;

const int num_samples = 16;
const int num_samples_light = 8;

bool get_sun_light(
	_in(ray_t) ray,
	_inout(float) optical_depthR,
	_inout(float) optical_depthM
){
	float t0 = 0.0;
	float t1 = 0.0;
	isect_sphere(ray, atmosphere, t0, t1);

	float march_pos = 0.;
	float march_step = t1 / float(num_samples_light);

	for (int i = 0; i < num_samples_light; i++) {
		vec3 s =
			ray.origin +
			ray.direction * (march_pos + 0.5 * march_step);
		float height = length(s) - earth_radius;
		if (height < 0.)
			return false;

		optical_depthR += exp(-height / hR) * march_step;
		optical_depthM += exp(-height / hM) * march_step;

		march_pos += march_step;
	}

	return true;
}

vec3 get_incident_light(_in(ray_t) ray)
{
	// "pierce" the atmosphere with the viewing ray
	float t0 = 0.0;
	float t1 = 0.0;
	if (!isect_sphere(
		ray, atmosphere, t0, t1)) {
		return vec3(0);
	}

	float march_step = t1 / float(num_samples);

	// cosine of angle between view and light directions
	float mu = dot(ray.direction, sun_dir);

	// Rayleigh and Mie phase functions
	// A black box indicating how light is interacting with the material
	// Similar to BRDF except
	// * it usually considers a single angle
	//   (the phase angle between 2 directions)
	// * integrates to 1 over the entire sphere of directions
	float phaseR = rayleigh_phase_func(mu);
	float phaseM =
#if 1
		henyey_greenstein_phase_func(mu);
#else
		schlick_phase_func(mu);
#endif

	// optical depth (or "average density")
	// represents the accumulated extinction coefficients
	// along the path, multiplied by the length of that path
	float optical_depthR = 0.;
	float optical_depthM = 0.;

	vec3 sumR = vec3(0);
	vec3 sumM = vec3(0);
	float march_pos = 0.;

	for (int i = 0; i < num_samples; i++) {
		vec3 s =
			ray.origin +
			ray.direction * (march_pos + 0.5 * march_step);
		float height = length(s) - earth_radius;

		// integrate the height scale
		float hr = exp(-height / hR) * march_step;
		float hm = exp(-height / hM) * march_step;
		optical_depthR += hr;
		optical_depthM += hm;

		// gather the sunlight
		ray_t light_ray = _begin(ray_t)
			s,
			sun_dir
		_end;
		float optical_depth_lightR = 0.;
		float optical_depth_lightM = 0.;
		bool overground = get_sun_light(
			light_ray,
			optical_depth_lightR,
			optical_depth_lightM);

		if (overground) {
			vec3 tau =
				betaR * (optical_depthR + optical_depth_lightR) +
				betaM * 1.1 * (optical_depthM + optical_depth_lightM);
			vec3 attenuation = exp(-tau);

			sumR += hr * attenuation;
			sumM += hm * attenuation;
		}

		march_pos += march_step;
	}

	return
		sun_power *
		(sumR * phaseR * betaR +
		sumM * phaseM * betaM);
}

vec3 ACESFilm( vec3 x )
{
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0, 1.0);
}

void main()
{
	vec2 aspect_ratio = vec2(u_res.x / u_res.y, 1);
	float fov = tan(radians(45.0));
	vec2 point_ndc = gl_FragCoord.xy / u_res.xy;
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, 1.0);
	
	vec4 worldPos = (inverse((MVPM))*vec4(point_cam, 1.0));
	worldPos.xyz /= worldPos.w;
	point_cam = normalize(worldPos.xyz);


	if (check_pos(gl_FragCoord.xy/4.0, 4.0)!=check){
		//reprojection from http://john-chapman-graphics.blogspot.ca/2013/01/what-is-motion-blur-motion-pictures-are.html
		//look into running all this on cpu
		discard;
		/*
		vec4 current = vec4((2.0 * point_ndc - 1.0) * aspect_ratio * fov, 1.0, 1.0);
    current = inverse(MVPM) * current;
    vec4 previous = LFMVPM * current;
    previous.xyz /= previous.w;
    previous.xy = previous.xy * 0.5 + 0.5;
    vec2 blurVec = previous.xy - TexCoords.xy;
		vec2 lookup = TexCoords.xy+blurVec;
		float mip = 0.0;
		if (lookup.x<0.0||lookup.x>1.0||lookup.y<0.0||lookup.y>1.0) {
			lookup = clamp(lookup, 0.0, 1.0);
			lookup = TexCoords.xy;
			mip = 1.0;
		}
		color = texture(lastFrame, lookup, mip);
		*/
	} else {

		vec3 col = vec3(0);
		if (point_cam.y>0.0) {
	// sun
			mat3 rot = rotate_around_x(-abs(sin(u_time / 20.)) * 90.);
			sun_dir *= rot;

      vec3 eye = vec3 (0, earth_radius + 1., 0);
      vec3 look_at = vec3 (0, earth_radius + 1.0, -1);
			
      //ray_t ray = get_primary_ray(point_cam, eye, look_at);
			ray_t ray = ray_t(vec3(0.0, earth_radius+1.0, 0.0), point_cam);
      col = get_incident_light(ray);
		}
		col = ACESFilm(col);
		color = vec4(col, 1);
	}
}
