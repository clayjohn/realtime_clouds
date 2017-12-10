#version 330 core

uniform vec2 resolution;
uniform float time;
uniform mat4 MVPM;

out vec4 color;

/*
Atmospheric scattering based off of: 
*/

#define PI 3.14159265359
#define BIAS 1e-4 // small offset to avoid self-intersections

// scattering coefficients at sea level (m)
const vec3 betaR = vec3(5.5e-6, 13.0e-6, 22.4e-6); // Rayleigh 
const vec3 betaM = vec3(21e-6); // Mie

// scale height (m)
// thickness of the atmosphere if its density were uniform
const float hR = 7994.0; // Rayleigh
const float hM = 1200.0; // Mie

const float earth_radius = 6360e3; // (m)
const float atmosphere_radius = 6420e3; // (m)

vec3 sun_dir = vec3(0, 1, 0);
const float sun_power = 30.0;


const int num_samples = 16;
const int num_samples_light = 8;



struct ray_t {
	vec3 origin;
	vec3 direction;
};

struct sphere_t {
	vec3 origin;
	float radius;
};

const sphere_t atmosphere = sphere_t(vec3(0, 0, 0), atmosphere_radius);


mat3 rotate_around_x(const in float angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(1, 0, 0, 0, _cos, -_sin, 0, _sin, _cos);
}

bool isect_sphere(const in ray_t ray, const in sphere_t sphere, inout float t0, inout float t1)
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
float henyey_greenstein_phase_func(float mu)
{
const float g = 0.76;
	return
						(1. - g*g)
	/ //---------------------------------------------
		((4. + PI) * pow(1. + g*g - 2.*g*mu, 1.5));
}

bool get_sun_light(
	const in ray_t ray,
	inout float optical_depthR,
	inout float optical_depthM
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

vec3 get_incident_light(const in ray_t ray)
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
	henyey_greenstein_phase_func(mu);

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
		ray_t light_ray = ray_t(s, sun_dir);
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

vec3 U2Tone(vec3 x) {
	const float A = 0.15;
	const float B = 0.50;
	const float C = 0.10;
	const float D = 0.20;
	const float E = 0.02;
	const float F = 0.30;

   return ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}

void main()
{
	vec2 aspect_ratio = vec2(resolution.x / resolution.y, 1);
	vec2 point_ndc = gl_FragCoord.xy / resolution.xy;
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio, 1.0);
	
	vec4 worldPos = (inverse((MVPM))*vec4(point_cam, 1.0));
	worldPos.xyz /= worldPos.w;
	point_cam = normalize(worldPos.xyz);

		vec3 col = vec3(0.4);
		if (point_cam.y>-0.05) {
			// sun
			mat3 rot = rotate_around_x(-abs(sin(time / 20.)) * 90.);
			sun_dir *= rot;
			
			ray_t ray = ray_t(vec3(0.0, earth_radius+1.0, 0.0), point_cam);
      col = get_incident_light(ray);
		}
		col = U2Tone(col);
		col /= U2Tone(vec3(2.5));
		col = sqrt(col);
		color = vec4(col, 1);
}
