#version 330 core
layout (location = 0) in vec2 position;

uniform vec3 sunPosition;
const float rayleigh = 2.0;
const float turbidity = 10.0;
const float mieCoefficient = 0.005;

out vec3 vSunDirection;
out float vSunfade;
out vec3 vBetaR;
out vec3 vBetaM;
out float vSunE;
out vec3 vSunColor;
out vec3 vAmbient;


/*
Atmospheric scattering based off of: https://www.shadertoy.com/view/XtBXDz
Author: valentingalea
*/

#define PI 3.14159265359

// scattering coefficients at sea level (m)
const vec3 betaR = vec3(5.5e-6, 13.0e-6, 22.4e-6); // Rayleigh 
const vec3 betaM = vec3(21e-6); // Mie

// scale height (m)
// thickness of the atmosphere if its density were uniform
const float hR = 7994.0; // Rayleigh
const float hM = 1200.0; // Mie

const float earth_radius = 6360e3; // (m)
const float atmosphere_radius = 6420e3; // (m)

const float sun_power = 30.0;

//These numbers can be lowered to acheive faster speed, but they are good enough here
const int num_samples = 16;
const int num_samples_light = 8;

struct ray_t {
	vec3 origin;
	vec3 direction;
};

bool isect_sphere(const in ray_t ray, inout float t0, inout float t1)
{
	vec3 rc = vec3(0.0) - ray.origin;
	float radius2 = atmosphere_radius * atmosphere_radius;
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

float HG(float costheta, float g) {
	const float k = 0.0795774715459; 
	return k*(1.0-g*g)/(pow(1.0+g*g-2.0*g*costheta, 1.5));
}

bool get_sun_light(
	const in ray_t ray,
	inout float optical_depthR,
	inout float optical_depthM
){
	float t0 = 0.0;
	float t1 = 0.0;
	isect_sphere(ray, t0, t1);

	float march_pos = 0.;
	float march_step = t1 / float(num_samples_light);

	for (int i = 0; i < num_samples_light; i++) {
		vec3 s =
			ray.origin +
			ray.direction * (march_pos + 0.5 * march_step);
		float height = length(s) - earth_radius;
		if (height < 0.0)
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
		ray, t0, t1)) {
		return vec3(0);
	}

	float march_step = t1 / float(num_samples);

	vec3 sun_dir = normalize(sunPosition);
	// cosine of angle between view and light directions
	float mu = dot(ray.direction, sun_dir);

	// Rayleigh and Mie phase functions
	float phaseR = rayleigh_phase_func(mu);
	float phaseM = HG(mu, 0.96);//0.76 more proper but this looks nice

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

/*
Mostly Preetham model stuff here
*/

const vec3 up = vec3( 0.0, 1.0, 0.0 );

// constants for atmospheric scattering
const float e = 2.71828182845904523536028747135266249775724709369995957;
const float pi = 3.141592653589793238462643383279502884197169;

// wavelength of used primaries, according to preetham
const vec3 lambda = vec3( 680E-9, 550E-9, 450E-9 );
// this pre-calcuation replaces older TotalRayleigh(vec3 lambda) function:
// (8.0 * pow(pi, 3.0) * pow(pow(n, 2.0) - 1.0, 2.0) * (6.0 + 3.0 * pn)) / (3.0 * N * pow(lambda, vec3(4.0)) * (6.0 - 7.0 * pn))
const vec3 totalRayleigh = vec3( 5.804542996261093E-6, 1.3562911419845635E-5, 3.0265902468824876E-5 );

// mie stuff
// K coefficient for the primaries
const float v = 4.0;
const vec3 K = vec3( 0.686, 0.678, 0.666 );
// MieConst = pi * pow( ( 2.0 * pi ) / lambda, vec3( v - 2.0 ) ) * K
const vec3 MieConst = vec3( 1.8399918514433978E14, 2.7798023919660528E14, 4.0790479543861094E14 );

// earth shadow hack
// cutoffAngle = pi / 1.95;
const float cutoffAngle = 1.6110731556870734;
const float steepness = 1.5;
const float EE = 1000.0;

float sunIntensity( float zenithAngleCos ) {
	zenithAngleCos = clamp( zenithAngleCos, -1.0, 1.0 );
	return EE * max( 0.0, 1.0 - pow( e, -( ( cutoffAngle - acos( zenithAngleCos ) ) / steepness ) ) );
}

vec3 totalMie( float T ) {
	float c = ( 0.2 * T ) * 10E-18;
	return 0.434 * c * MieConst;
}

void main() {
	gl_Position = vec4( position, -1.0, 1.0 );

	vSunDirection = normalize( sunPosition );

	vSunE = sunIntensity( dot( vSunDirection, up ) );

	vSunfade = 1.0 - clamp( 1.0 - exp( ( sunPosition.y / 450000.0 ) ), 0.0, 1.0 );

	float rayleighCoefficient = rayleigh - ( 1.0 * ( 1.0 - vSunfade ) );
	// extinction (absorbtion + out scattering)
	// rayleigh coefficients
	vBetaR = totalRayleigh * rayleighCoefficient;
	// mie coefficients
	vBetaM = totalMie( turbidity ) * mieCoefficient;

	ray_t ray = ray_t(vec3(0.0, earth_radius+1.0, 0.0), normalize(vSunDirection+vec3(0.01, 0.01, 0.0)));
	vSunColor = get_incident_light(ray);

	ray = ray_t(vec3(0.0, earth_radius+1.0, 0.0), normalize(vec3(0.4, 0.1, 0.0)));
	vAmbient = get_incident_light(ray);
}
	
