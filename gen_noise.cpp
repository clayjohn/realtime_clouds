#include <iostream>
#include <stdlib.h>
#include <math.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "include/glm/glm.hpp"
//#include "include/FastNoise/FastNoise.h"
#include "include/TileableVolumeNoise/TileableVolumeNoise.h"
#include "include/stb_image_write.h"

/*
 * g++ include/TileableVolumeNoise/tileableVolumeNoise.cpp gen_noise.cpp -o gen_noise -std=c++11
 */

// the remap function used in the shaders as described in Gpu Pro 7. It must match when using pre packed textures
float remap(float originalValue, float originalMin, float originalMax, float newMin, float newMax)
{
	return newMin + (((originalValue - originalMin) / (originalMax - originalMin)) * (newMax - newMin));
}

float smoothstep(float edge0, float edge1, float x) {
		float t = std::min(std::max((x - edge0) / (edge1 - edge0), 0.0f), 1.0f);
    return t * t * (3.0 - 2.0 * t);
}

int main() {
		//initialize perlin noise arrays for textures
		std::cout<<"Generating Perlin Noise for LUT's"<<std::endl;
//		FastNoise myNoise; // Create a FastNoise object
		std::cout<<"Generating weather Noise 512x512 RGB"<<std::endl;
		char *NoiseArray = new char[512*512*3];
		for (int i =0;i<512*512*3;i+=3) {
			
			glm::vec3 pos = glm::vec3(float((i/3)%512)/512.0, float((i/3)/512)/512.0, 0.051);
			glm::vec3 offset1 = glm::vec3(0.0, 0.0, 581.163);
			glm::vec3 offset2 = glm::vec3(0.0, 0.0, 1245.463);
			glm::vec3 offset3 = glm::vec3(0.0, 0.0, 2245.863);
			float perlinNoise = Tileable3dNoise::PerlinNoise(pos, 8, 3);
			float perlinNoise2 = Tileable3dNoise::PerlinNoise(pos+offset1, 8, 3);
			float perlinNoise3 = Tileable3dNoise::PerlinNoise(pos+offset2, 2, 3);
			//float perlinNoise4 = Tileable3dNoise::PerlinNoise(pos+offset3, 4, 3);
			perlinNoise3 = std::min(1.0, (smoothstep(0.45, 0.8, perlinNoise3)+smoothstep(0.25, 0.45, perlinNoise3)*0.5));
			NoiseArray[i] = char(perlinNoise*128.0+127.0);
			NoiseArray[i+1] = char(smoothstep(0.5, 0.7, perlinNoise2)*255.0);
			NoiseArray[i+2] = char(perlinNoise3*255.0);
		}
		stbi_write_bmp("assets/weather.bmp", 512, 512, 3, NoiseArray);
		delete NoiseArray;


	/*	
		std::cout<<"Generating Curl Noise 128x128 RGB"<<std::endl;
		char *curlNoiseArray = new char[128*128*3];
		myNoise.SetFrequency(0.08);
		for (int i =0;i<128*128*3;i+=3) {
			float *ttt = myNoise.GetDerivPerlin(100.013+(i/3)%128, 0.173+(i/3)/128, 0.55);
			float curl[3] = {ttt[3]-ttt[2], ttt[1]-ttt[3], ttt[2]-ttt[1]};
			delete ttt;
			curlNoiseArray[i] = char(curl[0]*50+128);
			curlNoiseArray[i+1] = char(curl[1]*50+128);
			curlNoiseArray[i+2] = char(curl[2]*50+128);
		}
		stbi_write_bmp("assets/curlnoise.bmp", 128, 128, 3, curlNoiseArray);
		delete curlNoiseArray;
*/
		
		//worley and perlin-worley are from github/sebh/TileableVolumeNoise
		//which is in turn based on noise described in 'real time rendering of volumetric cloudscapes for horizon zero dawn'
		std::cout<<"Generating Worley Noise 32x32x32 RGB"<<std::endl;
		char *worlNoiseArray = new char[32*32*32*3];
		for (int i=0;i<32*32*32*3;i+=3) {
			glm::vec3 pos = glm::vec3(float((i/3)%32)/32.0, float(((i/3)/32)%32)/32.0, float((i/3)/(32*32))/32.0);
			float cell0 = 1.0f - Tileable3dNoise::WorleyNoise(pos, 2);
			float cell1 = 1.0f - Tileable3dNoise::WorleyNoise(pos, 4);
			float cell2 = 1.0f - Tileable3dNoise::WorleyNoise(pos, 8);
			float cell3 = 1.0f - Tileable3dNoise::WorleyNoise(pos, 16);
			
			float cellFBM0 = cell0*0.5f + cell1*0.35f + cell2*0.15f;
			float cellFBM1 = cell1*0.5f + cell2*0.35f + cell3*0.15f;
			float cellFBM2 = cell2*0.75f + cell3*0.25f; // cellCount=4 -> worleyNoise4 is just noise due to sampling frequency=texel freque. So only take into account 2 frequenciM
			worlNoiseArray[i] = char(cellFBM0*255);
			worlNoiseArray[i+1] = char(cellFBM1*255);
			worlNoiseArray[i+2] = char(cellFBM2*255);
		}
		stbi_write_bmp("assets/worlnoise.bmp", 32*32, 32, 3, worlNoiseArray);
		delete worlNoiseArray;
		
		/*
		std::cout<<"Generating Perlin-Worley Noise 128x128x128 RGBA"<<std::endl;
		char *perlWorlNoiseArray = new char[128*128*128*4];
		for (int i=0;i<128*128*128*4;i+=4) {
			
				glm::vec3 pos = glm::vec3(float((i/4)%128)/128.0, float(((i/4)/128)%128)/128.0, float((i/4)/(128*128))/128.0);

				// Perlin FBM noise
				float perlinNoise = Tileable3dNoise::PerlinNoise(pos, 8, 3);
			
				const float worleyNoise00 = (1.0f - Tileable3dNoise::WorleyNoise(pos, 8));
				const float worleyNoise01 = (1.0f - Tileable3dNoise::WorleyNoise(pos, 32));
				const float worleyNoise02 = (1.0f - Tileable3dNoise::WorleyNoise(pos, 56));
					//const float worleyNoise3 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 80));
					//const float worleyNoise4 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 104));
					//const float worleyNoise5 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 128));	// half the frequency of texel, we should not go further (with cellCount = 32 and texture size = 64)

																														// PerlinWorley noise as described p.101 of GPU Pro 7
				float worleyFBM = worleyNoise00*0.625f + worleyNoise01*0.25f + worleyNoise02*0.125f;

				float PerlWorlNoise = remap(perlinNoise, 0.0, 1.0, worleyFBM, 1.0);
				

				//float worleyNoise0 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 4));
				//float worleyNoise1 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 8));
				float worleyNoise12 = (1.0f - Tileable3dNoise::WorleyNoise(pos, 16));
				//float worleyNoise3 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 32));
				float worleyNoise14 = (1.0f - Tileable3dNoise::WorleyNoise(pos, 64));

				// Three frequency of Worley FBM noise
				float worleyFBM0 = worleyNoise00*0.625f + worleyNoise12*0.25f + worleyNoise01*0.125f;
				float worleyFBM1 = worleyNoise12*0.625f + worleyNoise01*0.25f + worleyNoise14*0.125f;
				float worleyFBM2 = worleyNoise01*0.75f + worleyNoise14*0.25f; // cellCount=4 -> worleyNoise5 is just noise due to sampling frequency=texel frequency. So only take into account 2 frequencies for FBM


			perlWorlNoiseArray[i] = char(PerlWorlNoise*255);
			perlWorlNoiseArray[i+1] = char(worleyFBM0*255);
			perlWorlNoiseArray[i+2] = char(worleyFBM1*255);
			perlWorlNoiseArray[i+3] = char(worleyFBM2*255);
		}
		stbi_write_tga("assets/perlworlnoise.tga", 128*128, 128, 4, perlWorlNoiseArray);
		delete perlWorlNoiseArray;
*/

	return 0;
}
