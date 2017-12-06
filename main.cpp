#include <iostream>
#include <stdlib.h> 

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// Other includes
#include "include/shader.h"
#include "include/camera.h"

// glm
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#define STB_IMAGE_IMPLEMENTATION
#include "include/stb_image.h"


//callbacks
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
float remap(float originalValue, float originalMin, float originalMax, float newMin, float newMax);

//function prototypes
void Do_Movement();

// Camera
Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
bool keys[1024];
GLfloat lastX = 400, lastY = 300;
bool firstMouse = true;

GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;
GLuint frames = 0;
GLfloat timePassed = 0.0f;
GLfloat startTime = 0.0f;
glm::mat4 MVPM;
glm::mat4 LFMVPM;
int camera_dirty = 0;


// Window dimensions
const GLuint WIDTH = 600, HEIGHT = 600;
GLfloat ASPECT = float(WIDTH)/float(HEIGHT);

// The MAIN function, from here we start the application and run the game loop
int main()
{
    // Init GLFW
    glfwInit();
    // Set all the required options for GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
//    glfwWindowHint(GLFW_SAMPLES, 4);

    // Create a GLFWwindow object that we can use for GLFW's functions
    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Realtime Clouds", NULL, NULL);
    glfwMakeContextCurrent(window);
				
    // Set the required callback functions
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    //glfw options
    glfwSetWindowPos(window, 200, 17);//so you can see frame rate
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
		glfwSwapInterval(0);//turn off vsync
    // Set this to true so GLEW knows to use a modern approach to retrieving function pointers and extensions
    glewExperimental = GL_TRUE;
    // Initialize GLEW to setup the OpenGL Function pointers
    glewInit();

    // Define the viewport dimensions
    glViewport(0, 0, WIDTH, HEIGHT);

    //openGL options
    //glEnable(GL_MULTISAMPLE);
    //glEnable(GL_DEPTH_TEST);
    //glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    //glEnable(GL_PROGRAM_POINT_SIZE); 
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_STENCIL_TEST);

    // Build and compile our shader program
    Shader ourShader("sky.vert", "sky.frag");
		Shader postShader("tex.vert", "tex.frag");
		Shader upscaleShader("upscale.vert", "upscale.frag");
		Shader atmoShader("atmo.vert", "atmo.frag");

    // Set up vertex data (and buffer(s)) and attribute pointers
    GLfloat vertices[] = {
        // Positions         // Colors            //texcoords
        1.0f, -1.0f,  0.0f,   1.0f, 1.0f, 0.0f,    1.0f,  0.0f, 0.0f,
       -1.0f, -1.0f,  0.0f,   0.0f, 1.0f, 1.0f,    0.0f,  0.0f, 0.0f,
       -1.0f,  1.0f,  0.0f,   1.0f, 0.0f, 1.0f,    0.0f,  1.0f, 0.0f,
        1.0f, -1.0f,  0.0f,   1.0f, 1.0f, 0.0f,    1.0f,  0.0f, 0.0f,
       -1.0f,  1.0f,  0.0f,   1.0f, 0.0f, 1.0f,    0.0f,  1.0f, 0.0f,
        1.0f,  1.0f,  0.0f,   1.0f, 0.0f, 1.0f,    1.0f,  1.0f, 0.0f 
    };
	

    GLuint VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    // Position attribute
    GLint posAttrib = glGetAttribLocation(ourShader.Program, "position");
    glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(posAttrib);
    // Color attribute
    GLint colAttrib = glGetAttribLocation(ourShader.Program, "color");
    glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(colAttrib);

    GLint texAttrib = glGetAttribLocation(ourShader.Program, "texCoords");
    glVertexAttribPointer(texAttrib, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (GLvoid*)(6 * sizeof(GLfloat)));
    glEnableVertexAttribArray(texAttrib);

		GLuint fbo, fbotex;

		
    glGenFramebuffers(1, &fbo);
    glGenTextures(1, &fbotex);
    glBindTexture(GL_TEXTURE_2D, fbotex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 512, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbotex, 0);

		GLuint pongfbo, pongfbotex;

		
    glGenFramebuffers(1, &pongfbo);
    glGenTextures(1, &pongfbotex);
    glBindTexture(GL_TEXTURE_2D, pongfbotex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 512, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		//set texture to repeat
    glBindTexture(GL_TEXTURE_2D, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, pongfbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pongfbotex, 0);

		GLuint buffer1, buffertex1;
		
    glGenFramebuffers(1, &buffer1);
    glGenTextures(1, &buffertex1);
    glBindTexture(GL_TEXTURE_2D, buffertex1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 128, 128, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		//glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, buffer1);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, buffertex1, 0);

		GLuint skyfbo, skytex;
		
    glGenFramebuffers(1, &skyfbo);
    glGenTextures(1, &skytex);
    glBindTexture(GL_TEXTURE_2D, skytex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 512, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, skyfbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, skytex, 0);


		glBindFramebuffer(GL_FRAMEBUFFER, 0);

				//setup noise textures
		GLuint curltex, worltex, perlworltex, weathertex;


		//stbi_set_flip_vertically_on_load(true);
		int x, y, n;
		unsigned char *curlNoiseArray = stbi_load("assets/curlnoise.bmp", &x, &y, &n, 0);

		glGenTextures(1, &curltex);
    glBindTexture(GL_TEXTURE_2D, curltex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 128, 0, GL_RGB, GL_UNSIGNED_BYTE, curlNoiseArray);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
		stbi_image_free(curlNoiseArray);

		unsigned char *weatherNoiseArray = stbi_load("assets/weather.bmp", &x, &y, &n, 0);

		glGenTextures(1, &weathertex);
    glBindTexture(GL_TEXTURE_2D, weathertex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 512, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, weatherNoiseArray);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
		stbi_image_free(weatherNoiseArray);

		unsigned char *worlNoiseArray = stbi_load("assets/worlnoise.bmp", &x, &y, &n, 0);
		glGenTextures(1, &worltex);
    glBindTexture(GL_TEXTURE_3D, worltex);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB, 32, 32, 32, 0, GL_RGB, GL_UNSIGNED_BYTE, worlNoiseArray);
    //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glGenerateMipmap(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, 0);
		stbi_image_free(worlNoiseArray);

		unsigned char *perlWorlNoiseArray = stbi_load("assets/perlworlnoise.tga", &x, &y, &n, 4);

		glGenTextures(1, &perlworltex);
    glBindTexture(GL_TEXTURE_3D, perlworltex);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, 128, 128, 128, 0, GL_RGBA, GL_UNSIGNED_BYTE, perlWorlNoiseArray);
    //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glGenerateMipmap(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, 0);
		stbi_image_free(perlWorlNoiseArray);


        //set up Model-View-Projection matrix
    //this way you only update when camera moves
    glm::mat4 view;
    view = camera.GetViewMatrix();
    glm::mat4 projection; 
    projection = glm::perspective(glm::radians(camera.Zoom), (float)WIDTH/(float)HEIGHT, 0.1f, 1000.0f);
    MVPM = projection * view ;

		//setup shader info
    GLuint camdirtyu = glGetUniformLocation(ourShader.Program, "camera_dirty");
    GLuint uniformMatrix = glGetUniformLocation(ourShader.Program, "MVPM");
    GLuint viewMatrix = glGetUniformLocation(ourShader.Program, "VM");
    GLuint LFuniformMatrix = glGetUniformLocation(ourShader.Program, "LFMVPM");
    GLuint uniformMatrixpost = glGetUniformLocation(postShader.Program, "MVPM");
    GLuint atmouniformMatrix = glGetUniformLocation(atmoShader.Program, "MVPM");
    GLuint aspectUniform = glGetUniformLocation(ourShader.Program, "aspect");
		ourShader.Use();   
		glUniform1i(glGetUniformLocation(ourShader.Program, "perlworl"), 0); // set it manually
		glUniform1i(glGetUniformLocation(ourShader.Program, "worl"), 1); // set it manually
		glUniform1i(glGetUniformLocation(ourShader.Program, "curl"), 2); // set it manually
		glUniform1i(glGetUniformLocation(ourShader.Program, "lastFrame"), 3); // set it manually
		glUniform1i(glGetUniformLocation(ourShader.Program, "weather"), 4); // set it manually
		glUniform1i(glGetUniformLocation(ourShader.Program, "atmosphere"), 5); // set it manually

    GLuint uppos = glGetUniformLocation(upscaleShader.Program, "pos");
    GLuint upsize = glGetUniformLocation(upscaleShader.Program, "size");
    GLuint upuniformMatrix = glGetUniformLocation(upscaleShader.Program, "MVPM");
    GLuint upcheck = glGetUniformLocation(upscaleShader.Program, "check");
		glUniform1i(glGetUniformLocation(upscaleShader.Program, "buff"), 0); // set it manually

    GLuint texsize = glGetUniformLocation(postShader.Program, "size");
    GLuint texpos = glGetUniformLocation(postShader.Program, "pos");
    GLuint checku = glGetUniformLocation(ourShader.Program, "check");
    GLuint timeu = glGetUniformLocation(ourShader.Program, "time");
    GLuint atmotimeu = glGetUniformLocation(atmoShader.Program, "time");
    GLuint atmochecku = glGetUniformLocation(atmoShader.Program, "check");
		
    // Game loop
		int check = 0;
    while (!glfwWindowShouldClose(window))
    {
        // Set frame time
        GLfloat currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        timePassed = currentFrame;
        if (timePassed - startTime > 0.25 && frames > 10) {
          //frame rate
          //std::cout<<frames/(timePassed-startTime)<<std::endl;
          //time in milliseconds
          std::cout<<deltaTime*1000.0<<endl;
          startTime = timePassed;
          frames = 0;
        }
        lastFrame = currentFrame;
        frames++;
        // Check and call events
				LFMVPM = MVPM;
        glfwPollEvents();
        Do_Movement();

GLenum err;
				while((err = glGetError()) != GL_NO_ERROR)
{
	std::cout<<err<<std::endl;
}


				glBindFramebuffer(GL_FRAMEBUFFER, skyfbo);
				glViewport(0, 0, 512, 512);

        atmoShader.Use();

        glUniformMatrix4fv(atmouniformMatrix, 1, GL_FALSE, glm::value_ptr(MVPM));
				glUniform1f(atmotimeu, timePassed);
				glUniform1i(atmochecku, check%16);

        glBindVertexArray(VAO);
				
				glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
				
				glActiveTexture(GL_TEXTURE0);
    		glBindTexture(GL_TEXTURE_2D, skytex);
				glGenerateMipmap(GL_TEXTURE_2D);

				glBindFramebuffer(GL_FRAMEBUFFER, fbo);
				glViewport(0, 0, 512, 512);
			
        //glClearColor(0.0f, 0.0f, 0.4f, 1.0f);
        //glClear(GL_COLOR_BUFFER_BIT);
				
				ourShader.Use();
				glUniform1i(checku, check%16);
				glUniform1f(timeu, timePassed);
				glUniform1i(camdirtyu, camera_dirty);
				check++;
        // Pass the matrices to the shader
        glUniformMatrix4fv(uniformMatrix, 1, GL_FALSE, glm::value_ptr(MVPM));
        glUniformMatrix4fv(LFuniformMatrix, 1, GL_FALSE, glm::value_ptr(LFMVPM));
        glUniformMatrix4fv(viewMatrix, 1, GL_FALSE, glm::value_ptr(camera.GetViewMatrix()));
				glUniform1f(aspectUniform, ASPECT);

        glBindVertexArray(VAO);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_3D, perlworltex);
				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_3D, worltex);
				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, curltex);
				glActiveTexture(GL_TEXTURE3);
				glBindTexture(GL_TEXTURE_2D, pongfbotex);//last frame
				glActiveTexture(GL_TEXTURE4);
				glBindTexture(GL_TEXTURE_2D, weathertex);//last frame
				glActiveTexture(GL_TEXTURE5);
				glBindTexture(GL_TEXTURE_2D, skytex);//last frame

        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
				
				glActiveTexture(GL_TEXTURE0);
    		glBindTexture(GL_TEXTURE_2D, fbotex);
				glGenerateMipmap(GL_TEXTURE_2D);

        // Pass the matrices to the shader

				glBindFramebuffer(GL_FRAMEBUFFER, pongfbo);
				glViewport(0, 0, 512, 512);
				// Draw the triangle
        //postShader.Use();
				upscaleShader.Use();
				//glUniform1f(texsize, 1.0);
        // Pass the matrices to the shader
        glUniformMatrix4fv(upuniformMatrix, 1, GL_FALSE, glm::value_ptr(MVPM));
				glUniform1i(upcheck, (check+15)%16);
        //glClearColor(0.2f, 0.4f, 0.8f, 1.0f);

        //glClear(GL_COLOR_BUFFER_BIT);
        glBindVertexArray(VAO);
				
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, fbotex);

				glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);


				glBindFramebuffer(GL_FRAMEBUFFER, 0);
				
				glActiveTexture(GL_TEXTURE0);
    		//glBindTexture(GL_TEXTURE_2D, fbotex);
				//glGenerateMipmap(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, pongfbotex);
				glGenerateMipmap(GL_TEXTURE_2D);
				
				glViewport(0, 0, WIDTH, HEIGHT);


		   // Render
        // Clear the colorbuffer
        glClearColor(0.2f, 0.4f, 0.8f, 1.0f);
        //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        // Draw the triangle
        postShader.Use();

				glUniform1f(texsize, 1.0);
        // Pass the matrices to the shader
        glUniformMatrix4fv(uniformMatrixpost, 1, GL_FALSE, glm::value_ptr(MVPM));
				

        glBindVertexArray(VAO);
				
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, pongfbotex);

				glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
				
        // Swap the screen buffers
        glfwSwapBuffers(window);
				camera_dirty++;
    }
    // Properly de-allocate all resources once they've outlived their purpose
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    // Terminate GLFW, clearing any resources allocated by GLFW.
    glfwTerminate();
    return 0;
}

void Do_Movement()
{
    // Camera controls
    if(keys[GLFW_KEY_W])
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if(keys[GLFW_KEY_S])
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if(keys[GLFW_KEY_A])
        camera.ProcessKeyboard(LEFT, deltaTime);
    if(keys[GLFW_KEY_D])
        camera.ProcessKeyboard(RIGHT, deltaTime);

    glm::mat4 view;
    view = camera.GetViewMatrix();
    glm::mat4 projection; 
    projection = glm::perspective(glm::radians(camera.Zoom), (float)WIDTH/(float)HEIGHT, 0.1f, 1000.0f);
    MVPM = projection * view;
}

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    //cout << key << endl;
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
    if (key >= 0 && key < 1024)
    {
        if(action == GLFW_PRESS)
            keys[key] = true;
        else if(action == GLFW_RELEASE)
            keys[key] = false;  
    }
		camera_dirty = 0;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if(firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    GLfloat xoffset = xpos - lastX;
    GLfloat yoffset = lastY - ypos;  // Reversed since y-coordinates go from bottom to left
    
    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);

    glm::mat4 view;
    view = camera.GetViewMatrix();
    glm::mat4 projection; 
    projection = glm::perspective(glm::radians(camera.Zoom), (float)WIDTH/(float)HEIGHT, 0.1f, 1000.0f);
    MVPM = projection * view;
		camera_dirty = 0;
} 


void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);

    glm::mat4 view;
    view = camera.GetViewMatrix();
    glm::mat4 projection; 
    projection = glm::perspective(glm::radians(camera.Zoom), (float)WIDTH/(float)HEIGHT, 0.1f, 1000.0f);
    MVPM = projection * view;
		camera_dirty = true;
}

// the remap function used in the shaders as described in Gpu Pro 8. It must match when using pre packed textures
float remap(float originalValue, float originalMin, float originalMax, float newMin, float newMax)
{
	return newMin + (((originalValue - originalMin) / (originalMax - originalMin)) * (newMax - newMin));
}
