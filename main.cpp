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


// Window dimensions
const GLuint WIDTH = 512, HEIGHT = 512;
const GLuint downscale = 4; //4 is best//any more and the gains dont make up for the lag
GLuint downscalesq = downscale*downscale;
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

    // Build and compile our shader programs
    Shader ourShader("sky.vert", "sky.frag");
		Shader postShader("tex.vert", "tex.frag");
		Shader upscaleShader("upscale.vert", "upscale.frag");
		Shader atmoShader("atmo.vert", "atmo.frag");

    // Set up vertex data (and buffer(s)) and attribute pointers
    GLfloat vertices[] = {
        // Positions      
       -1.0f, -1.0f,
       -1.0f, 3.0f,
       3.0f,  -1.0f,
    };

    GLuint VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    // Position attribute
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);


		//our main full size framebuffer
		GLuint fbo, fbotex;

    glGenFramebuffers(1, &fbo);
    glGenTextures(1, &fbotex);
    glBindTexture(GL_TEXTURE_2D, fbotex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, WIDTH, HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbotex, 0);


		//our secondary full size framebuffer for copying and reading from the main framebuffer
		GLuint copyfbo, copyfbotex;
		
    glGenFramebuffers(1, &copyfbo);
    glGenTextures(1, &copyfbotex);
    glBindTexture(GL_TEXTURE_2D, copyfbotex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, WIDTH, HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);	
    glBindTexture(GL_TEXTURE_2D, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, copyfbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, copyfbotex, 0);


		//our downscaled buffer that we actually render to
		GLuint subbuffer, subbuffertex;
		
    glGenFramebuffers(1, &subbuffer);
    glGenTextures(1, &subbuffertex);
    glBindTexture(GL_TEXTURE_2D, subbuffertex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, WIDTH/downscale, HEIGHT/downscale, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, subbuffer);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, subbuffertex, 0);


		//small buffer for rendering sky background//should be removed and incorporated into sky
		GLuint skyfbo, skytex;
		
    glGenFramebuffers(1, &skyfbo);
    glGenTextures(1, &skytex);
    glBindTexture(GL_TEXTURE_2D, skytex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, WIDTH/downscale, HEIGHT/downscale, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
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
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glGenerateMipmap(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, 0);
		stbi_image_free(worlNoiseArray);

		unsigned char *perlWorlNoiseArray = stbi_load("assets/perlworlnoise.tga", &x, &y, &n, 4);

		glGenTextures(1, &perlworltex);
    glBindTexture(GL_TEXTURE_3D, perlworltex);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, 128, 128, 128, 0, GL_RGBA, GL_UNSIGNED_BYTE, perlWorlNoiseArray);
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
    ourShader.Use();   
    GLuint uniformMatrix = glGetUniformLocation(ourShader.Program, "MVPM");
    GLuint aspectUniform = glGetUniformLocation(ourShader.Program, "aspect");
    GLuint checku = glGetUniformLocation(ourShader.Program, "check");
    GLuint timeu = glGetUniformLocation(ourShader.Program, "time");
    GLuint resolutionu = glGetUniformLocation(ourShader.Program, "resolution");
    GLuint downscaleu = glGetUniformLocation(ourShader.Program, "downscale");
    GLuint perlworluniform = glGetUniformLocation(ourShader.Program, "perlworl");
    GLuint worluniform = glGetUniformLocation(ourShader.Program, "worl");
    GLuint curluniform = glGetUniformLocation(ourShader.Program, "curl");
    GLuint weatheruniform = glGetUniformLocation(ourShader.Program, "weather");
    GLuint atmosphereuniform = glGetUniformLocation(ourShader.Program, "atmosphere"); 
		
		upscaleShader.Use();
    GLuint upuniformMatrix = glGetUniformLocation(upscaleShader.Program, "MVPM");
    GLuint upLFuniformMatrix = glGetUniformLocation(upscaleShader.Program, "LFMVPM");
    GLuint upcheck = glGetUniformLocation(upscaleShader.Program, "check");
    GLuint upresolution = glGetUniformLocation(upscaleShader.Program, "resolution");
    GLuint updownscale = glGetUniformLocation(upscaleShader.Program, "downscale");
    GLuint buffuniform = glGetUniformLocation(upscaleShader.Program, "buff");
    GLuint ponguniform = glGetUniformLocation(upscaleShader.Program, "pong");

		atmoShader.Use();
    GLuint atmotimeu = glGetUniformLocation(atmoShader.Program, "time");
    GLuint atmochecku = glGetUniformLocation(atmoShader.Program, "check");
    GLuint atmoresolution = glGetUniformLocation(atmoShader.Program, "resolution");
    GLuint atmouniformMatrix = glGetUniformLocation(atmoShader.Program, "MVPM");
		
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
        LFMVPM = MVPM;
        // Check and call events
        glfwPollEvents();
        Do_Movement();

				GLenum err;
        while((err = glGetError()) != GL_NO_ERROR)
				{
    			std::cout<<err<<std::endl;
				}

/*
        glBindFramebuffer(GL_FRAMEBUFFER, skyfbo);
        glViewport(0, 0, WIDTH/downscale, HEIGHT/downscale);

        atmoShader.Use();

        glUniformMatrix4fv(atmouniformMatrix, 1, GL_FALSE, glm::value_ptr(MVPM));
        glUniform1f(atmotimeu, timePassed);
        glUniform2f(atmoresolution, float(WIDTH/downscale), float(HEIGHT/downscale));

        glBindVertexArray(VAO);        
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, skytex);
        glGenerateMipmap(GL_TEXTURE_2D);
        */        
				
				//Write to quarter scale buffer
        glBindFramebuffer(GL_FRAMEBUFFER, subbuffer);
        glViewport(0, 0, WIDTH/downscale, HEIGHT/downscale);
                    
        ourShader.Use();

        glUniform1f(timeu, timePassed);
        glUniformMatrix4fv(uniformMatrix, 1, GL_FALSE, glm::value_ptr(MVPM));
        glUniform1f(aspectUniform, ASPECT);
        glUniform1i(checku, (check)%(downscalesq));
        glUniform2f(resolutionu, GLfloat(WIDTH), GLfloat(HEIGHT));
        glUniform1f(downscaleu, GLfloat(downscale));

        glUniform1i(perlworluniform, 0);
        glUniform1i(worluniform, 1);
        glUniform1i(curluniform, 2);
        glUniform1i(weatheruniform, 3);
        glUniform1i(atmosphereuniform, 4);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_3D, perlworltex);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_3D, worltex);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, curltex);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, weathertex);
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, skytex);

        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);
                

				//upscale the buffer into full size framebuffer
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glViewport(0, 0, WIDTH, HEIGHT);

        upscaleShader.Use();
        glUniformMatrix4fv(upLFuniformMatrix, 1, GL_FALSE, glm::value_ptr(LFMVPM));
        glUniformMatrix4fv(upuniformMatrix, 1, GL_FALSE, glm::value_ptr(MVPM));
        glUniform1i(upcheck, (check)%downscalesq);
        glUniform2f(upresolution, GLfloat(WIDTH), GLfloat(HEIGHT));
        glUniform1f(updownscale, GLfloat(downscale));

        glUniform1i(buffuniform, 0);
        glUniform1i(ponguniform, 1);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, subbuffertex);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, copyfbotex);
                
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);


				//copy the full size buffer so it can be read from next frame
        glBindFramebuffer(GL_FRAMEBUFFER, copyfbo);

        postShader.Use();
                
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, fbotex);

        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);


				//copy to the screen
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, WIDTH, HEIGHT);

        postShader.Use();
                
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, fbotex);

        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);
                
        // Swap the screen buffers
        glfwSwapBuffers(window);
        check++;
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
} 


void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);

    glm::mat4 view;
    view = camera.GetViewMatrix();
    glm::mat4 projection; 
    projection = glm::perspective(glm::radians(camera.Zoom), (float)WIDTH/(float)HEIGHT, 0.1f, 1000.0f);
    MVPM = projection * view;
}

// the remap function used in the shaders as described in Gpu Pro 8. It must match when using pre packed textures
float remap(float originalValue, float originalMin, float originalMax, float newMin, float newMax)
{
	return newMin + (((originalValue - originalMin) / (originalMax - originalMin)) * (newMax - newMin));
}
