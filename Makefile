#OBJS specifies which files to compile as part of the project
OBJS = main.cpp #include/FastNoise/FastNoise.cpp include/TileableVolumeNoise/TileableVolumeNoise.cpp

#CC specifies which compiler we're using
CC = g++

#INCLUDE_PATHS specifies the additional include paths we'll need
INCLUDE_PATHS = -I -Iinclude #-Iinclude/TileableVolumeNoise/glm#-IC:\glew\include -IC:\SOIL\include

#LIBRARY_PATHS specifies the additional library paths we'll need
LIBRARY_PATHS = -L/usr/local/lib#-LC:\GLFW\lib #-LC:\glew\lib -LC:\SOIL\lib

#COMPILER_FLAGS specifies the additional compilation options we're using
COMPILER_FLAGS = -std=c++11 

#LINKER_FLAGS specifies the libraries we're linking against
LINKER_FLAGS =  -lXi -lGLEW -lGLU -lm -lGL -lm -lpthread -ldl -ldrm -lXdamage -lX11-xcb -lxcb-glx -lxcb-dri2 -lglfw3 -lrt -lm -ldl -lXrandr -lXinerama -lXxf86vm -lXext -lXcursor -lXrender -lXfixes -lX11 -lpthread -lxcb -lXau -lXdmcp
#-lGLEW -lglfw3 -lGL -lm -lXrandr -lXi -lX11 -lXxf86vm -lpthread #-lGL -lglfw3#-lglew32 -lmingw32 -lSOIL -lopengl32 -lglfw3  -lglu32 -lgdi32
#OBJ_NAME specifies the name of our exectuable
OBJ_NAME = main

#This is the target that compiles our executable
all : $(OBJS)
	$(CC) $(OBJS) $(INCLUDE_PATHS) $(LIBRARY_PATHS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)
