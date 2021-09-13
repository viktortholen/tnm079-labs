# TNM079 - Modeling and Animation
Hi and welcome to the course TNM079 - Modeling and Animation!

This repository contains the code for the application that we will use to implement the topics discussed in the lecture. Below you will find instructions on how to set everything up. If you still need help, feel free to contact the lab assistant.

If you are familiar with ```git``` we would recommend that you simply fork the repo so you have your own copy that you can work on in your team.

## Setting up the development environment
The user interface was implemented using wxWidgets. In addition to that, we will need GLUT and GLEW. These are external libraries that we need to set up in order to compile our code. Don't worry, it sounds more difficult than it is, especially if you are using Windows.

#### Windows
For Windows, assuming you use Visual Studio, the libraries are provided so you don't need to do anything. We have not been able to successfully test an environment with [CLion](https://www.jetbrains.com/clion/) and MinGW or Cygwin.

#### Mac
For Mac, we tested compiling the code with [CLion](https://www.jetbrains.com/clion/) (free for students) after installing the wxWidgets and glew libraries via brew. If you don't have [brew](https://brew.sh/) installed, you can do so typing
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
in the terminal. After the installation finished you can use it like this to install the libraries:
```
brew install wxmac
```
You will also need to install glew, e.g.
```
brew install glew
```
#### Linux
On Linux, we tested compiling the code with [CLion](https://www.jetbrains.com/clion/) (free for students). If you don't already have a working development environment for OpenGL, you'll need to install some libraries to begin with, i.e. GLUT and GLEW like so:
```
sudo apt-get install freeglut3 freeglut3-dev
```
```
sudo apt-get install libglew-dev
```
After this, we can install the wxWidgets libraries like this:
```
sudo apt-get install libwxgtk3.0-dev
```

## Getting the code
The source code for the application can be downloaded as zip [here](https://gitlab.itn.liu.se/TNM079/tnm079-student/repository/master/archive.zip) or you can clone the [repo](https://gitlab.itn.liu.se/TNM079/tnm079-student.git).

Note that
* if you want to fork/clone the repo, you will have to create an account on ITNs Gitlab server using your student credentials.
* it would be best if you work on your student network drive and not on the local one.

## Compiling the code
In the following, we will shortly describe how to compile the code on different platforms.

#### Windows
We will use CMake to create the project files for Visual Studio. If you do not have CMake installed, this is the first thing you will need to do. Start it and set the source code directory (first field on the top) to where you have unzipped the source code or cloned the repo, respectively, e.g. ```I:/tnm079/tnm079-student```. Subsequently, set the build directory (second field on the top) to something like ```I:/tnm079/build/```. Hit ```Configure```. A window will pop up, asking you what kind of project file you wish to generate. Please check, which version of Visual Studio is installed. At the time of writing this document, on the lab computers it is *__Visual Studio 14 2015 Win64__*. Press ```Finish```. If an error persists and you don't know how to fix it feel free to ask the lab assistant. After the configuration step terminated successfully you can hit ```Generate``` and open the generated solution.

##### Compilation targets in Visual Studio
You can compile the code in ```Debug```, ```Release```, or ```RelWithDebInfo``` (Release with debug information) mode. You can change this in the top toolbar. There is a dropdown menu that is per default set to ```Debug```. While you are working on the solution it might be advisable to use ```Debug``` or ```RelWithDebInfo``` so you can (you guessed it) debug the code, use breakpoints, etc. When you are done, you can switch to ```Release```. This has the effect that the application will run significantly faster.

#### Mac and Linux
If you have completed the setup and have CLion installed you can simply click on "Open project" on startup and everything should work just fine.
##### Compilation targets in CLion
You can compile the code in ```Debug```, ```Release```, or ```RelWithDebInfo``` (Release with debug information) mode. By default, only ```Debug``` is available (as you can see in the upper right corner). You can change this by going to ```File``` &rarr; ```Settings``` &rarr; ```Build, Execution, Deployment``` &rarr; ```CMake``` and hitting the plus button to add compilation targets. Per default there is only ```Debug```. Pressing the plus button twice will add the targets ```Release``` and ```RelWithDebInfo```. Now you will be able to switch between these targets in the upper right corner. While you are working on the solution it might be advisable to use ```Debug``` or ```RelWithDebInfo``` so you can (you guessed it) debug the code, use breakpoints, etc. When you are done, you can switch to ```Release```. This has the effect that the application will run significantly faster.

## General notes
Note that in the CMake options you can set which lab parts you want to compile. So if you are working on lab 2 you don't need to compile parts 3-6. When you move from one lab to the other later on, make sure you enable the next lab. Generally, you might as well enable all the parts right away, the code will just need a little longer to compile.

## Technical notes
wxWidgets Windows libraries compiled from git branch ```WX_2_8_BRANCH``` using the ```X64 Native Tools Command Prompt for VS 2017``` using the following arguments:

```nmake -f ./makefile.vc BUILD=debug TARGET_CPU=AMD64 SHARED=0 USE_OPENGL=1 "USER_C_FLAGS=/D_ITERATOR_DEBUG_LEVEL=0"```

## Lab notes
### Lab 1
- After you have finished lab 1, the object that you loaded may be black. Don't panic yet, you just need to set a colormap. You can do this by selecting the added object in the list
and then choosing one of the color maps.
- In the output log, it should say that 0 isolated vertices were found.
- If you want to compute the genus, remember that in the halfedge data structure, there are two halfedges stored for every real edge.
- When calculating the the genus, remember that `size_t` (e.g., returned by `GetNumVerts()`) is an unsigned type. This is problematic when intermediate results are negative (such as `V-E`). Make sure to use a signed type instead.

### Lab 2
- If you need the camera position for the other cost heurisitc, you can obtain it by including `GUI/GLViewer.h` and calling `GLViewer::GetCamera().GetPosition();`.

### Lab 3

### Lab 4

### Lab 5
- If your objects disappear when smoothing, check that the denominator of equation 16 is not equal to 0.

### Lab 6
- If you cannot see the fluid whatsoever after loading template 1 (i.e, if it doesn't look like the image below), check your ```Quadric.cpp``` implementation. For both, ```GetValue()``` and ```GetGradient()```, the computations need to be performed in object space.
- If the level-set doesn't move at all, check in the ```Evaluate()``` method in ```OperatorAdvect.h``` that you sample the vector field with ```x```, ```y```, ```z``` BUT calculate the gradient using ```i```, ```j```, ```k```.
- If your fluid falls through the bottom even though your boundary condition is correct, make sure you pay attention to the operator precedence when calculating the b vector, e.g., ```a/b*c != a/(b*c)```.
- If your water magically loses volume, make sure that in the ```ExternalForces()``` method you call ```GetValue()``` on ```mVelocityField``` and not on the ```FluidSolver``` itself.
- If your fluid starts to "crinkle", make sure you multiply the values for `A` with -1.
# TNM079-labs
