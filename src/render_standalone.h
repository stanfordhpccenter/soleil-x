
#ifndef __render_standalone_h__
#define __render_standalone_h__

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "OpenGL/glu.h"
#else
#include "GL/osmesa.h"
#include "GL/glu.h"
#endif


const int WIDTH = 1280;
const int HEIGHT = 720;

#endif // __render_standalone_h__

