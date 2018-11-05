
#ifndef __renderimage_h__
#define __renderimage_h__

#include "render_standalone.h"

const unsigned WIDTH = 1280;
const unsigned HEIGHT = 720;


void setCameraPosition(FieldData domainMin[3], FieldData domainMax[3]);


void renderImage(int numLines,
                 FieldData* rho,
                 FieldData* pressure,
                 FieldData* velocity,
                 FieldData* centerCoordinates,
                 FieldData* temperature,
                 FieldData domainMin[3],
                 FieldData domainMax[3]);


#endif // __renderimage_h__

