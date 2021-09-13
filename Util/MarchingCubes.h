/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas S�derstr�m (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#ifndef __marching_cubes_h__
#define __marching_cubes_h__

#include "Util/Util.h"
#include <vector>
#include <glm.hpp>

//! Method to triangulate a voxel
std::vector<glm::vec3 > Triangulate(float voxelValues[8], float i, float j,
                                        float k, float delta);

#endif
