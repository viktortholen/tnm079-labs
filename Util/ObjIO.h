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
#ifndef __obj_io_h__
#define __obj_io_h__

#include "Geometry/Mesh.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <glm.hpp>

class ObjIO {
public:
  ObjIO() {}

  bool Load(Mesh *, std::istream &is); // false return on error

protected:
  bool ReadHeader(std::istream &is);
  bool ReadData(std::istream &is);

  static glm::uvec3 ReadTri(std::istream &is);

  struct LoadData {
    std::vector<glm::vec3> verts;
    std::vector<glm::uvec3> tris;
  } loadData;
};

#endif
