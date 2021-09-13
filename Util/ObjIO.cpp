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
#include "ObjIO.h"
#include "Geometry/HalfEdgeMesh.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Util/Util.h>

bool ObjIO::Load(Mesh *mesh, std::istream &is) {
  // std::cerr << "Reading obj file.\nOutputting any skipped line(s) for
  // reference.\n";
  bool success = ReadHeader(is);
  if (!success) {
    return false;
  }

  success = ReadData(is);
  if (!success) {
    return false;
  }

  // Build mesh
  const size_t numTris = loadData.tris.size();
  for (size_t t = 0; t < numTris; t++) {
    glm::uvec3 &triangle = loadData.tris[t];
    std::vector<glm::vec3 > verts;
    verts.push_back(loadData.verts[triangle[0]]);
    verts.push_back(loadData.verts[triangle[1]]);
    verts.push_back(loadData.verts[triangle[2]]);

    mesh->AddFace(verts);
  }
  return true;
}

bool ObjIO::ReadHeader(std::istream &is) {
  std::string buf;
  // just read to the first line starting with a 'v'
  while (!is.eof() && is.peek() != 'v') {
    getline(is, buf);
    // std::cerr << "\"" << buf << "\"\n";
  }
  if (is.good())
    return true;
  else
    return false;
}

bool ObjIO::ReadData(std::istream &is) {
  std::string lineBuf;
  int c;
  int i = 0;
  while (!is.eof()) {
    c = is.peek();
    switch (c) {
    case 'V':
    case 'v': {
      std::string startBuf;
      is >> startBuf;       // get the start of the line
      getline(is, lineBuf); // get the rest of the line
      if (startBuf == "v") {
        loadData.verts.push_back(vec3<float>(lineBuf));
      }
    } break;
    case 'F':
    case 'f': {
      std::stringstream buf;
      is.get(*buf.rdbuf(), '\n'); // read a line into buf
      is.get();                   // read the not extracted \n
      buf << "\n";                // and add it to the string stream

      std::string tmp;
      buf >> tmp; // get the first f or F (+ whitespace)

      // count the number of faces, delimited by whitespace
      int count = 0;
      while (buf >> tmp) {
        count++;
      }
      // reset stream
      buf.clear();
      buf.seekg(0, std::ios::beg);

      // Determine wheter we have a triangle or a quad
      if (count == 3) {
        loadData.tris.push_back(ReadTri(buf));
      } else {
        std::cerr << "Encountered polygon with " << count
                  << " faces. Bailing out.\n";
        return false;
      }
    } break;
    default:
      // otherwise just skip the row
      getline(is, lineBuf);
      // output it so we see what we miss :)
      // std::cerr << "\"" << lineBuf << "\"\n";
      break;
    }
    i++;
  }
  return true;
}

glm::uvec3 ObjIO::ReadTri(std::istream &is) {
  //  This is a simplified version of an obj reader that can't read normal and
  //  texture indices
  std::string buf, v;
  is >> buf;
  assert(buf == "f" || buf == "F");

  getline(is, v); // read indices
  return vec3<unsigned int>(v) -
         glm::uvec3(1, 1, 1); // obj file format is 1-based
}
