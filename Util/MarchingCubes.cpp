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
#include "MarchingCubes.h"

/*!
 * Grabbed from:
 * http://astronomy.swin.edu.au/~pbourke/modelling/polygonise/
 * Modified slightly.
 *
 * NB! Uses clockwise orientation.
 */
std::vector<glm::vec3> Triangulate(float voxelValues[8], float i, float j,
                                        float k, float delta) {
#include "MarchingCubesTable.h"
  int cubeindex = 0;
  static glm::vec3 vertlist[12];
  std::vector<glm::vec3 > verts;

  if (voxelValues[0] < 0.f)
    cubeindex |= 1;
  if (voxelValues[1] < 0.f)
    cubeindex |= 2;
  if (voxelValues[2] < 0.f)
    cubeindex |= 4;
  if (voxelValues[3] < 0.f)
    cubeindex |= 8;
  if (voxelValues[4] < 0.f)
    cubeindex |= 16;
  if (voxelValues[5] < 0.f)
    cubeindex |= 32;
  if (voxelValues[6] < 0.f)
    cubeindex |= 64;
  if (voxelValues[7] < 0.f)
    cubeindex |= 128;

  /* Cube is entirely in/out of the surface */
  if (edgeTable[cubeindex] == 0)
    return verts;

  /* Find the vertices where the surface intersects the cube */
  if (edgeTable[cubeindex] & 1)
    // 0-1
    vertlist[0] =
        glm::vec3(i + Root(voxelValues[0], voxelValues[1]) * delta, j, k);
  if (edgeTable[cubeindex] & 2)
    // 1-2
    vertlist[1] = glm::vec3(
        i + delta, j + Root(voxelValues[1], voxelValues[2]) * delta, k);
  if (edgeTable[cubeindex] & 4)
    // 3-2
    vertlist[2] = glm::vec3(
        i + Root(voxelValues[3], voxelValues[2]) * delta, j + delta, k);
  if (edgeTable[cubeindex] & 8)
    // 0-3
    vertlist[3] =
        glm::vec3(i, j + Root(voxelValues[0], voxelValues[3]) * delta, k);
  if (edgeTable[cubeindex] & 16)
    //  4-5
    vertlist[4] = glm::vec3(
        i + Root(voxelValues[4], voxelValues[5]) * delta, j, k + delta);
  if (edgeTable[cubeindex] & 32)
    //  5-6
    vertlist[5] = glm::vec3(
        i + delta, j + Root(voxelValues[5], voxelValues[6]) * delta, k + delta);
  if (edgeTable[cubeindex] & 64)
    //  7-6
    vertlist[6] = glm::vec3(
        i + Root(voxelValues[7], voxelValues[6]) * delta, j + delta, k + delta);
  if (edgeTable[cubeindex] & 128)
    //  4-7
    vertlist[7] = glm::vec3(
        i, j + Root(voxelValues[4], voxelValues[7]) * delta, k + delta);
  if (edgeTable[cubeindex] & 256)
    // 0-4
    vertlist[8] =
        glm::vec3(i, j, k + Root(voxelValues[0], voxelValues[4]) * delta);
  if (edgeTable[cubeindex] & 512)
    // 1-5
    vertlist[9] = glm::vec3(
        i + delta, j, k + Root(voxelValues[1], voxelValues[5]) * delta);
  if (edgeTable[cubeindex] & 1024)
    // 2-6
    vertlist[10] = glm::vec3(
        i + delta, j + delta, k + Root(voxelValues[2], voxelValues[6]) * delta);
  if (edgeTable[cubeindex] & 2048)
    // 3-7
    vertlist[11] = glm::vec3(
        i, j + delta, k + Root(voxelValues[3], voxelValues[7]) * delta);

  /* Create the triangle */
  for (size_t m = 0; triTable[cubeindex][m] != -1; m += 3) {
    // Use counter clockwise orientation
    verts.push_back(vertlist[triTable[cubeindex][m]]);
    verts.push_back(vertlist[triTable[cubeindex][m + 2]]);
    verts.push_back(vertlist[triTable[cubeindex][m + 1]]);
  }

  return verts;
}
