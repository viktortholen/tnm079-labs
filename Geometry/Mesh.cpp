/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 * Acknowledgements for original code base:
 * - Gunnar Johansson
 * - Ken Museth
 * - Michael Bang Nielsen
 * - Ola Nilsson
 * - Andreas Soderstrom
 *
 * Code updated in the period 2017-2018 by Jochen Jankowai
 *
 *************************************************************************************************/
#include <Geometry/Mesh.h>
#include <iostream>

const Mesh::VisualizationMode Mesh::CurvatureVertex =
    NewVisualizationMode("Vertex curvature");
const Mesh::VisualizationMode Mesh::CurvatureFace =
    NewVisualizationMode("Face curvature");

float Mesh::Area() const {
  std::cerr << "Error: area() not implemented for this Mesh" << std::endl;
  return -1;
}

float Mesh::Volume() const {
  std::cerr << "Error: volume() not implemented for this Mesh" << std::endl;
  return -1;
}

size_t Mesh::Genus() const {
  std::cerr << "Error: genus() not implemented for this Mesh" << std::endl;
  return -1;
}
