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
#ifndef _SIMPLE_DECIMATION_MESH
#define _SIMPLE_DECIMATION_MESH

#include "Decimation/DecimationMesh.h"

class SimpleDecimationMesh : public virtual DecimationMesh {
public:
  SimpleDecimationMesh() {}
  virtual ~SimpleDecimationMesh() {}

protected:
  virtual void computeCollapse(EdgeCollapse *collapse);
};

#endif
