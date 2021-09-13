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
#ifndef _DECIMATION_INTERFACE
#define _DECIMATION_INTERFACE

#include <cstddef>

class DecimationInterface
{
public :
  DecimationInterface() { }
  virtual ~DecimationInterface() { }

  virtual bool decimate() = 0;

  virtual bool decimate(size_t targetFaces) = 0;


};

#endif
