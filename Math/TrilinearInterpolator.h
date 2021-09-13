#pragma once

#include <Math/Volume.h>

class TrilinearInterpolator  {
public:
  TrilinearInterpolator();
  ~TrilinearInterpolator();

  template <typename T>
      T Interpolate(float x, float y, float z,
          const Volume<T>& grid) {
      auto i = static_cast<size_t>(x);
      auto j = static_cast<size_t>(y);
      auto k = static_cast<size_t>(z);

      float bx = x - static_cast<float>(i);
      float by = y - static_cast<float>(j);
      float bz = z - static_cast<float>(k);

      return (grid.GetValue(i, j, k) * (1 - bx) * (1 - by) * (1 - bz) +
          grid.GetValue(i + 1, j, k) * bx * (1 - by) * (1 - bz) +
          grid.GetValue(i + 1, j + 1, k) * bx * by * (1 - bz) +
          grid.GetValue(i, j + 1, k) * (1 - bx) * by * (1 - bz) +
          grid.GetValue(i, j, k + 1) * (1 - bx) * (1 - by) * bz +
          grid.GetValue(i + 1, j, k + 1) * bx * (1 - by) * bz +
          grid.GetValue(i + 1, j + 1, k + 1) * bx * by * bz +
          grid.GetValue(i, j + 1, k + 1) * (1 - bx) * by * bz);
  }
};

