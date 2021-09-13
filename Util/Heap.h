/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#ifndef _HEAP
#define _HEAP

#include <iostream>
#include <limits>
#include <vector>

class Heap {
public:
  Heap();
  ~Heap();

  struct Heapable {
    Heapable() : cost(0), position((std::numeric_limits<size_t>::max)()) {}
    Heapable(float val)
        : cost(val), position((std::numeric_limits<size_t>::max)()) {}

    float cost;
    size_t position;

    bool operator<(const Heapable &h) const { return this->cost < h.cost; }
  };

  void push(Heapable *h);

  Heapable *peek();

  Heapable *pop();

  Heapable *remove(Heapable *h);

  inline size_t size() { return mNodes.size() - 1; }

  inline bool isEmpty() { return size() == 0; }

  void update(Heapable *h);

  void print(std::ostream &os);

protected:
  inline size_t parent(size_t i) { return i / 2; }
  inline size_t leftChild(size_t i) { return 2 * i; }
  inline size_t rightChild(size_t i) { return 2 * i + 1; }

  void percolateUp(size_t hole);
  void percolateDown(size_t hole);

  std::vector<Heapable *> mNodes;
};

#endif
