// #include <dirt/box.h>
// #include <dirt/vec.h>

// #include <algorithm>
// #include <vector>
// #include <dirt/kdtree.h>

// std::vector<Photon> _nodes;

// void recursiveTreeBuild(uint32_t dst, uint32_t start, uint32_t end)
// {
//     if (end == start) {
//         // Leaf node
//         _nodes[dst].setSplitInfo(0, 0, 0);
//         return;
//     } else if (end - start == 1) {
//         // Single child only. Special case
//         if (_nodes[dst].pos.x < _nodes[start].pos.x)
//             std::swap(_nodes[dst], _nodes[start]);
//         _nodes[dst].setSplitInfo(start, 0, 1);
//         _nodes[start].setSplitInfo(0, 0, 0);
//         return;
//     }

//     Box3f bounds(_nodes[dst].pos);
//     for (uint32_t i = start; i < end; ++i)
//         bounds.enclose(_nodes[i].pos);
//     uint32_t splitDim = maxDim(bounds.diagonal());

//     std::sort(_nodes.begin() + start, _nodes.begin() + end, [&](const Photon &a, const Photon &b) {
//         return a.pos[splitDim] < b.pos[splitDim];
//     });

//     uint32_t splitIdx = start + (end - start + 1)/2;
//     float rightPlane = _nodes[splitIdx].pos[splitDim];
//     float  headPlane = _nodes[dst].pos[splitDim];
//     float  leftPlane = _nodes[splitIdx - 1].pos[splitDim];

//     if (headPlane < leftPlane || headPlane > rightPlane) {
//         uint32_t swapIdx = headPlane > rightPlane ? splitIdx : splitIdx - 1;
//         std::swap(_nodes[dst], _nodes[swapIdx]);
//     }

//     uint32_t childIdx = start;
//     if (splitIdx > childIdx + 1)
//         std::swap(_nodes[childIdx + 1], _nodes[splitIdx]);

//     recursiveTreeBuild(childIdx + 0, start + 2, splitIdx + 1);
//     recursiveTreeBuild(childIdx + 1, splitIdx + 1, end);

//     _nodes[dst].setSplitInfo(childIdx, splitDim, 2);
// }

// void KdTree::addPhoton(Photon p)
// {
//     _nodes.emplace_back(p);
// }

// void KdTree::buildTree()
// {
//     if (!_nodes.empty())
//         recursiveTreeBuild(0, 1, _nodes.size());
// }

// void KdTree::nearestNeighbours(Vec3f pos, std::vector<const Photon *> &photons, std::vector<float> &distances, int k, const float maxDist) const
// {
//     photons.clear();
//     distances.clear();
//     photons.reserve(k);
//     distances.reserve(k);

//     if (_nodes.empty())
//         return;

//     float maxDistSq = maxDist*maxDist;

//     const Photon *stack[28];
//     const Photon **stackPtr = stack;

//     const Photon *current = &_nodes[0];
//     while (true) {
//         float dSq = length2(current->pos - pos);
//         if (dSq < maxDistSq) {
//             if (int(photons.size()) < k) {
//                 photons.push_back(current);
//                 distances.push_back(dSq);

//                 if (int(photons.size()) == k) {
//                     // Build max heap
//                     const int halfK = k/2;
//                     for (int i = halfK - 1; i >= 0; --i) {
//                         int parent = i;
//                         const Photon *reloc = photons[i];
//                         float relocDist = distances[i];
//                         while (parent < halfK) {
//                             int child = parent*2 + 1;
//                             if (child < k - 1 && distances[child] < distances[child + 1])
//                                 child++;
//                             if (relocDist >= distances[child])
//                                 break;
//                             photons[parent] = photons[child];
//                             distances[parent] = distances[child];
//                             parent = child;
//                         }
//                         photons[parent] = reloc;
//                         distances[parent] = relocDist;
//                     }
//                     maxDistSq = distances[0];
//                 }
//             } else {
//                 const int halfK = k/2;
//                 int parent = 0;
//                 while (parent < halfK) {
//                     int child = parent*2 + 1;
//                     if (child < k - 1 && distances[child] < distances[child + 1])
//                         child++;
//                     if (dSq >= distances[child])
//                         break;
//                     photons[parent] = photons[child];
//                     distances[parent] = distances[child];
//                     parent = child;
//                 }
//                 photons[parent] = current;
//                 distances[parent] = dSq;
//                 maxDistSq = distances[0];
//             }
//         }

//         uint32_t splitDim = current->splitDim();
//         float planeDist = pos[splitDim] - current->pos[splitDim];
//         bool traverseLeft  = current->hasLeftChild () && (planeDist <= 0.0f || planeDist*planeDist < maxDistSq);
//         bool traverseRight = current->hasRightChild() && (planeDist >= 0.0f || planeDist*planeDist < maxDistSq);

//         uint32_t childIdx = current->childIdx();
//         if (traverseLeft && traverseRight) {
//             if (planeDist <= 0.0f) {
//                 *stackPtr++ = &_nodes[childIdx + 1];
//                 current = &_nodes[childIdx];
//             } else {
//                 *stackPtr++ = &_nodes[childIdx];
//                 current = &_nodes[childIdx + 1];
//             }
//         } else if (traverseLeft) {
//             current = &_nodes[childIdx];
//         } else if (traverseRight) {
//             current = &_nodes[childIdx + 1];
//         } else {
//             if (stackPtr == stack)
//                 break;
//             else
//                 current = *--stackPtr;
//         }
//     }

//     for (auto &f : distances)
//         f = std::sqrt(f);
// }