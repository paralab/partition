#ifndef SFC_H
#define SFC_H

void SortMorton(std::vector<double> &coords, uint64_t count, std::vector<uint64_t> &sorted_indices_out);
void ConvertToIntegerCoords(std::vector<double> &coords,uint64_t count, std::vector<uint64_t>& coords_integer_out);
bool MortonCompare(uint64_t x1, uint64_t y1, uint64_t z1, uint64_t x2, uint64_t y2, uint64_t z2 );
#endif