#include "MinCost.h"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif


template class MinCost<int,int>;
template class MinCost<int,double>;
template class MinCost<int,int16_t>;
template class MinCost<int,int8_t>;

template class DualMinCost<int>;
template class DualMinCost<double>;
template class DualMinCost<int16_t>;
template class DualMinCost<int8_t>;
