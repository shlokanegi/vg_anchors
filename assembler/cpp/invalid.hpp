#ifndef SHASTA_INVALID_HPP
#define SHASTA_INVALID_HPP

// In many contexts, we use invalid<T>
// to indicate a value that is invalid, uninitialized, or unknown.

#include <numeric>

namespace shasta {
    template<class T> static const T invalid = std::numeric_limits<T>::max();
    template<class T> static const T unlimited = std::numeric_limits<T>::max();
}

#endif
