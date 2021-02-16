#ifndef GEOM_UTILS_HPP
#define GEOM_UTILS_HPP

#include "fields.hpp"

#include <cmath>
#include <type_traits>

namespace GeomUtils
{

// computes the squared 3D cartesian distance from origin
static inline coord_t
hypotsq (coord_t x, coord_t y, coord_t z)
{// {{{
    return x*x + y*y + z*z;
}// }}}

// computes the signed distance x2-x1, taking into account the periodicity
static inline coord_t
periodic_dist (coord_t x1, coord_t x2, coord_t periodicity)
{// {{{
    // TODO double check this function!
    auto out = x2 - x1;

    if (out > (coord_t)(0.5) * periodicity)
        out -= periodicity;
    else if (out < (coord_t)(-0.5) * periodicity)
        out += periodicity;

    return out;
}// }}}

// computes the unsigned distance |x2-x1|, taking into account the periodicity
static inline coord_t
abs_periodic_dist( coord_t x1, coord_t x2, coord_t periodicity)
{// {{{
    auto out = std::abs(x2 - x1);

    if (out > (coord_t)(0.5) * periodicity)
        out = periodicity - out;

    return out;
}// }}}

// replaces r2 with the signed distance r2-r1, taking into account the periodicity
static inline void
periodic_dist (const coord_t *r1, coord_t *r2, coord_t periodicity)
{// {{{
    for (size_t ii=0; ii != 3; ++ii)
        r2[ii] = periodic_dist(r1[ii], r2[ii], periodicity);
}// }}}

// computes the squared 3D cartesian distance between two points,
// taking into acount the periodicity
static inline coord_t
periodic_hypotsq (const coord_t *r1, const coord_t *r2, coord_t periodicity)
{// {{{
    auto dx0 = abs_periodic_dist(r1[0], r2[0], periodicity);
    auto dx1 = abs_periodic_dist(r1[1], r2[1], periodicity);
    auto dx2 = abs_periodic_dist(r1[2], r2[2], periodicity);

    return hypotsq(dx0, dx1, dx2);
}// }}}

} // namespace GeomUtils


#endif // GEOM_UTILS_HPP
