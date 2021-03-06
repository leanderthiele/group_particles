#ifndef GEOM_UTILS_HPP
#define GEOM_UTILS_HPP

#include "fields.hpp"

#include <cmath>
#include <array>
#include <type_traits>

namespace grp_prt_detail {

namespace GeomUtils {

// wraps idx with period N
static inline size_t
periodic_idx (int idx, int N)
{// {{{
    return (N + idx%N) % N;
}// }}}

// computes the squared 3D cartesian distance from origin
static inline coord_t
hypotsq (coord_t x, coord_t y, coord_t z)
{// {{{
    return x*x + y*y + z*z;
}// }}}

static inline coord_t
hypotsq (coord_t *r)
{// {{{
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}// }}}

// computes the unsigned distance |x2-x1|, taking into account the periodicity
static inline coord_t
abs_periodic_dist (coord_t x1, coord_t x2, coord_t periodicity)
{// {{{
    auto out = std::abs(x2 - x1);

    if (out > (coord_t)(0.5) * periodicity)
        out = periodicity - out;

    return out;
}// }}}

// computes the signed distance x2-x1, taking into account the periodicity
static inline coord_t
periodic_dist (coord_t x1, coord_t x2, coord_t periodicity)
{
    auto out = x2 - x1;

    if (out > (coord_t)(0.5) * periodicity)
        out -= periodicity;
    else if (out < (coord_t)(-0.5) * periodicity)
        out += periodicity;

    return out;
}

// equal action to previous function, but uses precomputed periodicity
// (and does not use absolute value)
__attribute__((hot))
static inline coord_t
periodic_dist_whint (coord_t x1, coord_t x2, coord_t periodicity, int periodic_to_add)
{// {{{
    coord_t dx = x2 - x1;

    return dx + periodicity * (coord_t)periodic_to_add;
}// }}}

// computes the squared 3D cartesian distance between two points,
// taking into acount the periodicity
static inline coord_t
periodic_hypotsq (const coord_t *__restrict__ r1, const coord_t *__restrict__ r2,
                  coord_t periodicity)
{// {{{
    auto dx0 = abs_periodic_dist(r1[0], r2[0], periodicity);
    auto dx1 = abs_periodic_dist(r1[1], r2[1], periodicity);
    auto dx2 = abs_periodic_dist(r1[2], r2[2], periodicity);

    return hypotsq(dx0, dx1, dx2);
}// }}}

// same as the above, but with precomputed periodicity hint
__attribute__((hot))
static inline coord_t
periodic_hypotsq (const coord_t *__restrict__ r1, const coord_t *__restrict__ r2,
                  coord_t periodicity, const std::array<int,3> &periodic_to_add)
{// {{{
    coord_t dx[3];
    
    for (size_t ii=0; ii != 3; ++ii)
        dx[ii] = periodic_dist_whint(r1[ii], r2[ii], periodicity, periodic_to_add[ii]);

    return hypotsq(dx);
}// }}}

} // namespace GeomUtils

} // namespace grp_prt_detail


#endif // GEOM_UTILS_HPP
