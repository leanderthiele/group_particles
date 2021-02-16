#ifndef GEOM_UTILS_HPP
#define GEOM_UTILS_HPP

namespace GeomUtils
{

// computes the squared 3D cartesian distance from origin
template<typename T>
static inline T
hypotsq (T x, T y, T z)
{
    return x*x + y*y + z*z;
}

// computes the signed distance x2-x1, taking into account the periodicity
template<typename Tout, typename T1, typename T2, typename T3>
static inline Tout
periodic_dist (T1 x1, T2 x2, T3 periodicity)
{
    // TODO double check this function!
    Tout out = x2 - x1;

    if (out > (T3)(0.5) * periodicity)
        out -= periodicity;
    else if (out < (T3)(-0.5) * periodicity)
        out += periodicity;

    return out;
}

// computes the squared 3D cartesian distance between two points,
// taking into acount the periodicity
template<typename Tout, typename T1, typename T2, typename T3>
static inline Tout
periodic_hypotsq (T1 *r1, T2 *r2, T3 periodicity)
{
    auto dx0 = periodic_dist<Tout>(r1[0], r2[0], periodicity);
    auto dx1 = periodic_dist<Tout>(r1[1], r2[1], periodicity);
    auto dx2 = periodic_dist<Tout>(r1[2], r2[2], periodicity);

    return hypotsq(dx0, dx1, dx2);
}

} // namespace GeomUtils


#endif // GEOM_UTILS_HPP
