/*! @file hdf5_utils.hpp
 *
 * @brief some utility functions to read attributes from hdf5 files
 */

#ifndef READ_FIELDS_HPP
#define READ_FIELDS_HPP

#include <cassert>
#include <memory>
#include <string>
#include <cstdlib>

#include "H5Cpp.h"

/*! @brief some utility functions to read attributes from hdf5 files */
namespace hdf5Utils {

/*! @brief reads scalar attribute from an hdf5 group
 *
 * @tparam TH5      the attribute's type in the hdf5 file
 * @tparam Trequ    the type the attribute should be returned as
 *
 * @param header    the group the attribute is attached to
 * @param name      the attribute's name
 *
 * @return the attribute
 */
template<typename TH5, typename Trequ, typename Tobj=H5::Group>
static Trequ
read_scalar_attr (Tobj &header, const std::string &name)
{// {{{
    TH5 out;
    auto attr = header.openAttribute(name);
    assert (attr.getDataType().getSize() == sizeof(TH5));
    attr.read(attr.getDataType(), &out);
    return (Trequ)(out);
}// }}}

/*! @brief reads element of vector attribute from an hdf5 group
 *
 * @tparam TH5      the attribute's element type in the hdf5 file
 * @tparam Trequ    the type the attribute element should be returned as
 *
 * @param header    the group the vector attribute is attached to
 * @param name      the attribute's name
 * @param idx       element index in the vector attribute
 *
 * @return the attribute element
 */
template<typename TH5, typename Trequ>
static Trequ
read_vector_attr (H5::Group &header, const std::string &name, size_t idx)
{// {{{
    auto attr = header.openAttribute(name);
    auto aspace = attr.getSpace();
    hsize_t dim_lenghts[16];
    auto Ndims = aspace.getSimpleExtentDims(dim_lenghts);
    assert(Ndims == 1);
    assert(dim_lenghts[0] > idx);
    assert(attr.getDataType().getSize() == sizeof(TH5));
    TH5 out[dim_lenghts[0]];
    attr.read(attr.getDataType(), out);
    return (Trequ)(out[idx]);
}// }}}

/*! @brief reads entire vector attribute from an hdf5 group
 *
 * @tparam TH5      the attribute's element type in the hdf5 file
 * @tparam Trequ    the type the attribute's elements should be returned as
 *
 * @param header    the group the vector attribute is attached to
 * @param name      the attribute's name
 * @param out       output buffer (sufficient space must be already allocated)
 *
 * @return number of elements in the vector attribute
 */
template<typename TH5, typename Trequ>
static size_t
read_vector_attr (H5::Group &header, const std::string &name, Trequ *out)
{// {{{
    auto attr = header.openAttribute(name);
    auto aspace = attr.getSpace();
    hsize_t dim_lenghts[16];
    auto Ndims = aspace.getSimpleExtentDims(dim_lenghts);
    assert(Ndims == 1);
    assert(attr.getDataType().getSize() == sizeof(TH5));
    
    if constexpr (std::is_same_v<TH5, Trequ>)
        attr.read(attr.getDataType(), out);
    else
    {
        TH5 out1[dim_lenghts[0]];
        attr.read(attr.getDataType(), out1);
        for (size_t ii=0; ii != dim_lenghts[0]; ++ii)
            out[ii] = (Trequ)out1[ii];
    }

    return dim_lenghts[0];
}// }}}

} // namespace hdf5Utils

#endif // READ_FIELDS_HPP
