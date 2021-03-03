#ifndef READ_FIELDS_HPP
#define READ_FIELDS_HPP

#include <cassert>
#include <memory>
#include <string>
#include <type_traits>
#include <cstdlib>
#include <cstdio>

#include "fields.hpp"
#include "callback.hpp"

namespace hdf5Utils {

template<typename TH5, typename Trequ>
static Trequ
read_scalar_attr (H5::Group &header, const std::string &name)
{// {{{
    TH5 out;
    auto attr = header.openAttribute(name);
    assert (attr.getDataType().getSize() == sizeof(TH5));
    attr.read(attr.getDataType(), &out);
    return (Trequ)(out);
}// }}}

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

// It is assumed that data is already allocated storage of the required size
static void
read_field (std::shared_ptr<H5::H5File> fptr, const std::string &name,
            // these are only for debugging purposes
            size_t element_size, size_t Nitems, size_t dim,
            void * data)
{// {{{
    auto dset   = fptr->openDataSet(name);
    auto dspace = dset.getSpace();

    hsize_t dim_lengths[16];
    auto Ndims  = dspace.getSimpleExtentDims(dim_lengths);
    auto Dtype  = dset.getDataType();

    // some easy consistency checks
    assert(Dtype.getSize() == element_size);
    assert(Nitems == dim_lengths[0]);
    assert((Ndims==1 && dim==1) || (Ndims==2 && dim_lengths[1]==dim));

    // read into memory
    auto memspace = H5::DataSpace(Ndims, dim_lengths);
    dset.read(data, Dtype, memspace, dspace);
}// }}}

// it is assumed that data is already of the correct size
// and the individual pointers are already allocated
// T is one of GroupFields, ParticleFields
template<typename AFields, typename T>
static void
read_fields (const Callback<AFields> &callback,
             std::shared_ptr<H5::H5File> fptr, size_t Nitems, void **data)
{// {{{
    // where to find our data sets in the hdf5 file
    std::string name_prefix;
    if constexpr (T::field_type == FieldTypes::GrpFld)
        name_prefix = callback.grp_name();
    else
        name_prefix = callback.prt_name();

    // loop over the fields
    for (size_t ii=0; ii != T::Nfields; ++ii)
        // read from disk
        read_field(fptr, name_prefix + T::names[ii],
                   T::sizes[ii], Nitems, T::dims[ii],
                   data[ii]);
}// }}}

} // namespace hdf5Utils

#endif // READ_FIELDS_HPP
