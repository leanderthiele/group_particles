#ifndef READ_FIELDS_HPP
#define READ_FIELDS_HPP

#include <memory>
#include <string>
#include <cstdlib>

#include "fields.hpp"
#include "callback.hpp"

// It is assumed that data is already allocated storage of the required size
void
read_field (std::shared_ptr<H5::H5File> fptr, const std::string &name,
            // these are only for debugging purposes
            size_t element_size, size_t Nitems, size_t dim,
            void * data)
{// {{{
    auto dset   = fptr->openDataSet(name);
    auto dspace = dset.getSpace();

    hsize_t dim_lengths[16];
    auto Ndims   = dspace.getSimpleExtentDims(dim_lengths);
    auto Npoints = dspace.getSimpleExtentNpoints();
    auto Dtype   = dset.getDataType();

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
template<FieldTypes field_type, typename Fields>
void
read_fields (const Callback &callback,
             std::shared_ptr<H5::H5File> fptr, size_t Nitems, void **data)
{
    // where to find our data sets in the hdf5 file
    std::string name_prefix;
    if constexpr (field_type == FieldTypes::GrpField)
        name_prefix = callback.grp_name();
    else
        name_prefix = callback.prt_name();

    // loop over the fields
    for (size_t ii=0; ii != Fields::Nfields; ++ii)
    {
        // allocate storage TODO get rid of this
        data[ii] = std::malloc(Nitems * Fields::strides[ii]);

        // read from disk
        read_field(fptr, name_prefix + Fields::names[ii],
                   Fields::sizes[ii], Nitems, Fields::dim[ii],
                   data[ii]);
    }
}

#endif // READ_FIELDS_HPP
