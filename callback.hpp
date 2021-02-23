#ifndef CALLBACK_HPP
#define CALLBACK_HPP

#include <cstddef>
#include <string>
#include <memory>
#include <type_traits>

#include "H5Cpp.h"

#include "fields.hpp"

// base class for all user-defined call-back functions
// User should inherit from this class and override the
// member functions
//
// AFields is a composite type constructed using the AllFields
// template from a GrpFields and a PrtFields type
template<typename AFields>
struct Callback
{
    // T is either AFields::GroupFields or AFields::ParticleFields
    template<typename T>
    class BaseProperties
    {// {{{
        // declare as char so we can advance byte-wise
        // (char is guaranteed to be one byte wide)
        char *data[T::Nfields];

        template<typename Field>
        static constexpr size_t get_field_idx ();

    public :
        // optional offset initializes at offset particles inside the data_in_memory pointers
        BaseProperties (void *data_in_memory[T::Nfields], size_t offset=0UL);

        // advances the internal pointers by the respective strides
        // User should not call this!
        void advance ();

        // retrieves the internal pointers
        // User should not need to call this
        void *operator[] (size_t idx) const;

        // retrieves the position
        // User should not need to call this
        auto coord () const;
        
        // The user function: returns the value of Field.
        // If Field is 1dimensional, the value itself is returned.
        // Otherwise, a pointer to the first element is returned.
        template<typename Field>
        auto get () const;
    };// }}}

    typedef BaseProperties<typename AFields::GroupFields>
        GrpProperties;
    typedef BaseProperties<typename AFields::ParticleFields>
        PrtProperties;

    // these functions write the chunk file name corresponding
    // to chunk_idx into return argument and return true.
    // If no chunk corresponds to chunk_idx, false is returned
    virtual bool grp_chunk (size_t chunk_idx, std::string &fname) const = 0;
    virtual bool prt_chunk (size_t chunk_idx, std::string &fname) const = 0;

    // Where to find the group fields in the hdf5 file
    // (for Illustris: "Group/")
    virtual std::string grp_name () const = 0;

    // Where to find the particle fields in the hdf5 file
    // (for Illustris: "PartType<>/")
    virtual std::string prt_name () const = 0;

    // allows the user to read metadata for each group chunk
    // (thus this function is not marked as constant)
    // This function is required to write the number of groups in this chunk into the
    // return value.
    virtual void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                                size_t &Ngroups) = 0;

    // allows the user to read metadata for each particle chunk
    // (thus this function is not marked as constant)
    // This function is required to write the box size and number of particles in this
    // chunk into the return values.
    virtual void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                                coord_t &Bsize, size_t &Npart) = 0;

    // returns whether the group described by the argument should be
    // considered
    virtual bool grp_select (const GrpProperties &grp) const = 0;

    // this function will be called for all groups for which grp_select
    // returns true.
    // It will be called in order of groups encountered
    // (thus, it is NOT required to be thread safe!)
    // It is not marked as const since it probably stores some data in the child class.
    virtual void grp_action (const GrpProperties &grp) = 0;

    // returns the group radius. Only particles falling within the radius
    // should be considered
    virtual coord_t grp_radius (const GrpProperties &grp) const = 0;

    // this function will be called for all particles for which prt_select
    // returns true.
    // It will not be called in order of particles encountered
    // and should thus be thread safe.
    // grp_idx corresponds to the order in which groups were passed to grp_action.
    // R is the radius from the group center (as defined by the user)
    // It is not marked as const since it probably stores some data in the child class.
    virtual void prt_action (size_t grp_idx, const GrpProperties &grp,
                             const PrtProperties &prt, coord_t Rsq) = 0;
};


// --- Implementation of the BaseProperties struct ---

template<typename AFields>
template<typename T>
Callback<AFields>::BaseProperties<T>::BaseProperties (void *data_in_memory[T::Nfields], size_t offset)
{
    for (size_t ii=0; ii != T::Nfields; ++ii)
        data[ii] = (char *)(data_in_memory[ii]) + offset * T::strides_fcoord[ii];
}

template<typename AFields>
template<typename T>
inline void
Callback<AFields>::BaseProperties<T>::advance ()
{
    for (size_t ii=0; ii != T::Nfields; ++ii)
        data[ii] += T::strides_fcoord[ii];
}

template<typename AFields>
template<typename T>
inline void *
Callback<AFields>::BaseProperties<T>::operator[] (size_t idx) const
{
    return data[idx];
}

template<typename AFields>
template<typename T>
inline auto
Callback<AFields>::BaseProperties<T>::coord () const
{
    return (coord_t *)data[0];
}

template<typename AFields>
template<typename T>
template<typename Field>
constexpr size_t
Callback<AFields>::BaseProperties<T>::get_field_idx ()
{
    static_assert(T::field_type == Field::type);

    if constexpr (T::field_type == FieldTypes::GrpFld)
        return AFields::GroupFields::template idx<Field>;
    else
        return AFields::ParticleFields::template idx<Field>;
}

template<typename AFields>
template<typename T>
template<typename Field>
inline auto
Callback<AFields>::BaseProperties<T>::get () const
{
    constexpr size_t idx = get_field_idx<Field> ();

    if constexpr (Field::dim == 1)
        return (typename Field::value_type) *(typename Field::value_type *)data[idx];
    else
        if constexpr (Field::coord)
        // need to use the global coord_t that this field has been converted to
            return (coord_t *)data[idx];
        else
            return (typename Field::value_type *)data[idx];
}

#endif // CALLBACK_HPP
