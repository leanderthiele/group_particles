#ifndef WORKSPACE_MEMORY_HPP
#define WORKSPACE_MEMORY_HPP

#include <cassert>
#include <cstdlib>
#include <algorithm>

#include "workspace.hpp"


template<typename AFields>
Workspace<AFields>::Workspace (Callback<AFields> &callback_) :
    callback(callback_)
{// {{{
    for (size_t ii=0; ii != AFields::GroupFields::Nfields; ++ii)
    {
        tmp_grp_properties[ii] = nullptr;
        grp_properties[ii] = nullptr;
    }
    for (size_t ii=0; ii != AFields::ParticleFields::Nfields; ++ii)
        tmp_prt_properties[ii] = nullptr;
    grp_radii = nullptr;
}// }}}

template<typename AFields>
Workspace<AFields>::~Workspace ()
{// {{{
    for (size_t ii=0; ii != AFields::GroupFields::Nfields; ++ii)
    {
        if (grp_properties[ii])
            std::free(grp_properties[ii]);
        if (tmp_grp_properties[ii])
            std::free(tmp_grp_properties[ii]);
    }

    for (size_t ii=0; ii != AFields::ParticleFields::Nfields; ++ii)
        if (tmp_prt_properties[ii])
            std::free(tmp_prt_properties[ii]);

    if (grp_radii)
        std::free(grp_radii);
}// }}}

template<typename AFields>
void Workspace<AFields>::realloc_grp_storage (size_t new_size)
{// {{{
    for (size_t ii=0; ii != AFields::GroupFields::Nfields; ++ii)
        grp_properties[ii] = std::realloc(grp_properties[ii],
                                          new_size * AFields::GroupFields::strides[ii]);

    grp_radii = (float *)std::realloc(grp_radii, new_size * sizeof(float));
}// }}}

template<typename AFields>
template<typename T>
void Workspace<AFields>::realloc_tmp_storage (size_t new_size, void **buf)
{// {{{
    for (size_t ii=0; ii != T::Nfields; ++ii)
    {
        if (buf[ii]) std::free(buf[ii]);
        buf[ii] = std::malloc(new_size * T::strides[ii]);
    }
}// }}}

template<typename AFields>
void Workspace<AFields>::realloc_grp_storage_if_necessary ()
{// {{{
    // check if alloc is necessary
    if (Ngrp < alloced_grp)
        return;

    assert(Ngrp == alloced_grp);

    alloced_grp = std::max(128UL, 2 * alloced_grp);
    realloc_grp_storage(alloced_grp);
}// }}}

template<typename AFields>
void Workspace<AFields>::shrink_grp_storage ()
{// {{{
    alloced_grp = Ngrp + 1UL;
    realloc_grp_storage(alloced_grp);
}// }}}

template<typename AFields>
template<typename T>
inline void
Workspace<AFields>::collect_properties
    (void **dest, void **src, size_t idx)
{// {{{
    for (size_t ii=0; ii != T::Nfields; ++ii)
        dest[ii] = (char *)(src[ii]) + idx * T::strides[ii];
}// }}}


#endif // WORKSPACE_MEMORY_HPP
