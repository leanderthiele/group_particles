#ifndef WORKSPACE_MEMORY_HPP
#define WORKSPACE_MEMORY_HPP

#include <cassert>
#include <cstdlib>
#include <algorithm>

#include "workspace.hpp"

template<typename GroupFields, typename ParticleFields>
void Workspace::realloc_grp_storage (size_t new_size)
{// {{{
    for (size_t ii=0; ii != GroupFields::Nfields; ++ii)
        grp_properties[ii] = std::realloc(grp_properties[ii],
                                          new_size * GroupFields::strides[ii]);

    grp_radii = (float *)std::realloc(grp_radii, new_size * sizeof(float));
}// }}}

template<typename GroupFields, typename ParticleFields>
void Workspace::realloc_tmp_grp_storage (size_t new_size)
{// {{{
    // we do this with free and malloc instead of realloc to avoid overhead of copying
    // the old data which we don't need anymore
    for (size_t ii=0; ii != GroupFields::Nfields; ++ii)
    {
        if (tmp_grp_properties[ii])
            free(tmp_grp_properties[ii]);
        tmp_grp_properties[ii] = std::malloc(new_size * GroupFields::strides[ii]);
    }
}// }}}

template<typename GroupFields, typename ParticleFields>
void Workspace::realloc_grp_storage_if_necessary ()
{// {{{
    // check if alloc is necessary
    if (Ngrp < alloced_grp)
        return;

    assert(Ngrp == alloced_grp);

    alloced_grp = std::max(128UL, 2 * alloced_grp);
    realloc_grp_storage(alloced_grp);
}// }}}

template<typename GroupFields, typename ParticleFields>
void Workspace::shrink_grp_storage ()
{// {{{
    alloced_grp = Ngrp + 1UL;
    realloc_grp_storage(alloced_grp);
}// }}}

template<typename GroupFields, typename ParticleFields>
Workspace::~Workspace ()
{// {{{
    for (size_t ii=0; ii != GroupFields::Nfields; ++ii)
    {
        if (grp_properties[ii])
            std::free(grp_properties[ii]);
        if (tmp_grp_properties[ii])
            std::free(tmp_grp_properties[ii]);
    }
    if (grp_radii)
        std::free(grp_radii[ii]);
}// }}}

template<typename GroupFields, typename ParticleFields>
Workspace::Workspace ()
{// {{{
    for (size_t ii=0; ii != GroupFields::Nfields; ++ii)
    {
        tmp_grp_properties[ii] = nullptr;
        grp_properties[ii] = nullptr;
    }
    grp_radii = nullptr;
}// }}}



#endif // WORKSPACE_MEMORY_HPP
