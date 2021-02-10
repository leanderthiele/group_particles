#ifndef WORKSPACE_HPP
#define WORKSPACE_HPP

#include "callback.hpp"
#include "fields.hpp"

template<typename GroupFields, typename ParticleFields>
class Workspace
{
    // metadata we need to store permanently
    float Bsize;

    // data we need to store permanently
    // (acoording to user-defined selection and radius calculation)
    size_t Ngrp = 0UL;
    size_t alloced_grp = 0UL;
    void *grp_properties[GroupFields::Nfields];
    float *grp_radii;

    // temporary buffers
    void *tmp_grp_properties[GroupFields::Nfields];
    void *tmp_prt_properties[ParticleFields::Nfields];

    void realloc_grp_storage (size_t new_size);

    template<typename Fields>
    void realloc_tmp_storage (size_t new_size, void **buf);

    // we can also use this function for the initial malloc
    void realloc_grp_storage_if_necessary ();

    void shrink_grp_storage ();

public :
    void grp_loop (Callback &callback);
    void prt_loop (Callback &callback);

    Workspace ();

    ~Workspace ();

    // Note : we do not implement copy and move constructors,
    //        despite explicit memory management.
    //        The caller should simply not use these.
};

#endif // WORKSPACE_HPP
