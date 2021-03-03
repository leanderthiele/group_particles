#ifndef WORKSPACE_META_INIT_HPP
#define WORKSPACE_META_INIT_HPP

#include "H5Cpp.h"

#include "workspace.hpp"

namespace grp_prt_detail {

template<typename AFields>
void
Workspace<AFields>::meta_init ()
{
    std::string fname;

    // look at the first group chunk
    callback.grp_chunk(0, fname);
    auto fptr_grp = std::make_shared<H5::H5File>(fname, H5F_ACC_RDONLY);
    callback.read_grp_meta_init(fptr_grp);
    fptr_grp->close();

    // look at the first group chunk
    callback.prt_chunk(0, fname);
    auto fptr_prt = std::make_shared<H5::H5File>(fname, H5F_ACC_RDONLY);
    callback.read_prt_meta_init(fptr_prt);
    fptr_prt->close();
}

} // namespace grp_prt_detail

#endif // WORKSPACE_META_INIT_HPP
