#ifndef CALLBACK_UTILS_HPP
#define CALLBACK_UTILS_HPP

// This file contains some subclasses of the Callback base class,
// by inheriting a subset from this list the user can achieve a lot
// of functionality out of the box without writing much new code

#include <string>
#include <memory>
#include <numeric_limits>
#include <cstdio>
#include <vector>
#include <functional>

#include "H5Cpp.h"

#include "callbacks.hpp"

namespace callback_utils
{

// Some common ways to have chunks organized
namespace chunk_fmt
{// {{{
    class SingleGrp : virtual public Callback
    {
        const std::string grp_fname;
    public :
        SingleGrp (const std::string &grp_fname_) :
            grp_fname(grp_fname_) { }
        bool grp_chunk (size_t chunk_idx, std::string &fname) const override
        { fname = grp_fname; return (bool)(chunk_idx); }
    };

    class SinglePrt : virtual public Callback
    {
        const std::string prt_fname;
    public :
        SinglePrt (const std::string &prt_fname_) :
            prt_fname(prt_fname_) { }
        bool prt_chunk (size_t chunk_idx, std::string &fname) const override
        { fname = prt_fname; return (bool)(chunk_idx); }
    };

    class MultiGrp : virtual public Callback
    {
        const std::string grp_fname;
        const size_t max_idx;
    public :
        MultiGrp (const std::string &grp_fname_, size_t max_idx_) :
            grp_fname(grp_fname_), max_idx(max_idx_) { }
        bool grp_chunk (size_t chunk_idx, std::string &fname) const override
        { char buf[grp_fname.size()+10]; std::sprintf(buf, grp_fname.c_str(), chunk_idx);
          fname = std::string(buf); return chunk_idx <= max_idx; }
    };

    class MultiPrt : virtual public Callback
    {
        const std::string prt_fname;
        const size_t max_idx;
    public :
        MultiPrt (const std::string &prt_fname_, size_t max_idx_) :
            prt_fname(prt_fname_), max_idx(max_idx_) { }
        bool prt_chunk (size_t chunk_idx, std::string &fname) const override
        { char buf[prt_fname.size()+10]; std::sprintf(buf, prt_fname.c_str(), chunk_idx);
          fname = std::string(buf); return chunk_idx <= max_idx; }
    };

    struct Single : virtual public Callback, virtual public SingleGrp, virtual public SinglePrt
    {
        Single (const std::string &grp_fname,
                const std::string &prt_fname) :
            SingleGrp(grp_fname), SinglePrt(prt_fname) { }
    };

    struct Multi : virtual public Callback, virtual public MultiGrp, virtual public MultiPrt
    {
        Multi (const std::string &grp_fname, size_t grp_max_idx,
               const std::string &prt_fname, size_t prt_max_idx) :
            MultiGrp(grp_fname, grp_max_idx), MultiPrt(prt_fname, prt_max_idx) { }
    };
} // namespace chunk_fmt }}}

// Some helpers for the Illustris file format
namespace illustris
{// {{{
    template<uint8_t PartType>
    struct Naming : virtual public Callback
    {
    public :
        Naming (uint8_t part_type_) :
            part_type(part_type_) { }
        std::string grp_name () const override
        { return "Group/"; }
        std::string prt_name () const override
        { return "PartType"+std::to_string(PartType)+"/"; }
    };

    template<uint8_t PartType>
    class Meta : virtual public Callback
    {
        template<typename TH5, typename Trequ>
        Trequ read_scalar_attr (H5::Group &header, const std::string &name) const
        {// {{{
            TH5 out;
            auto attr = header.openAttribute(name);
            assert (attr.getDataType().getSize() == sizeof(TH5));
            attr.read(attr.getDataType(), &out);
            return (Trequ)(out);
        }// }}}
        template<typename TH5, typename Trequ>
        Trequ read_vector_attr (H5::Group &header, const std::string &name, size_t idx) const
        {// {{{
            auto attr = header.openAttribute(name);
            auto aspace = attr.getSpace();
            hsize_t dim_lenghts[16];
            auto Ndims = aspace.getSimpleExtentDims(dim_lenghts);
            assert(Ndims == 1)
            assert(dim_lenghts[0] > idx);
            TH5 out[dim_lenghts[0]];
            attr.read(attr.getDataType(), out);
            return (Trequ)(out[idx]);
        }// }}}
    protected :
        // user can override these functions if it is desired to do more than the required
        // job with the hdf5 file metadata
        // These functions are called by the exposed functions.
        virtual void read_grp_meta_custom (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr)
        { return; }
        virtual void read_prt_meta_custom (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr)
        { return; }
    public :
        // TODO figure out what types the Illustris entries have!
        void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr, size_t &Ngroups) override
        {
            auto header = fptr->openGroup("/Header");
            Ngroups = read_scalar_attr</*TODO*/,size_t>(header, "Ngroups_ThisFile");
            read_grp_meta_custom(chunk_idx, fptr);
        }
        void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr, float &Bsize, size_t &Npart) override
        {
            auto header = fptr->openGroup("/Header");
            Bsize = read_scalar_attr</*TODO*/,float>(header, "BoxSize");
            Npart = read_vector_attr</*TODO*/,size_t>(header, "NumPart_ThisFile", PartType);
            read_prt_meta_custom(chunk_idx, fptr);
        }
    };

    template<uint8_t PartType>
    struct Conventional : virtual public Callback, virtual public Naming<PartType>, virtual public Meta<PartType>
    { };
}// namespace illustris }}}

// Some common cases for the selection functions
namespace select
{// {{{
    template<size_t Midx, typename TMass=float>
    class MassWindow : virtual public Callback
    {
        Tmass Mmin, Mmax;
    public :
        MassWindow (TMass Mmin_, Tmass Mmax_) :
            Mmin(Mmin_), Mmax(Mmax_) { }
        bool grp_select (void **grp_properties) const override
        { TMass M = *(TMass *)(grp_properties[Midx]); return M>Mmin && M<Mmax; }
    };
    
    template<size_t Midx, typename TMass=float>
    struct MassHighCutoff : virtual public Callback, virtual public MassWindow<Midx,TMass>
    {
        MassHighCutoff (TMass Mmin) :
            MassWindow<Midx,TMass>(Mmin, std::numeric_limits<TMass>::max()) { }
    };
    
    template<size_t Midx, typename TMass=float>
    struct MassLowCutoff : virtual public Callback, virtual public MassWindow<Midx,TMass>
    {
        MassLowCutoff (TMass Mmax) :
            MassWindow<Midx,TMass>(std::numeric_limits<TMass>::min(), Mmax) { }
    };

    struct AllParticles : virtual public Callback
    {
        bool prt_select (void **grp_properties, void **prt_properties, float R) const override
        { return true; }
    }
}// namespace select }}}

// Some common cases for the radius function
namespace radius
{// {{{
    template<size_t Ridx, typename TR=float>
    class Simple
    {
        TR scaling;
    public :
        Simple (TR scaling_) : scaling(scaling_) { }
        float grp_radius (void **grp_properties)
        { return scaling * *(TR *)(grp_properties[Ridx]); }
    };
}// namespace radius }}}

namespace actions
{// {{{
    class MultiGrp : virtual public Callback
    {
    protected :
        std::vector<std::function<void(void **)>> grp_actions;
        void add_grp_action (std::function<void(void **)> fct)
        { grp_actions.push_back(fct); }
    public :
        void grp_action (void **grp_properties) override
        { for (auto fct : grp_actions) fct(grp_properties); }
    };

    template<typename Tdata, typename Tfromprt=Tdata>
    class StoreHomogeneous : virtual public Callback, virtual public MultiGrp
    {
        std::shared_ptr<std::vector<Tdata>> data;
        void enlarge_data (void **grp_properties)
        { data->resize(data->size() + 1); }
    protected :
        // how to reduce data about a particle to a Tfromprt object
        virtual Tfromprt prt_reduce (size_t grp_idx, void **grp_properties, void **prt_properties, float R) const = 0;
        // how to combine the reduced particle data with the stored data
        virtual void prt_combine (size_t grp_idx, Tdata &data_item, Tfromprt prt_item) = 0;
    public :
        StoreHomogeneous (std::shared_ptr<std::vector<Tdata>> data_) :
            data(data_)
        { add_grp_action(enlarge_data); }
        void prt_action (size_t grp_idx, void **grp_properties, void **prt_properties, float R) override
        {
            auto prt_item = prt_reduce(grp_idx, grp_properties, prt_properties, R);
            prt_combine(grp_idx, data[grp_idx], prt_item);
        }
    };

    template<typename Tdata>
    class StoreGrpProperties : virtual public Callback, virtual public MultiGrp
    {
        std::shared_ptr<std::vector<Tdata>> data:
        void append_data (void **grp_properties)
        { data->push_back(grp_reduce(grp_properties)); }
    protected :
        virtual Tdata grp_reduce (void **grp_properties) = 0;
    public :
        StoreGrpProperties (std::shared_ptr<std::vector<Tdata>> data_) :
            data(data_)
        { add_grp_action(append_data); }
    };

}// namespace actions }}}

} // namespace callback_utils

#endif // CALLBACK_UTILS_HPP
