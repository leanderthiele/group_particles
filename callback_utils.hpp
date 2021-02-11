#ifndef CALLBACK_UTILS_HPP
#define CALLBACK_UTILS_HPP

// This file contains some subclasses of the Callback base class,
// by inheriting a subset from this list the user can achieve a lot
// of functionality out of the box without writing much new code

#include <cassert>
#include <string>
#include <memory>
#include <limits>
#include <cstdio>
#include <vector>
#include <functional>
#include <utility>
#include <cstdio>

#include "H5Cpp.h"

#include "callback.hpp"

namespace CallbackUtils
{

// Some common ways to have chunks organized
namespace chunk_fmt
{// {{{
    class SingleGrp : virtual public Callback
    {
        const std::string grp_fname;
    public :
        SingleGrp (const std::string &grp_fname_) :
            grp_fname(grp_fname_)
        {
            #ifndef NDEBUG
            std::fprintf(stderr, "Initialized CallbackUtils::chunk_fmt::SingleGrp with grp_fname=%s\n", grp_fname.c_str());
            #endif
        }
        bool grp_chunk (size_t chunk_idx, std::string &fname) const override
        { fname = grp_fname; return (bool)(chunk_idx); }
    };

    class SinglePrt : virtual public Callback
    {
        const std::string prt_fname;
    public :
        SinglePrt (const std::string &prt_fname_) :
            prt_fname(prt_fname_)
        {
            #ifndef NDEBUG
            std::fprintf(stderr, "Initialized CallbackUtils::chunk_fmt::SinglePrt with prt_fname=%s\n", prt_fname.c_str());
            #endif
        }
        bool prt_chunk (size_t chunk_idx, std::string &fname) const override
        { fname = prt_fname; return (bool)(chunk_idx); }
    };

    class MultiGrp : virtual public Callback
    {
        const std::string grp_fname;
        const size_t max_idx;
    public :
        MultiGrp (const std::string &grp_fname_, size_t max_idx_) :
            grp_fname(grp_fname_), max_idx(max_idx_)
        {
            #ifndef NDEBUG
            std::fprintf(stderr, "Initialized CallbackUtils::chunk_fmt::MultiGrp with grp_fname=%s and max_idx=%lu\n", grp_fname.c_str(), max_idx);
            #endif
        }
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
            prt_fname(prt_fname_), max_idx(max_idx_)
        {
            #ifndef NDEBUG
            std::fprintf(stderr, "Initialized CallbackUtils::chunk_fmt::MultiPrt with prt_fname=%s and max_idx=%lu\n", prt_fname.c_str(), max_idx);
            #endif
        }
        bool prt_chunk (size_t chunk_idx, std::string &fname) const override
        { char buf[prt_fname.size()+10]; std::sprintf(buf, prt_fname.c_str(), chunk_idx);
          fname = std::string(buf); return chunk_idx <= max_idx; }
    };

    struct Single : virtual public Callback, public SingleGrp, public SinglePrt
    {
        Single (const std::string &grp_fname,
                const std::string &prt_fname) :
            SingleGrp(grp_fname), SinglePrt(prt_fname) { }
    };

    struct Multi : virtual public Callback, public MultiGrp, public MultiPrt
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
            assert(Ndims == 1);
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
        void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr, size_t &Ngroups) override
        {
            auto header = fptr->openGroup("/Header");
            Ngroups = read_scalar_attr<int32_t,size_t>(header, "Ngroups_ThisFile");
            read_grp_meta_custom(chunk_idx, fptr);
        }
        void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr, float &Bsize, size_t &Npart) override
        {
            auto header = fptr->openGroup("/Header");
            Bsize = read_scalar_attr<double,float>(header, "BoxSize");
            Npart = read_vector_attr<int32_t,size_t>(header, "NumPart_ThisFile", PartType);
            read_prt_meta_custom(chunk_idx, fptr);
        }
    };

    template<uint8_t PartType>
    struct Conventional : virtual public Callback, public Naming<PartType>, public Meta<PartType>
    { };
}// namespace illustris }}}

// Some common cases for the selection functions
namespace select
{// {{{
    template<size_t Midx, typename TMass=float>
    class GrpMassWindow : virtual public Callback
    {
        TMass Mmin, Mmax;
    public :
        GrpMassWindow (TMass Mmin_, TMass Mmax_) :
            Mmin(Mmin_), Mmax(Mmax_)
        {
            #ifndef NDEBUG
            std::fprintf(stderr, "Initialized CallbackUtils::select::GrpMassWindow with Mmin=%f and Mmax=%f.\n", Mmin, Mmax);
            #endif
        }
        bool grp_select (void **grp_properties) const override
        { TMass M = *(TMass *)(grp_properties[Midx]); return M>Mmin && M<Mmax; }
    };
    
    template<size_t Midx, typename TMass=float>
    struct GrpMassLowCutoff : virtual public Callback, public GrpMassWindow<Midx,TMass>
    {
        GrpMassLowCutoff (TMass Mmin) :
            GrpMassWindow<Midx,TMass>(Mmin, std::numeric_limits<TMass>::max()) { }
    };
    
    template<size_t Midx, typename TMass=float>
    struct GrpMassHighCutoff : virtual public Callback, public GrpMassWindow<Midx,TMass>
    {
        GrpMassHighCutoff (TMass Mmax) :
            GrpMassWindow<Midx,TMass>(std::numeric_limits<TMass>::min(), Mmax) { }
    };

    struct PrtAll : virtual public Callback
    {
        bool prt_select (void **grp_properties, void **prt_properties, float R) const override
        { return true; }
    };
}// namespace select }}}

// Some common cases for the radius function
namespace radius
{// {{{
    template<size_t Ridx, typename TR=float>
    class Simple : virtual public Callback
    {
        TR scaling;
    public :
        Simple (TR scaling_) : scaling(scaling_)
        {
            #ifndef NDEBUG
            std::fprintf(stderr, "Initialized CallbackUtils::radius::Simple with scaling=%f.\n", scaling);
            #endif
        }
        float grp_radius (void **grp_properties) const override
        { return scaling * *(TR *)(grp_properties[Ridx]); }
    };
}// namespace radius }}}

namespace actions
{// {{{
    class MultiGrpAction : virtual public Callback
    {
    protected :
        std::vector<std::pair<void *,std::function<void(void *, void **)>>> grp_actions;
    public :
        void grp_action (void **grp_properties) override
        { for (auto fct : grp_actions) fct.second(fct.first, grp_properties); }
    };

    template<typename Tdata, typename Tfromprt=Tdata>
    class StoreHomogeneous : virtual public Callback, virtual public MultiGrpAction
    {
        std::vector<Tdata> *data;
        static void enlarge_data (void *obj, void **grp_properties)
        {
            StoreHomogeneous *p = (StoreHomogeneous *)obj;
            p->data->resize(p->data->size() + 1);
        }
    protected :
        // how to reduce data about a particle to a Tfromprt object
        virtual Tfromprt prt_reduce (size_t grp_idx, void **grp_properties, void **prt_properties, float R) const = 0;
        // how to combine the reduced particle data with the stored data
        virtual void prt_combine (size_t grp_idx, Tdata &data_item, Tfromprt prt_item) const = 0;
    public :
        StoreHomogeneous (std::vector<Tdata> *data_) :
            data(data_)
        { MultiGrpAction::grp_actions.push_back(std::make_pair(this, enlarge_data)); }
        void prt_action (size_t grp_idx, void **grp_properties, void **prt_properties, float R) override
        {
            auto prt_item = prt_reduce(grp_idx, grp_properties, prt_properties, R);
            prt_combine(grp_idx, (*data)[grp_idx], prt_item);
        }
    };

    template<typename Tdata>
    class StoreGrpProperties : virtual public Callback, virtual public MultiGrpAction
    {
        std::vector<Tdata> *data;
        static void append_data (void *obj, void **grp_properties)
        {
            StoreGrpProperties *p = (StoreGrpProperties *)obj;
            p->data->push_back(p->grp_reduce(grp_properties));
        }
    protected :
        virtual Tdata grp_reduce (void **grp_properties) const = 0;
    public :
        StoreGrpProperties (std::vector<Tdata> *data_) :
            data(data_)
        { MultiGrpAction::grp_actions.push_back(std::make_pair(this, append_data)); }
    };

}// namespace actions }}}

} // namespace CallbackUtils

#endif // CALLBACK_UTILS_HPP
