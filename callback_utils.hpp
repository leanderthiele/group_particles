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
#include <type_traits>

#include "H5Cpp.h"

#include "callback.hpp"
#include "fields.hpp"

namespace CallbackUtils
{

// Some common ways to have chunks organized
namespace chunk
{// {{{
    template<typename AFields>
    class SingleGrp : virtual public Callback<AFields>
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
        { fname = grp_fname; return !chunk_idx; }
    };

    template<typename AFields>
    class SinglePrt : virtual public Callback<AFields>
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
        { fname = prt_fname; return !chunk_idx; }
    };

    template<typename AFields>
    class MultiGrp : virtual public Callback<AFields>
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

    template<typename AFields>
    class MultiPrt : virtual public Callback<AFields>
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

    template<typename AFields>
    struct Single : virtual public Callback<AFields>, public SingleGrp<AFields>, public SinglePrt<AFields>
    {
        Single (const std::string &grp_fname,
                const std::string &prt_fname) :
            SingleGrp<AFields>(grp_fname), SinglePrt<AFields>(prt_fname) { }
    };

    template<typename AFields>
    struct Multi : virtual public Callback<AFields>, public MultiGrp<AFields>, public MultiPrt<AFields>
    {
        Multi (const std::string &grp_fname, size_t grp_max_idx,
               const std::string &prt_fname, size_t prt_max_idx) :
            MultiGrp<AFields>(grp_fname, grp_max_idx), MultiPrt<AFields>(prt_fname, prt_max_idx) { }
    };
} // namespace chunk }}}

// Some common ways to have hdf5 file contents organized
namespace name
{// {{{
    template<typename AFields, uint8_t PartType>
    struct Illustris : virtual public Callback<AFields>
    {
    public :
        std::string grp_name () const override
        { return "Group/"; }
        std::string prt_name () const override
        { return "PartType"+std::to_string(PartType)+"/"; }
    };
} // namespace name }}}

// Some common ways to get meta data
namespace meta
{// {{{
    template<typename AFields, uint8_t PartType>
    class Illustris : virtual public Callback<AFields>
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
}// namespace meta }}}

// Some common cases for the selection functions
namespace select
{// {{{
    template<typename AFields>
    struct GrpAll : virtual public Callback<AFields>
    {
        bool grp_select (void **grp_properties)
        { return true; }
    };

    template<typename AFields, typename MassField>
    class GrpMassWindow : virtual public Callback<AFields>
    {
        static_assert(MassField::dim == 1);
        static_assert(MassField::type == FieldTypes::GrpFld);
        static_assert(std::is_floating_point_v<typename MassField::value_type>);
        typename MassField::value_type Mmin, Mmax;
    public :
        GrpMassWindow (typename MassField::value_type Mmin_, typename MassField::value_type Mmax_) :
            Mmin(Mmin_), Mmax(Mmax_)
        {
            #ifndef NDEBUG
            std::fprintf(stderr, "Initialized CallbackUtils::select::GrpMassWindow with "
                                 "Mmin=%.3e and Mmax=%.3e.\n", Mmin, Mmax);
            #endif
        }
        bool grp_select (void **grp_properties) const override
        {
            auto M = Callback<AFields>::template get_property<MassField>(grp_properties);
            return M>Mmin && M<Mmax;
        }
    };
    
    template<typename AFields, typename MassField>
    struct GrpMassLowCutoff : virtual public Callback<AFields>, public GrpMassWindow<AFields, MassField>
    {
        GrpMassLowCutoff (typename MassField::value_type Mmin) :
            GrpMassWindow<AFields,MassField>(Mmin, std::numeric_limits<typename MassField::value_type>::max()) { }
    };
    
    template<typename AFields, typename MassField>
    struct GrpMassHighCutoff : virtual public Callback<AFields>, public GrpMassWindow<AFields,MassField>
    {
        GrpMassHighCutoff (typename MassField::value_type Mmax) :
            GrpMassWindow<AFields,MassField>(std::numeric_limits<typename MassField::value_type>::min(), Mmax) { }
    };

    template<typename AFields>
    struct PrtAll : virtual public Callback<AFields>
    {
        bool prt_select (size_t grp_idx, void **grp_properties, void **prt_properties, float R) const override
        { return true; }
    };
}// namespace select }}}

// Some common cases for the radius function
namespace radius
{// {{{
    template<typename AFields, typename RField>
    class Simple : virtual public Callback<AFields>
    {
        static_assert(RField::dim == 1);
        static_assert(RField::type == FieldTypes::GrpFld);
        static_assert(std::is_floating_point_v<typename RField::value_type>);
        typename RField::value_type scaling;
    public :
        Simple (typename RField::value_type scaling_) : scaling { scaling_ }
        {
            #ifndef NDEBUG
            std::fprintf(stderr, "Initialized CallbackUtils::radius::Simple with scaling=%f.\n", scaling);
            #endif
        }
        Simple () : scaling { (typename RField::value_type)1.0 } { }
        float grp_radius (void **grp_properties) const override
        { return scaling * Callback<AFields>::template get_property<RField>(grp_properties); }
    };
}// namespace radius }}}

// Some common cases for how we may want to treat the data
namespace action
{// {{{
    // An abstract base class collecting different actions to be taken when a new group
    // is encountered.
    // The user should not inherit from this class directly
    // but rather through MultiGrpAction below.
    template<typename AFields>
    class MultiGrpActionBase : virtual public Callback<AFields>
    {
    protected :
        std::vector<std::pair<void *,std::function<void(void *, void **)>>> grp_actions;
    public :
        void grp_action (void **grp_properties) override final
        { for (auto fct : grp_actions) fct.second(fct.first, grp_properties); }
    };

    // User should inherit from this class to register an action with the MultiGrpActionBase
    // The constructor arguments are a pointer to the child class instance
    // and a pointer to a static member function of the child class that
    // takes a child class instance as its first argument and grp_properties as its second.
    template<typename AFields>
    struct MultiGrpAction : virtual public Callback<AFields>,
                            virtual private MultiGrpActionBase<AFields>
    {
        MultiGrpAction (void *obj, void (*fct)(void *, void**))
        { MultiGrpActionBase<AFields>::grp_actions.push_back(std::make_pair(obj, fct)); }
    };

    // Abstract base class that implements storage of homogeneous data for each group.
    // (data is particle dependent)
    // Each group gets an entry in the data vector, of type Tdata.
    // For each particle that should be considered part of a group,
    // the corresponding entry in the data vector is modified.
    // To this end, the user needs to override two functions :
    //      -- prt_reduce -> how to reduce a particle's properties to an item of a given type
    //      -- prt_combine -> how to combine the outcome of the previous function with
    //                        the existing entry in the data vector
    template<typename AFields, typename Tdata, typename Tfromprt=Tdata>
    class StorePrtHomogeneous : virtual public Callback<AFields>,
                                public MultiGrpAction<AFields>
    {
        std::vector<Tdata> *data;
        static void enlarge_data (void *obj, void **grp_properties)
        {
            StorePrtHomogeneous *p = (StorePrtHomogeneous *)obj;
            p->data->resize(p->data->size() + 1UL);
        }
    protected :
        virtual Tfromprt prt_reduce (size_t grp_idx, void **grp_properties, void **prt_properties, float R) const = 0;
        virtual void prt_combine (size_t grp_idx, Tdata &data_item, Tfromprt prt_item) const = 0;
    public :
        StorePrtHomogeneous (std::vector<Tdata> *data_) :
            data(data_),
            MultiGrpAction<AFields> { this, enlarge_data }
        { }
        void prt_action (size_t grp_idx, void **grp_properties, void **prt_properties, float R) override final
        {
            auto prt_item = prt_reduce(grp_idx, grp_properties, prt_properties, R);
            prt_combine(grp_idx, (*data)[grp_idx], prt_item);
        }
    };

    // Abstract base class that implements storage of homogeneous data for each group.
    // (data is group dependent)
    // Since each group is only read once, in contrast to StorePrtHomogeneous only
    // the grp_reduce function is required to be overriden by the subclass.
    template<typename AFields, typename Tdata>
    class StoreGrpHomogeneous : virtual public Callback<AFields>,
                                public MultiGrpAction<AFields>
    {
        std::vector<Tdata> *data;
        static void append_data (void *obj, void **grp_properties)
        {
            StoreGrpHomogeneous *p = (StoreGrpHomogeneous *)obj;
            p->data->push_back(p->grp_reduce(grp_properties));
        }
    protected :
        virtual Tdata grp_reduce (void **grp_properties) const = 0;
    public :
        StoreGrpHomogeneous (std::vector<Tdata> *data_) :
            data(data_),
            MultiGrpAction<AFields> { this, append_data }
        { }
    };

    // Implements storage of a single group property named Field
    // in a data vector.
    // Multiple template instantiations can be used to store multiple
    // group properties.
    template<typename AFields, typename Field, typename storeasT = typename Field::value_type>
    class StoreGrpProperty : virtual public Callback<AFields>,
                             public StoreGrpHomogeneous<AFields, storeasT>
    {
        static_assert(Field::dim == 1, "Currently not implemented, could probably be done");
        storeasT grp_reduce (void **grp_properties) const override final
        {
            return Callback<AFields>::template get_property<Field>(grp_properties);
        }
    public :
        StoreGrpProperty (std::vector<storeasT> *data) :
            StoreGrpHomogeneous<AFields, storeasT>(data)
        { }
    };
}// namespace action }}}

} // namespace CallbackUtils

#endif // CALLBACK_UTILS_HPP
