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
#include "hdf5_utils.hpp"

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
        std::string grp_name () const override
        { return "Group/"; }
        std::string prt_name () const override
        { return "PartType"+std::to_string(PartType)+"/"; }
    };

    template<typename AFields>
    struct Gadget : virtual public Callback<AFields>
    {
        std::string grp_name () const override
        { return "Group/"; }
        std::string prt_name () const override
        { return "PartType1/"; }
    };
} // namespace name }}}

// Some common things you'd want to store in terms of meta data
namespace meta_init
{// {{{
    template<typename AFields>
    class MultiPrtMetaInitBase : virtual public Callback<AFields>
    {
        static constexpr size_t buf_size = 64UL;
        size_t N_inits = 0UL;
        std::pair<void *, std::function<void(void *, std::shared_ptr<H5::H5File>)>> inits[buf_size];
    public :
        void read_prt_meta_init (std::shared_ptr<H5::H5File> fptr) override final
        { for (size_t ii=0; ii != N_inits; ++ii) inits[ii].second(inits[ii].first, fptr); }
        void register_init (void *obj, std::function<void(void *, std::shared_ptr<H5::H5File>)> fct)
        {
            inits[N_inits++] = std::make_pair(obj, fct);
            assert(N_inits < buf_size);
        }
    };

    template<typename AFields>
    struct MultiPrtMetaInit : virtual public Callback<AFields>,
                              virtual private MultiPrtMetaInitBase<AFields>
    {
        MultiPrtMetaInit (void *obj, std::function<void(void *, std::shared_ptr<H5::H5File>)> fct)
        { MultiPrtMetaInitBase<AFields>::register_init(obj, fct); }
    };

    template<typename AFields>
    class IllustrisCosmology : virtual public Callback<AFields>,
                               public MultiPrtMetaInit<AFields>
    {
        static void read_cosmology (void *obj, std::shared_ptr<H5::H5File> fptr)
        {
            IllustrisCosmology *p = (IllustrisCosmology *)obj;

            auto header = fptr->openGroup("/Header");
            #define READ(x) p->x = hdf5Utils::read_scalar_attr<double,double>(header, #x)
            READ(HubbleParam);
            READ(Omega0);
            READ(OmegaLambda);
            READ(OmegaBaryon);
            READ(Redshift);
            READ(Time);
            #undef READ
        }
    public :
        // data members the user may want to use
        double HubbleParam, Omega0, OmegaLambda, OmegaBaryon, Redshift, Time;

        IllustrisCosmology () :
            MultiPrtMetaInit<AFields> { this, read_cosmology }
        { }
    };

    template<typename AFields>
    class IllustrisMassTable : virtual public Callback<AFields>,
                               public MultiPrtMetaInit<AFields>
    {
        static void read_mass_table (void *obj, std::shared_ptr<H5::H5File> fptr)
        {
            IllustrisMassTable *p = (IllustrisMassTable *)obj;

            auto header = fptr->openGroup("/Header");
            p->Ntypes = hdf5Utils::read_vector_attr<double,double>(header, "MassTable", p->MassTable);
        }
    public :
        // data members the user may want to use
        // use a large enough buffer for all cases
        double MassTable[16];
        size_t Ntypes;

        IllustrisMassTable () :
            MultiPrtMetaInit<AFields> { this, read_mass_table }
        { }

    };
} // namespace meta_init }}}

// Some common ways to get meta data
namespace meta
{// {{{
    template<typename AFields, uint8_t PartType>
    class Illustris : virtual public Callback<AFields>
    {
    public :
        void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr, size_t &Ngroups) const override
        {
            auto header = fptr->openGroup("/Header");
            Ngroups = hdf5Utils::read_scalar_attr<int32_t,size_t>(header, "Ngroups_ThisFile");
        }
        void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr, coord_t &Bsize, size_t &Npart) const override
        {
            auto header = fptr->openGroup("/Header");
            Bsize = hdf5Utils::read_scalar_attr<double,coord_t>(header, "BoxSize");
            Npart = hdf5Utils::read_vector_attr<int32_t,size_t>(header, "NumPart_ThisFile", PartType);
        }
    };
}// namespace meta }}}

// Some common cases for the selection functions
namespace select
{// {{{
    template<typename AFields>
    struct GrpAll : virtual public Callback<AFields>
    {
        bool grp_select (const typename Callback<AFields>::GrpProperties &grp) const
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
        bool grp_select (const typename Callback<AFields>::GrpProperties &grp) const override
        {
            auto M = grp.template get<MassField>();
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
        coord_t grp_radius (const typename Callback<AFields>::GrpProperties &grp) const override
        { return scaling * grp.template get<RField>(); }
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
        // very very unlikely that any user would want to use this more than for 64 different actions
        // If that case happens, we have the assert below.
        static constexpr size_t buf_size = 64UL;
        size_t N_actions = 0UL;
        std::pair<void *, std::function<void(void *, const typename Callback<AFields>::GrpProperties &)>> grp_actions[buf_size];
    public :
        void grp_action (const typename Callback<AFields>::GrpProperties &grp) override final
        { for (size_t ii=0; ii != N_actions; ++ii) grp_actions[ii].second(grp_actions[ii].first, grp); }
        void register_fct (void *obj, std::function<void(void *, const typename Callback<AFields>::GrpProperties &)> fct)
        {
            grp_actions[N_actions++] = std::make_pair(obj, fct);
            assert(N_actions < buf_size);
        }
    };

    // User should inherit from this class to register an action with the MultiGrpActionBase
    // The constructor arguments are a pointer to the child class instance
    // and a pointer to a static member function of the child class that
    // takes a child class instance as its first argument and grp_properties as its second.
    template<typename AFields>
    struct MultiGrpAction : virtual public Callback<AFields>,
                            virtual private MultiGrpActionBase<AFields>
    {
        MultiGrpAction (void *obj, std::function<void(void *, const typename Callback<AFields>::GrpProperties &)> fct)
        { MultiGrpActionBase<AFields>::register_fct(obj, fct); }
    };

    // Abstract base class that implements storage of homogeneous data for each group.
    // (data is particle dependent)
    //
    // Each group gets an entry in the data vector, of type Tdata.
    // If Tdata has a constructor of the signature
    //      Tdata (const typename Callback<AFields::GrpProperties &),
    // this constructor will be called. In this case, the default constructor should be
    // deleted (this ensures that we have an extra check that the custom constructor
    // was declared correctly).
    // Otherwise, the default constructor will be called.
    //
    // For each particle that should be considered part of a group,
    // the corresponding entry in the data vector is modified.
    // To this end, the user needs to override one function :
    //      prt_insert -> specifies how a particle is supposed to be combined with the
    //                    existing element of type Tdata in the data vector
    template<typename AFields, typename Tdata>
    class StorePrtHomogeneous : virtual public Callback<AFields>,
                                public MultiGrpAction<AFields>
    {
        std::vector<Tdata> &data;
        static void enlarge_data (void *obj, const typename Callback<AFields>::GrpProperties &grp)
        {
            StorePrtHomogeneous *p = (StorePrtHomogeneous *)obj;

            if constexpr (std::is_constructible_v<Tdata, const typename Callback<AFields>::GrpProperties &>)
            {
                static_assert(!std::is_default_constructible_v<Tdata>);
                p->data.emplace_back(grp);
            }
            else
                p->data.emplace_back();
        }
    protected :
        virtual void prt_insert (size_t grp_idx,
                                 const typename Callback<AFields>::GrpProperties &grp,
                                 const typename Callback<AFields>::PrtProperties &prt,
                                 coord_t Rsq,
                                 Tdata &data_item) = 0;
    public :
        StorePrtHomogeneous (std::vector<Tdata> &data_) :
            data(data_),
            MultiGrpAction<AFields> { this, enlarge_data }
        { }
        void prt_action (size_t grp_idx,
                         const typename Callback<AFields>::GrpProperties &grp,
                         const typename Callback<AFields>::PrtProperties &prt,
                         coord_t Rsq) override final
        {
            prt_insert(grp_idx, grp, prt, Rsq, data[grp_idx]);
        }
    };

    // Abstract base class that implements storage of homogeneous data for each group.
    // (data is group dependent)
    template<typename AFields, typename Tdata>
    class StoreGrpHomogeneous : virtual public Callback<AFields>,
                                public MultiGrpAction<AFields>
    {
        std::vector<Tdata> &data;
        static void append_data (void *obj, const typename Callback<AFields>::GrpProperties &grp)
        {
            StoreGrpHomogeneous *p = (StoreGrpHomogeneous *)obj;
            p->data.push_back(p->grp_reduce(grp));
        }
    protected :
        virtual Tdata grp_reduce (const typename Callback<AFields>::GrpProperties &grp) const = 0;
    public :
        StoreGrpHomogeneous (std::vector<Tdata> &data_) :
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
    private :
        storeasT grp_reduce (const typename Callback<AFields>::GrpProperties &grp) const override final
        {
            return grp.template get<Field>();
        }
    public :
        StoreGrpProperty (std::vector<storeasT> &data) :
            StoreGrpHomogeneous<AFields, storeasT>(data)
        { }
    };
}// namespace action }}}

} // namespace CallbackUtils

#endif // CALLBACK_UTILS_HPP
