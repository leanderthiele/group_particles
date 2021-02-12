#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <initializer_list>
#include <type_traits>

#ifndef NDEBUG
#   include <cstdio>
#endif // NDEBUG

enum class FieldTypes { GrpFld, PrtFld };

#define FIELD(name_, dim_, value_type_, type_, coord_)                \
    struct name_                                                      \
    {                                                                 \
        name_ () = delete;                                            \
        static constexpr const char name[] = #name_;                  \
        static constexpr const size_t size = sizeof(value_type_);     \
        static constexpr const size_t dim  = dim_;                    \
        static constexpr const size_t stride = dim_                   \
                                               * sizeof(value_type_); \
        static constexpr const FieldTypes type = type_;               \
        static constexpr const bool coord = coord_;                   \
        typedef value_type_ value_type;                               \
        static_assert(!coord_ || dim_==3,                             \
                      "Non-3dimensional coordinate field "#name_);    \
        static_assert(!coord_                                         \
                      || std::is_same_v<float,value_type_>            \
                      || std::is_same_v<double,value_type_>,          \
                      "Non-float coordinate field "#name_             \
                      " not supported");                              \
    }

template<FieldTypes field_type_, typename... Fields>
class FieldCollection
{//{{{
    static constexpr bool all_valid (std::initializer_list<FieldTypes> types)
    {// {{{
        for (auto t : types)
            if (t != field_type_) return false;
        return true;
    }// }}}

    template<typename first_field, typename... other_fields>
    static constexpr bool first_field_is_coord ()
    {// {{{
        return first_field::coord;
    }// }}}

    static constexpr bool str_equ (const char *str1, const char *str2)
    {// {{{
        // this function simply does the same as C strcmp, but as a constexpr
        return *str1==*str2 && ( *str1=='\0' || str_equ(str1+1, str2+1) );
    }// }}}

    template<typename first_field, typename second_field, typename... other_fields>
    static constexpr bool all_unequal ()
    {// {{{
        constexpr bool first_two = !std::is_same_v<first_field,second_field>;
        if constexpr (sizeof...(other_fields))
            return first_two
                   && all_unequal<first_field, other_fields...>()
                   && all_unequal<second_field, other_fields...>();
        else
            return first_two;
    }// }}}

    template<typename to_find, typename first_field, typename... other_fields>
    static constexpr size_t find_idx (size_t idx_of_first = 0UL)
    {// {{{
        if constexpr (std::is_same_v<to_find, first_field>)
            return idx_of_first;
        else
        {
            static_assert( sizeof...(other_fields),
                           "Requested index of field that is not contained in collection." );
            return find_idx<to_find, other_fields...>(idx_of_first+1UL);
        }
    }// }}}

    template<typename first_field, typename... other_fields>
    struct extract_coord_type
    {// {{{
        typedef typename first_field::value_type value_type;
    };// }}}
public :
    // this is only a type, no instances of it can be created
    FieldCollection () = delete;

    static constexpr const size_t Nfields   = sizeof...(Fields);
    static constexpr const char  *names[]   = { Fields::name ... };
    static constexpr const size_t sizes[]   = { Fields::size ... };
    static constexpr const size_t dims[]    = { Fields::dim ... };
    static constexpr const size_t strides[] = { Fields::stride ... };

    // we need to do arithmetic with the coordinate values, so we need to know their type
    typedef typename extract_coord_type<Fields...>::value_type coord_t;

    // store this information so we can use it to check order
    static constexpr const FieldTypes field_type = field_type_;

    // helper for the user so they don't have to remember in which order
    // they passed the fields
    template<typename T>
    static constexpr size_t idx = find_idx<T, Fields...>();

    static_assert( Nfields,
                   "empty field collection not allowed, need at least the "
                   "coordinate field" );

    // check that all fields belong here, i.e. are Particle or Group fields
    static_assert( all_valid({Fields::type...}),
                   "There is a Particle field in a Group type or vice versa." );

    // check that the first field is indeed a coordinate field
    static_assert( first_field_is_coord<Fields...>(),
                   "The first field in one type is not a coordinate field.");

    // check that there is no duplication (which is not a problem per se but likely indicates a bug)
    static_assert( all_unequal<Fields...>(),
                   "Duplicate field, this is likely not what you intended to do.");
};//}}}

template<typename... Fields>
struct GrpFields : FieldCollection<FieldTypes::GrpFld, Fields...>
{ };

template<typename... Fields>
struct PrtFields : FieldCollection<FieldTypes::PrtFld, Fields...>
{ };

template<typename GroupFields_, typename ParticleFields_>
struct AllFields
{
    typedef GroupFields_    GroupFields;
    typedef ParticleFields_ ParticleFields;

    static_assert( GroupFields::field_type == FieldTypes::GrpFld,
                   "First template parameter for AllFields must be a GrpFields type" );
    static_assert( ParticleFields::field_type == FieldTypes::PrtFld,
                   "Second template parameter for AllFields must be a PrtFields type" );

    static void print_field_info ()
    {
        print_field_info_fct<GroupFields>("GroupFields");
        print_field_info_fct<ParticleFields>("ParticleFields");
    }

private :
    template<typename Fields>
    static void print_field_info_fct (const char *FieldsName)
    {
        std::fprintf(stderr, "In the FieldsCollection %s are contained :\n", FieldsName);
        for (size_t ii=0; ii != Fields::Nfields; ++ii)
            std::fprintf(stderr, "\t[%2lu] %-20s   stride : %2lu byte\n",
                                 ii, Fields::names[ii], Fields::strides[ii]);
    }

};

#endif // FIELDS_HPP
