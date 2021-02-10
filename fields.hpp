#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <initializer_list>

#ifndef NDEBUG
#   include <cstdio>
#endif // NDEBUG

enum class FieldTypes { GrpFld, PrtFld };

template<FieldTypes field_type_, typename... Fields>
struct FieldCollection
{//{{{
    static constexpr const size_t Nfields   = sizeof...(Fields);
    static constexpr const char  *names[]   = { Fields::name ... };
    static constexpr const size_t strides[] = { Fields::stride ... };

    // store this information so we can use it to check order
    static constexpr const FieldTypes field_type = field_type_;

    // check that all fields belong here, i.e. are Particle or Group fields
    static constexpr bool all_valid (std::initializer_list<FieldTypes> types)
    {
        for (auto t : types)
            if (t != field_type_) return false;
        return true;
    }
    static_assert( all_valid({Fields::type...}),
                   "There is a Particle field in a Group type or vice versa." );

    // check that the first field is indeed a coordinate field
    template<typename first_field, typename... other_fields>
    static constexpr bool first_field_is_coord ()
    {
        return first_field::coord;
    }
    static_assert( first_field_is_coord<Fields...>(),
                   "The first field in one type is not a coordinate field.");

    // check that there is no duplication (which is not a problem per se but likely indicates a bug)
    static constexpr bool str_equ (const char *str1, const char *str2)
    {
        // this function simply does the same as C strcmp, but as a constexpr
        return *str1==*str2 && ( *str1=='\0' || str_equ(str1+1, str2+1) );
    }
    static constexpr bool all_unequal (std::initializer_list<const char *> str)
    {
        for (auto str1=str.begin(); str1 != str.end(); ++str1)
            for (auto str2=str1+1; str2 != str.end(); ++str2)
                if (str_equ(*str1, *str2))
                    return false;
        return true;
    }
    static_assert( all_unequal({ Fields::name ... }),
                   "Duplicate field, this is likely not what you intended to do.");
};//}}}

template<typename... Fields>
struct GrpFields : FieldCollection<FieldTypes::GrpFld, Fields...>
{ };

template<typename... Fields>
struct PrtFields : FieldCollection<FieldTypes::PrtFld, Fields...>
{ };

#ifndef NDEBUG
template<typename Fields>
void print_templ_arg_info_fct (const char *FieldsName)
{
    std::fprintf(stderr, "In the FieldsCollection %s are contained :\n", FieldsName);
    for (size_t ii=0; ii != Fields::Nfields; ++ii)
        std::fprintf(stderr, "\t[%2lu] %-20s   stride : %2lu byte\n",
                             ii, Fields::names[ii], Fields::strides[ii]);
}
#define print_templ_arg_info(name) print_templ_arg_info_fct<name>(#name)
#endif // NDEBUG

// A compile-time helper to see whether the user is correctly instatiating the templates
template<typename GroupFields, typename ParticleFields>
void template_checks ()
{
    static_assert(GroupFields::field_type == FieldTypes::GrpFld,
                  "First template parameter is not a GroupFields type");
    static_assert(ParticleFields::field_type == FieldTypes::PrtFld,
                  "First template parameter is not a ParticleFields type");

    #ifndef NDEBUG
    print_templ_arg_info(GroupFields);
    print_templ_arg_info(ParticleFields);
    #endif // NDEBUG
}

#endif // FIELDS_HPP
