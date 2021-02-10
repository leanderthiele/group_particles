#ifndef ILLUSTRIS_FIELDS_HPP
#define ILLUSTRIS_FIELDS_HPP

#include <type_traits>

#include "fields.hpp"

#define FIELD(name_, dim_, value_type_, type_, coord_)                \
    struct name_                                                      \
    {                                                                 \
        static constexpr const char name[] = #name_;                  \
        static constexpr const size_t stride = dim_                   \
                                               * sizeof(value_type_); \
        static constexpr const FieldTypes type = type_;               \
        static constexpr const bool coord = coord_;                   \
        static_assert(!coord_ || dim_==3,                             \
                      "Non-3dimensional coordinate field "#name_);    \
        static_assert(!coord_ || std::is_same_v<float,value_type_>,   \
                      "Non-float coordinate field "#name_             \
                      " not supported");                              \
    }

// Particle fields
FIELD(CenterOfMass, 3, float, FieldTypes::PrtFld, true);
FIELD(Coordinates, 3, float, FieldTypes::PrtFld, true);
FIELD(Density, 1, float, FieldTypes::PrtFld, false);
FIELD(ElectronAbundance, 1, float, FieldTypes::PrtFld, false);
// TODO

// Group fields
FIELD(GroupCM, 3, float, FieldTypes::GrpFld, true);
FIELD(GroupPos, 3, float, FieldTypes::GrpFld, true);
FIELD(GroupMass, 1, float, FieldTypes::GrpFld, false);
FIELD(GroupVel, 3, float, FieldTypes::GrpFld, false);
// TODO



#endif // ILLUSTRIS_FIELDS_HPP
