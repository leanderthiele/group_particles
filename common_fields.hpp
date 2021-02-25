#ifndef ILLUSTRIS_FIELDS_HPP
#define ILLUSTRIS_FIELDS_HPP

#include <type_traits>

#include "fields.hpp"

namespace IllustrisFields
{
    // Particle fields
    FIELD(CenterOfMass, 3, float, FieldTypes::PrtFld, true);
    FIELD(Coordinates, 3, double, FieldTypes::PrtFld, true);
    FIELD(Density, 1, float, FieldTypes::PrtFld, false);
    FIELD(ElectronAbundance, 1, float, FieldTypes::PrtFld, false);
    FIELD(Masses, 1, float, FieldTypes::PrtFld, false);
    FIELD(InternalEnergy, 1, float, FieldTypes::PrtFld, false);
    // TODO

    // Group fields
    FIELD(GroupCM, 3, float, FieldTypes::GrpFld, true);
    FIELD(GroupPos, 3, float, FieldTypes::GrpFld, true);
    FIELD(GroupMass, 1, float, FieldTypes::GrpFld, false);
    FIELD(GroupVel, 3, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Crit200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Crit200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Crit500, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Crit500, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Mean200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Mean200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_TopHat200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_TopHat200, 1, float, FieldTypes::GrpFld, false);
    // TODO
} // namespace IllustrisFields

namespace GadgetFields
{
    // Particle fields
    FIELD(Coordinates, 3, float, FieldTypes::PrtFld, true);
    FIELD(Velocities, 3, float, FieldTypes::PrtFld, false);

    // Group fields
    FIELD(GroupPos, 3, float, FieldTypes::GrpFld, true);
    FIELD(GroupVel, 3, float, FieldTypes::GrpFld, false);
    FIELD(GroupMass, 1, float, FieldTypes::GrpFld, false);
    FIELD(GroupAscale, 1, float, FieldTypes::GrpFld, false);
} // namespace GadgetFields



#endif // ILLUSTRIS_FIELDS_HPP
