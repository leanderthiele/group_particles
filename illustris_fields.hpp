#ifndef ILLUSTRIS_FIELDS_HPP
#define ILLUSTRIS_FIELDS_HPP

#include <type_traits>

#include "fields.hpp"

namespace IllustrisFields
{

// Particle fields
// -- Note that the Coordinates field is double for some reason,
//    at the moment the code can't deal with this
//    The problem is that while we have CenterOfMass for gas particles,
//    only Coordinates (which are double) are available for DM particles.
//    FIXME
//    Probably need to use if constexpr.
//    The only problem really is how to do calculations with the numbers,
//    but if constexpr should do the trick
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
// TODO


} // IllustrisFields



#endif // ILLUSTRIS_FIELDS_HPP
