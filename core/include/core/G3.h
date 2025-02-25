#ifndef _G3_H
#define _G3_H

#include <cereal/cereal.hpp>
#include <cereal/types/polymorphic.hpp>

#define G3_POINTERS(x) \
typedef std::shared_ptr<x> x##Ptr; \
typedef std::shared_ptr<const x> x##ConstPtr; \

#define G3_SERIALIZABLE(x, version) \
CEREAL_CLASS_VERSION(x, version);  \
CEREAL_REGISTER_TYPE_WITH_NAME(x, #x);

#endif
