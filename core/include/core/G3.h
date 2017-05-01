#ifndef _G3_H
#define _G3_H

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/polymorphic.hpp>

#include <G3Units.h>

#define G3_POINTERS(x) \
typedef boost::shared_ptr<x> x##Ptr; \
typedef boost::shared_ptr<const x> x##ConstPtr; \

#define G3_SERIALIZABLE(x, version) \
CEREAL_CLASS_VERSION(x, version);  \
CEREAL_REGISTER_TYPE_WITH_NAME(x, #x);

#endif
