#ifndef _G3_H
#define _G3_H

#include <cereal/cereal.hpp>
#include <cereal/types/polymorphic.hpp>

#define G3_POINTERS(x) \
typedef std::shared_ptr<x> x##Ptr; \
typedef std::shared_ptr<const x> x##ConstPtr; \

#define G3_SERIALIZABLE(x, version) \
CEREAL_CLASS_VERSION(x, version);  \
CEREAL_REGISTER_TYPE_WITH_NAME(x, #x); \
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(x, cereal::specialization::member_serialize);

#define G3_SPLIT_SERIALIZABLE(x, version) \
CEREAL_CLASS_VERSION(x, version);  \
CEREAL_REGISTER_TYPE_WITH_NAME(x, #x); \
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(x, cereal::specialization::member_load_save);

template <class T>
inline uint32_t
g3_class_version(T *)
{
	return cereal::detail::Version<T>::version;
}

#define G3_CHECK_VERSION(v) \
	if ((uint32_t)v > g3_class_version(this)) \
		log_fatal("Trying to read newer class version (%d) than " \
		    "supported (%d). Please upgrade your software.", v, \
		    g3_class_version(this));

#endif
