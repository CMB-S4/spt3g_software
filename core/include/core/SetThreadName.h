#ifndef SETTHREADNAME_H
#define SETTHREADNAME_H

#include <string>

#include <pthread.h>
#ifdef __FreeBSD__
#include <pthread_np.h>
#endif

/// Set the name for the currently-running thread
/// \param name the thread name to set, which should not exceed 15 characters
inline void setThreadName(std::string name){
	if(name.size()>15) // Linux limits thread names to 15 bytes
		name=name.substr(0,15);
#if defined(__linux__) || defined(__FreeBSD__) \
		|| defined(__NetBSD__) || defined(__OpenBSD__)
	pthread_setname_np(pthread_self(), name.c_str());
#endif
#ifdef __APPLE__
	// Apple refuses to allow setting the name on anything but the current
	// thread, for some reason
	pthread_setname_np(name.c_str());
#endif
}

#endif //SETTHREADNAME_H
