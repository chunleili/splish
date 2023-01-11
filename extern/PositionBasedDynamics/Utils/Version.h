#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "a102b34fe7054920f510b2980b8df13a37f2d012"
#define GIT_REFSPEC "refs/heads/nonNewton-decouple"
#define GIT_LOCAL_STATUS "DIRTY"

#ifdef DL_OUTPUT
#pragma WARNING(Local changes not committed.)
#endif

#endif
