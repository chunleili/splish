#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "7c1e5e3c6f590a7679f314a1f2b161521fee5401"
#define GIT_REFSPEC "refs/heads/bugfix"
#define GIT_LOCAL_STATUS "DIRTY"

#ifdef DL_OUTPUT
#pragma WARNING(Local changes not committed.)
#endif

#endif
