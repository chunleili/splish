#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "34fc5aa2f654db23fe5bad8bbe311b12118003ef"
#define GIT_REFSPEC "refs/heads/nonNewton-decouple"
#define GIT_LOCAL_STATUS "CLEAN"

#ifdef DL_OUTPUT

#endif

#endif
