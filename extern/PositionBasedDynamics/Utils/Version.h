#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "0340ab17c69a6ce708c488d6c27c292168f30c3a"
#define GIT_REFSPEC "refs/heads/chunlei-dev"
#define GIT_LOCAL_STATUS "CLEAN"

#ifdef DL_OUTPUT

#endif

#endif
