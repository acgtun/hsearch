#ifndef SDK_H_
#define SDK_H_

#include <time.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define MAX_LINE_LEN 6000

inline void MemoryAllocateCheck(void* pointer, const char* file, int line) {
  if (pointer == NULL) {
    printf("Memory allocate error in %s at line %d\n", file, line);
    exit (EXIT_FAILURE);
  }
}

inline void FileOpenCheck(FILE* pfile, const char* file, int line) {
  if (pfile == NULL) {
    printf("File open error in %s at line %d\n", file, line);
    exit (EXIT_FAILURE);
  }
}

#define HANDLE_ERROR(err) (HandleError( err, __FILE__, __LINE__ ))
#define LOG_INFO printf("--- %s:%s:%d\n",  __FILE__, __func__, __LINE__)
#define FILE_OPEN_CHECK(pfile) (FileOpenCheck( pfile, __FILE__, __LINE__))
#define MEMORY_ALLOCATE_CHECK(pointer)  (MemoryAllocateCheck(pointer, __FILE__, __LINE__))

#define FREAD_CHECK(func, size) { \
	uint32_t s = func; \
	if(s != size) { \
		printf("read file error. --- %s:%s:%d\n", __FILE__, __func__, __LINE__); \
		exit(EXIT_FAILURE); \
	} \
}

#define ERROR_INFO(msg) { \
	printf("--ERROR INFO-- %s --- %s:%s:%d\n", msg, __FILE__, __func__, __LINE__); \
	exit(EXIT_FAILURE); \
}

#define TIME_INFO(func, msg) { \
	clock_t start_t, end_t; \
	start_t = clock(); \
	func; \
	end_t = clock(); \
	printf("--INFO-- %s takes %.3lf seconds.\n", msg, (double) ((end_t - start_t) / CLOCKS_PER_SEC )); \
}

inline void INFO(const char* msg) {
  printf("--INFO-- %s\n", msg);
}
inline void INFO(const char* msg, const char* val) {
  printf("--INFO-- %s %s\n", msg, val);
}
inline void INFO(const char* msg, const uint64_t &val) {
  printf("--INFO-- %s %" PRIu64 "\n", msg, val);
}
inline void INFO(const char* msg, const uint32_t &val) {
  printf("--INFO-- %s" "%" PRIu32 "\n", msg, val);
}
inline void INFO(const char* msg, const int &val) {
  printf("--INFO-- %s %d\n", msg, val);
}
inline void INFO(const char* msg, const double &val) {
  printf("--INFO-- %s %lf\n", msg, val);
}
inline void INFO(const char* msg, const uint32_t & val1, const char* val2) {
  printf("--INFO-- %s" "%" PRIu32 " %s\n", msg, val1, val2);
}
inline void INFO(const char* msg, const uint64_t & val1, const char* val2) {
  printf("--INFO-- %s" "%" PRIu64 " %s\n", msg, val1, val2);
}

#endif /* SDK_H_ */
