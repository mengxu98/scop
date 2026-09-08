#ifndef SCOP_DYNLOAD_H
#define SCOP_DYNLOAD_H

// Process-local symbol lookup for optional BLAS entry points. Windows CI
// does not provide <dlfcn.h>; LoadLibrary / GetProcAddress is the equivalent.

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

inline void* scop_dlopen(const char* path) {
  if (path == nullptr) {
    return static_cast<void*>(GetModuleHandleA(nullptr));
  }
  HMODULE existing = GetModuleHandleA(path);
  if (existing != nullptr) {
    return static_cast<void*>(existing);
  }
  return static_cast<void*>(LoadLibraryA(path));
}

inline void* scop_dlsym(void* handle, const char* name) {
  if (handle == nullptr || name == nullptr) {
    return nullptr;
  }
  // FARPROC is a function pointer; MSVC/MinGW reject static_cast to void*.
  return reinterpret_cast<void*>(GetProcAddress(static_cast<HMODULE>(handle), name));
}

#else
#include <dlfcn.h>

inline void* scop_dlopen(const char* path) {
  return dlopen(path, RTLD_LAZY | RTLD_LOCAL);
}

inline void* scop_dlsym(void* handle, const char* name) {
  return dlsym(handle, name);
}
#endif

#endif
