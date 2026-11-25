#ifndef SCOP_DYNLOAD_H
#define SCOP_DYNLOAD_H


#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

inline void* open_dynamic_library(const char* path) {
  if (path == nullptr) {
    return static_cast<void*>(GetModuleHandleA(nullptr));
  }
  HMODULE existing = GetModuleHandleA(path);
  if (existing != nullptr) {
    return static_cast<void*>(existing);
  }
  return static_cast<void*>(LoadLibraryA(path));
}

inline void* find_dynamic_symbol(void* handle, const char* name) {
  if (handle == nullptr || name == nullptr) {
    return nullptr;
  }
  return reinterpret_cast<void*>(GetProcAddress(static_cast<HMODULE>(handle), name));
}

#else
#include <dlfcn.h>

inline void* open_dynamic_library(const char* path) {
  return dlopen(path, RTLD_LAZY | RTLD_LOCAL);
}

inline void* find_dynamic_symbol(void* handle, const char* name) {
  return dlsym(handle, name);
}
#endif

#endif
