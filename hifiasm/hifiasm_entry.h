#ifndef HIFIASM_ENTRY_H
#define HIFIASM_ENTRY_H

#if defined(__GNUC__)
#define HIFIASM_API __attribute__((visibility("default")))
#else
#define HIFIASM_API
#endif

HIFIASM_API int hifiasm_main(int argc, char *argv[]);

#endif
