
// -----------------  LOGGING  -----------------------

//#define DEBUG_PRINT_ENABLED

#define INFO(fmt, args...) \
    do { \
        printf("INFO: " fmt, ## args); \
    } while (0)
#define WARN(fmt, args...) \
    do { \
        printf("WARN: " fmt, ## args); \
    } while (0)
#define ERROR(fmt, args...) \
    do { \
        printf("ERROR: " fmt, ## args); \
    } while (0)
#define FATAL(fmt, args...) \
    do { \
        printf("FATAL: " fmt, ## args); \
        exit(1); \
    } while (0)

#ifdef DEBUG_PRINT_ENABLED
#define DEBUG(fmt, args...) \
    do { \
        printf("DEBUG: " fmt, ## args); \
    } while (0)
#else
#define DEBUG(fmt, args...)
#endif
