#include <sys/resource.h>
#include <sys/time.h>
// #include <sys/types.h>
// #include <sys/sysinfo.h>
#include "yak-priv.h"

int yak_verbose = 3;

static double yak_realtime0;

double yak_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}


// size_t getCurrentRSS() {
// #ifdef _WIN32
//     PROCESS_MEMORY_COUNTERS pmc;
//     if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
//         return (size_t)(pmc.WorkingSetSize / 1024);
//     }
//     return 0;

// #elif __linux__
//     size_t rss = 0;
//     std::ifstream stat_stream("/proc/self/status", std::ios_base::in);
//     std::string line;

//     while (std::getline(stat_stream, line)) {
//         if (line.compare(0, 6, "VmRSS:") == 0) {
//             // 找到 VmRSS 行，格式通常为: "VmRSS:     1234 KB"
//             size_t first_digit = line.find_first_of("0123456789");
//             size_t last_digit = line.find_last_of("0123456789");
//             rss = std::stoul(line.substr(first_digit, last_digit - first_digit + 1));
//             break;
//         }
//     }
//     return rss;
// #else
//     return 0; // 不支持的平台
// #endif
// }

static inline double yak_realtime_core(void)
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void yak_reset_realtime(void)
{
	yak_realtime0 = yak_realtime_core();
}

double yak_realtime(void)
{
	return yak_realtime_core() - yak_realtime0;
}

long yak_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}
