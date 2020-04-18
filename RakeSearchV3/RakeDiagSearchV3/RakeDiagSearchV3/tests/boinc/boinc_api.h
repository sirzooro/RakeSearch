#ifndef BOINC_API_H
#define BOINC_API_H

#include <string>

inline void boinc_checkpoint_completed() {}

inline void boinc_fraction_done(double /*fraction*/) {}

inline int boinc_time_to_checkpoint() { return 0; }

inline void boinc_init() {}

inline void boinc_set_min_checkpoint_period(int) {}

inline int boinc_resolve_filename_s(const char* s1, std::string& s2)
{
    s2 = s1;
    return 0;
}

inline void boinc_finish(int n) { exit(n); }

#endif // BOINC_API_H
