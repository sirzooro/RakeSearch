#ifndef BOINC_API_H
#define BOINC_API_H

inline void boinc_checkpoint_completed() {}

inline void boinc_fraction_done(double /*fraction*/) {}

inline int boinc_time_to_checkpoint() { return 0; }

#endif // BOINC_API_H
