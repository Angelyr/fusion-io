#ifdef PCMS_ENABLED

#ifndef PCMS_SOURCE_H
#define PCMS_SOURCE_H

#include "fusion_io.h"
#include "pcms_lib.h"

class pcms_source {
  public:

  fio_source* src;
  PCMS_Library* lib;

  pcms_source() {}
  pcms_source(const int type, const char* filename, PCMS_Library& lib);
  ~pcms_source();
};

#endif //PCMS_SOURCE_H

#endif //PCMS_ENABLED