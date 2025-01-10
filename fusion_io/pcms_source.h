#ifdef PCMS_ENABLED

#ifndef PCMS_SOURCE_H
#define PCMS_SOURCE_H

#include "fusion_io.h"
#include "pcms_lib.h"

inline int fio_open_source(fio_source** src, const int type, const char* filename, fusion_io::Library lib) {
  int result = fio_open_source(src, type, filename);
  (*src)->add_pcms_fields(lib);
  return result;
}

inline Omega_h::HostRead<double> evaluate(fio_field* field, double* x, double* values) {
  field->eval(x, values);
  Omega_h::HostWrite<double> host_write(field->dimension());
  for (int i=0; i<field->dimension(); i++) {
    host_write[i] = values[i];
  }
  Omega_h::HostRead<double> host_read(host_write.write());
  return host_read;
}

#endif //PCMS_SOURCE_H

#endif //PCMS_ENABLED