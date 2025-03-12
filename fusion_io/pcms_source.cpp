#include "pcms_source.h"

pcms_source::pcms_source(const int type, const char* filename, PCMS_Library& lib) {
  this->lib = &lib;
  int result = fio_open_source(&src, type, filename);

  if(result != FIO_SUCCESS) {
    std::cerr << "Error opening file" << std::endl;
    exit(1);
  };
  src->lib = &lib;
  src->add_pcms_fields();
}

pcms_source::~pcms_source() {
  // fio_close_source(&src);
}