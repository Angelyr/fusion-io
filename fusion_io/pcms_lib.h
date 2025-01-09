#ifndef PCMS_LIB_H
#define PCMS_LIB_H

#include <pcms/server.h>
#include <pcms/client.h>
#include <pcms/field.h>
#include <pcms/types.h>

namespace {
  template <typename T>
  struct IndirectionLevel {
      static const int value = 0;
  };
  template <typename T>
  struct IndirectionLevel<T*> {
      static const int value = IndirectionLevel<T>::value + 1;
  };
}

namespace fusion_io {
  class Library {
    public:
    MPI_Comm comm;
    pcms::CouplerClient* client;

    Library(int argc, char** argv) {
      MPI_Init(&argc, &argv);
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);
      client = new pcms::CouplerClient("pcms_client", comm);
    }

    ~Library() {
      MPI_Finalize();
      delete client;
    }
  };

  template <typename fieldtype>
  class FusionFieldAdapter {
    public:
      using value_type = std::remove_pointer_t<std::remove_pointer_t<fieldtype>>;
      fieldtype* field;
      int size;
      std::vector<pcms::GO> gids_;

      FusionFieldAdapter(fieldtype* fieldIn, int sizeIn=0) {
        field = fieldIn;
        size = sizeIn;
        gids_ = std::vector<pcms::GO>(size);
        std::iota(gids_.begin(), gids_.end(), static_cast<pcms::GO>(0));
      }

      [[nodiscard]] std::vector<pcms::GO> GetGids() const { 
        return gids_;
      }

      [[nodiscard]] pcms::ReversePartitionMap GetReversePartitionMap(const redev::Partition& partition) const {
        return {};
      }

      template < typename = typename std::enable_if< IndirectionLevel<fieldtype>::value==1 > >
      int Serialize(pcms::ScalarArrayView<fieldtype, pcms::HostMemorySpace> buffer,
                    pcms::ScalarArrayView<const pcms::LO, pcms::HostMemorySpace> permutation) const
      {
        if (buffer.size() > 0) {
          for (pcms::LO i = 0; i < size; i++) {
            buffer[i] = field[permutation[i]];
          }
        }
        return size;
      }

      template < typename = typename std::enable_if< IndirectionLevel<fieldtype>::value==1 > >
      void Deserialize(pcms::ScalarArrayView<const fieldtype, pcms::HostMemorySpace> buffer,
                       pcms::ScalarArrayView<const pcms::LO, pcms::HostMemorySpace> permutation) const
      {
        REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
        for (size_t i = 0; i < buffer.size(); ++i) {
          field[permutation[i]] = buffer[i];
        }
      }
      
      template <typename T1, typename T2>
      int Serialize(T1, T2) const noexcept { return 0;}

      template <typename T1, typename T2>
      void Deserialize(T1, T2) const noexcept {}
  };
}

#endif //PCMS_LIB_H