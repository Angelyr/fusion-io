#ifdef PCMS_ENABLED

#ifndef PCMS_LIB_H
#define PCMS_LIB_H

#include <pcms/server.h>
#include <pcms/client.h>
#include <pcms/field.h>
#include <pcms/types.h>

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

  struct GetRank
  {
    GetRank(const int& id, const int& dim) : id_(id), dim_(dim) {}
    auto operator()(const redev::ClassPtn& ptn) const
    {
      PCMS_FUNCTION_TIMER;
      const auto ent = redev::ClassPtn::ModelEnt({dim_, id_});
      return ptn.GetRank(ent);
    }
    auto operator()(const redev::RCBPtn& /*unused*/) const
    {
      PCMS_FUNCTION_TIMER;
      std::cerr << "RCB partition not handled yet\n";
      std::terminate();
      return 0;
    }
    const int& id_;
    const int& dim_;
  };

  template <typename fieldType>
  class FieldAdapter {
    public:
      using value_type = std::remove_pointer_t<std::remove_pointer_t<std::remove_pointer_t<fieldType>>>;
      using memory_space = pcms::HostMemorySpace;

      FieldAdapter(fieldType fieldIn, int sizeX=1, int sizeY=1, int sizeZ=1) :
        field(fieldIn), 
        size{sizeX, sizeY, sizeZ} {
        int totalSize = sizeX*sizeY*sizeZ;
        gids = std::vector<pcms::GO>(totalSize);
        std::iota(gids.begin(), gids.end(), static_cast<pcms::GO>(0));
      }

      [[nodiscard]] std::vector<pcms::GO> GetGids() const { 
        return gids;
      }

      [[nodiscard]] pcms::ReversePartitionMap GetReversePartitionMap(const redev::Partition& partition) const {
        pcms::ReversePartitionMap reverse_partition;
        int local_index = 0;
        for (int i=0; i < totalSize; i++) {
          auto dr = std::visit(GetRank{i, 0}, partition);
          reverse_partition[dr].emplace_back(local_index++);
        }
        return reverse_partition;
      }

      template<typename T>
      using remove1Pointer = typename std::remove_pointer<T>::type;
      template<typename T>
      using remove2Pointers = remove1Pointer<remove1Pointer<T>>;
      template<typename T>
      using remove3Pointers = remove1Pointer<remove1Pointer<remove1Pointer<T>>>;

      //Serialize 1D Pointers
      template <typename U = fieldType>
      typename std::enable_if<std::is_pointer<U>::value 
                          && !std::is_pointer<remove1Pointer<U>>::value, int>::type
      Serialize(pcms::ScalarArrayView<value_type, memory_space> buffer,
                    pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
      {
        REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
        if (buffer.size() > 0) {
          for (int i = 0; i < totalSize; i++) {
            buffer[i] = field[permutation[i]];
          }
        }
        return totalSize;
      }

      //Deserialize 1D Pointers
      template <typename U = fieldType>
      typename std::enable_if<std::is_pointer<U>::value 
                          && !std::is_pointer<remove1Pointer<U>>::value, void>::type
      Deserialize(pcms::ScalarArrayView<const value_type, memory_space> buffer,
                       pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
      {
        REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
        for (int i = 0; i < buffer.size(); ++i) {
          field[permutation[i]] = buffer[i];
        }
      }

      //Serialize 2D Pointers
      template <typename U = fieldType>
      typename std::enable_if<std::is_pointer<remove1Pointer<U>>::value 
                          && !std::is_pointer<remove2Pointers<U>>::value, int>::type
      Serialize(pcms::ScalarArrayView<value_type, memory_space> buffer,
                    pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
      {
        REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
        if (buffer.size() > 0) {
          for (int i = 0; i < totalSize; i++) {
            int x = permutation[i] % size[0];
            int y = permutation[i] / size[0];
            buffer[i] = field[x][y];
          }
        }
        return totalSize;
      }

      //Deserialize 2D Pointers
      template <typename U = fieldType>
      typename std::enable_if<std::is_pointer<remove1Pointer<U>>::value 
                          && !std::is_pointer<remove2Pointers<U>>::value, void>::type
      Deserialize(pcms::ScalarArrayView<const value_type, memory_space> buffer,
                       pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
      {
        REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
        for (int i = 0; i < buffer.size(); ++i) {
          int x = permutation[i] % size[0];
          int y = permutation[i] / size[0];
          field[x][y] = buffer[i];
        }
      }

      //Serialize 3D Pointers
      template <typename U = fieldType>
      typename std::enable_if<std::is_pointer<remove2Pointers<U>>::value 
                          && !std::is_pointer<remove3Pointers<U>>::value, int>::type
      Serialize(pcms::ScalarArrayView<value_type, memory_space> buffer,
                    pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
      {
        REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
        if (buffer.size() > 0) {
          for (int i = 0; i < totalSize; i++) {
            int z = permutation[i] / (size[0] * size[1]);
            int y = (permutation[i] / size[0]) % size[1];
            int x = permutation[i] % size[0];
            buffer[i] = field[x][y][z];
          }
        }
        return totalSize;
      }

      //Deserialize 3D Pointers
      template <typename U = fieldType>
      typename std::enable_if<std::is_pointer<remove2Pointers<U>>::value 
                          && !std::is_pointer<remove3Pointers<U>>::value, void>::type
      Deserialize(pcms::ScalarArrayView<const value_type, memory_space> buffer,
                       pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
      {
        REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
        for (int i = 0; i < buffer.size(); ++i) {
          int z = permutation[i] / (size[0] * size[1]);
          int y = (permutation[i] / size[0]) % size[1];
          int x = permutation[i] % size[0];
          field[x][y][z] = buffer[i];
        }
      }

    private:
      fieldType field;
      int size[3];
      int totalSize;
      std::vector<pcms::GO> gids;
  };
}

#endif //PCMS_LIB_H

#endif //PCMS_ENABLED