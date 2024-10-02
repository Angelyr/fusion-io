#pragma once
#include "fusion_io_source.h"
#include "fusion_io.h"
#include <pcms/server.h>
#include <pcms/client.h>
#include <pcms/field.h>
#include <pcms/types.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <string.h>


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

namespace fusion_io
{
  class Library
  {
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
  class FusionFieldAdapter
  {
    public:
      using value_type = std::remove_pointer_t<std::remove_pointer_t<fieldtype>>;
      fieldtype* field;
      int* size;

      FusionFieldAdapter(fieldtype* fieldIn) {
        field = fieldIn;
      }

      template <typename T1, typename T2>
      int Serialize(T1, T2) const noexcept { return 0;}

      template <typename T1, typename T2>
      void Deserialize(T1, T2) const noexcept {}

      [[nodiscard]] std::vector<pcms::GO> GetGids() const { return {}; }

      [[nodiscard]] pcms::ReversePartitionMap GetReversePartitionMap(const redev::Partition& partition) const {
        return {};
      }

      template < typename = typename std::enable_if< IndirectionLevel<fieldtype>::value==1 > >
      int Serialize(pcms::ScalarArrayView<fieldtype, pcms::HostMemorySpace> buffer,
                    pcms::ScalarArrayView<const pcms::LO, pcms::HostMemorySpace> permutation) const
      {
        if (buffer.size() > 0) {
          for (pcms::LO i = 0; i < *size; i++) {
            buffer[i] = field[permutation[i]];
          }
        }
        return *size;
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

  };

  Omega_h::HostRead<double> evaluate(fio_field* field, double* x, double* values) {
    field->eval(x, values);
    Omega_h::HostWrite<double> host_write(field->dimension());
    for (int i=0; i<field->dimension(); i++) {
      host_write[i] = values[i];
    }
    Omega_h::HostRead<double> host_read(host_write.write());
    return host_read;
  }

  int fio_open_source(fio_source** src, const int type, const char* filename, Library lib) {
    int ierr;
    *src = 0;

    switch(type) {
    case(FIO_GATO_SOURCE):
    {
      gato_source* tgt_src = new gato_source();
      ierr = tgt_src->open(filename);
      *src = tgt_src;

      lib.client->AddField("psival", FusionFieldAdapter(tgt_src->psival));
      lib.client->AddField("pressure", FusionFieldAdapter(tgt_src->pressure));
      lib.client->AddField("ftor", FusionFieldAdapter(tgt_src->ftor));
      lib.client->AddField("pprime", FusionFieldAdapter(tgt_src->pprime));
      lib.client->AddField("ffprime", FusionFieldAdapter(tgt_src->ffprime));
      lib.client->AddField("rcc", FusionFieldAdapter(tgt_src->rcc));
      lib.client->AddField("zcc", FusionFieldAdapter(tgt_src->zcc));
      lib.client->AddField("psimesh", FusionFieldAdapter(tgt_src->psimesh));
      lib.client->AddField("dpsidr", FusionFieldAdapter(tgt_src->dpsidr));
      lib.client->AddField("dpsidz", FusionFieldAdapter(tgt_src->dpsidz));
    }
    break;

    case(FIO_GEQDSK_SOURCE):
    {
      geqdsk_source* tgt_src = new geqdsk_source();
      ierr = tgt_src->open(filename);
      *src = tgt_src;
      
      lib.client->AddField("rmaxis", FusionFieldAdapter(&(tgt_src->rmaxis)));
      lib.client->AddField("zmaxis", FusionFieldAdapter(&(tgt_src->zmaxis)));
      lib.client->AddField("nw", FusionFieldAdapter(&(tgt_src->nw)));
      lib.client->AddField("nh", FusionFieldAdapter(&(tgt_src->nh)));
      lib.client->AddField("dx", FusionFieldAdapter(&(tgt_src->dx)));
      lib.client->AddField("dz", FusionFieldAdapter(&(tgt_src->dz)));
      lib.client->AddField("rleft", FusionFieldAdapter(&(tgt_src->rleft)));
      lib.client->AddField("zmid", FusionFieldAdapter(&(tgt_src->zmid)));
      lib.client->AddField("zbottom", FusionFieldAdapter(&(tgt_src->zbottom)));
      lib.client->AddField("simag", FusionFieldAdapter(&(tgt_src->simag)));
      lib.client->AddField("sibry", FusionFieldAdapter(&(tgt_src->sibry)));
      lib.client->AddField("psi", FusionFieldAdapter(tgt_src->psi));
      lib.client->AddField("psirz", FusionFieldAdapter(tgt_src->psirz));
      lib.client->AddField("fpol", FusionFieldAdapter(tgt_src->fpol));
      lib.client->AddField("ffprime", FusionFieldAdapter(tgt_src->ffprime));
      lib.client->AddField("press", FusionFieldAdapter(tgt_src->press));
    }
    break;

    case(FIO_GPEC_SOURCE):
    {
      gpec_source* tgt_src = new gpec_source();
      ierr = tgt_src->open(filename);
      *src = tgt_src;

      lib.client->AddField("b0_ntor", FusionFieldAdapter(&(tgt_src->b0.ntor)));
      lib.client->AddField("b0_nr", FusionFieldAdapter(&(tgt_src->b0.nr)));
      lib.client->AddField("b0_nz", FusionFieldAdapter(&(tgt_src->b0.nz)));
      lib.client->AddField("b0_n_comp", FusionFieldAdapter(&(tgt_src->b0.n_comp)));
      lib.client->AddField("b0_r", FusionFieldAdapter(tgt_src->b0.r));
      lib.client->AddField("b0_z", FusionFieldAdapter(tgt_src->b0.z));
      lib.client->AddField("b0_v_real", FusionFieldAdapter(tgt_src->b0.v_real));
      lib.client->AddField("b0_v_imag", FusionFieldAdapter(tgt_src->b0.v_imag));

      lib.client->AddField("b1_ntor", FusionFieldAdapter(&(tgt_src->b1.ntor)));
      lib.client->AddField("b1_nr", FusionFieldAdapter(&(tgt_src->b1.nr)));
      lib.client->AddField("b1_nz", FusionFieldAdapter(&(tgt_src->b1.nz)));
      lib.client->AddField("b1_n_comp", FusionFieldAdapter(&(tgt_src->b1.n_comp)));
      lib.client->AddField("b1_r", FusionFieldAdapter(tgt_src->b1.r));
      lib.client->AddField("b1_z", FusionFieldAdapter(tgt_src->b1.z));
      lib.client->AddField("b1_v_real", FusionFieldAdapter(tgt_src->b1.v_real));
      lib.client->AddField("b1_v_imag", FusionFieldAdapter(tgt_src->b1.v_imag));

      lib.client->AddField("bx_ntor", FusionFieldAdapter(&(tgt_src->bx.ntor)));
      lib.client->AddField("bx_nr", FusionFieldAdapter(&(tgt_src->bx.nr)));
      lib.client->AddField("bx_nz", FusionFieldAdapter(&(tgt_src->bx.nz)));
      lib.client->AddField("bx_n_comp", FusionFieldAdapter(&(tgt_src->bx.n_comp)));
      lib.client->AddField("bx_r", FusionFieldAdapter(tgt_src->bx.r));
      lib.client->AddField("bx_z", FusionFieldAdapter(tgt_src->bx.z));
      lib.client->AddField("bx_v_real", FusionFieldAdapter(tgt_src->bx.v_real));
      lib.client->AddField("bx_v_imag", FusionFieldAdapter(tgt_src->bx.v_imag));
    } 
    break;

    case(FIO_M3DC1_SOURCE):
    {
      m3dc1_source* tgt_src = new m3dc1_source();
      ierr = tgt_src->open(filename);
      *src = tgt_src;

      lib.client->AddField("bzero", FusionFieldAdapter(&(tgt_src->bzero)));
      lib.client->AddField("rzero", FusionFieldAdapter(&(tgt_src->rzero)));
      lib.client->AddField("z_ion", FusionFieldAdapter(&(tgt_src->z_ion)));
      lib.client->AddField("ion_mass", FusionFieldAdapter(&(tgt_src->ion_mass)));
      lib.client->AddField("period", FusionFieldAdapter(&(tgt_src->period)));
      lib.client->AddField("n0", FusionFieldAdapter(&(tgt_src->n0)));
      lib.client->AddField("L0", FusionFieldAdapter(&(tgt_src->L0)));
      lib.client->AddField("B0", FusionFieldAdapter(&(tgt_src->B0)));
      lib.client->AddField("p0", FusionFieldAdapter(&(tgt_src->p0)));
      lib.client->AddField("t0", FusionFieldAdapter(&(tgt_src->t0)));
      lib.client->AddField("v0", FusionFieldAdapter(&(tgt_src->v0)));
      lib.client->AddField("J0", FusionFieldAdapter(&(tgt_src->J0)));
      lib.client->AddField("Phi0", FusionFieldAdapter(&(tgt_src->Phi0)));
      lib.client->AddField("temp0", FusionFieldAdapter(&(tgt_src->temp0)));
      lib.client->AddField("eta0", FusionFieldAdapter(&(tgt_src->eta0)));
      lib.client->AddField("linear", FusionFieldAdapter(&(tgt_src->linear)));
      lib.client->AddField("eqsubtract", FusionFieldAdapter(&(tgt_src->eqsubtract)));
      lib.client->AddField("extsubtract", FusionFieldAdapter(&(tgt_src->extsubtract)));
      lib.client->AddField("icomplex", FusionFieldAdapter(&(tgt_src->icomplex)));
      lib.client->AddField("i3d", FusionFieldAdapter(&(tgt_src->i3d)));
      lib.client->AddField("version", FusionFieldAdapter(&(tgt_src->version)));
      lib.client->AddField("itor", FusionFieldAdapter(&(tgt_src->itor)));
      lib.client->AddField("ntor", FusionFieldAdapter(&(tgt_src->ntor)));
      lib.client->AddField("ntime", FusionFieldAdapter(&(tgt_src->ntime)));
      lib.client->AddField("ifprime", FusionFieldAdapter(&(tgt_src->ifprime)));
      lib.client->AddField("kprad_z", FusionFieldAdapter(&(tgt_src->kprad_z)));
      lib.client->AddField("numvar", FusionFieldAdapter(&(tgt_src->numvar)));
    }
    break;

    case(FIO_MARS_SOURCE):
    {
      mars_source* tgt_src = new mars_source();
      ierr = tgt_src->open(filename);
      *src = tgt_src;

      lib.client->AddField("nrp1", FusionFieldAdapter(&(tgt_src->nrp1)));
      lib.client->AddField("nchi", FusionFieldAdapter(&(tgt_src->nchi)));
      lib.client->AddField("aspect", FusionFieldAdapter(&(tgt_src->aspect)));
      lib.client->AddField("r0exp", FusionFieldAdapter(&(tgt_src->r0exp)));
      lib.client->AddField("b0exp", FusionFieldAdapter(&(tgt_src->b0exp)));
      lib.client->AddField("bfac", FusionFieldAdapter(&(tgt_src->bfac)));
      lib.client->AddField("vfac", FusionFieldAdapter(&(tgt_src->vfac)));
      lib.client->AddField("psiiso", FusionFieldAdapter(tgt_src->psiiso));

      lib.client->AddField("eqdata_cse", FusionFieldAdapter(tgt_src->eqdata->cse));
      lib.client->AddField("eqdata_peq", FusionFieldAdapter(tgt_src->eqdata->peq));
      lib.client->AddField("eqdata_t", FusionFieldAdapter(tgt_src->eqdata->t));
      lib.client->AddField("eqdata_ttp", FusionFieldAdapter(tgt_src->eqdata->ttp));
      lib.client->AddField("eqdata_ppeq", FusionFieldAdapter(tgt_src->eqdata->ppeq));
      lib.client->AddField("eqdata_dpsids", FusionFieldAdapter(tgt_src->eqdata->dpsids));
      lib.client->AddField("eqdata_req", FusionFieldAdapter(tgt_src->eqdata->req));
      lib.client->AddField("eqdata_zeq", FusionFieldAdapter(tgt_src->eqdata->zeq));
      lib.client->AddField("eqdata_rja", FusionFieldAdapter(tgt_src->eqdata->rja));
      lib.client->AddField("eqdata_g11l", FusionFieldAdapter(tgt_src->eqdata->g11l));
      lib.client->AddField("eqdata_g22l", FusionFieldAdapter(tgt_src->eqdata->g22l));
      lib.client->AddField("eqdata_g33l", FusionFieldAdapter(tgt_src->eqdata->g33l));
      lib.client->AddField("eqdata_g12l", FusionFieldAdapter(tgt_src->eqdata->g12l));
      lib.client->AddField("eqdata_rdcdz", FusionFieldAdapter(tgt_src->eqdata->rdcdz));
      lib.client->AddField("eqdata_rdsdz", FusionFieldAdapter(tgt_src->eqdata->rdsdz));
      lib.client->AddField("eqdata_rbz", FusionFieldAdapter(tgt_src->eqdata->rbz));
      lib.client->AddField("eqdata_tp", FusionFieldAdapter(tgt_src->eqdata->tp));

      lib.client->AddField("eqdatam_cse", FusionFieldAdapter(tgt_src->eqdatam->cse));
      lib.client->AddField("eqdatam_peq", FusionFieldAdapter(tgt_src->eqdatam->peq));
      lib.client->AddField("eqdatam_t", FusionFieldAdapter(tgt_src->eqdatam->t));
      lib.client->AddField("eqdatam_ttp", FusionFieldAdapter(tgt_src->eqdatam->ttp));
      lib.client->AddField("eqdatam_ppeq", FusionFieldAdapter(tgt_src->eqdatam->ppeq));
      lib.client->AddField("eqdatam_dpsids", FusionFieldAdapter(tgt_src->eqdatam->dpsids));
      lib.client->AddField("eqdatam_req", FusionFieldAdapter(tgt_src->eqdatam->req));
      lib.client->AddField("eqdatam_zeq", FusionFieldAdapter(tgt_src->eqdatam->zeq));
      lib.client->AddField("eqdatam_rja", FusionFieldAdapter(tgt_src->eqdatam->rja));
      lib.client->AddField("eqdatam_g11l", FusionFieldAdapter(tgt_src->eqdatam->g11l));
      lib.client->AddField("eqdatam_g22l", FusionFieldAdapter(tgt_src->eqdatam->g22l));
      lib.client->AddField("eqdatam_g33l", FusionFieldAdapter(tgt_src->eqdatam->g33l));
      lib.client->AddField("eqdatam_g12l", FusionFieldAdapter(tgt_src->eqdatam->g12l));
      lib.client->AddField("eqdatam_rdcdz", FusionFieldAdapter(tgt_src->eqdatam->rdcdz));
      lib.client->AddField("eqdatam_rdsdz", FusionFieldAdapter(tgt_src->eqdatam->rdsdz));
      lib.client->AddField("eqdatam_rbz", FusionFieldAdapter(tgt_src->eqdatam->rbz));
      lib.client->AddField("eqdatam_tp", FusionFieldAdapter(tgt_src->eqdatam->tp));

      lib.client->AddField("bplasma_maxm", FusionFieldAdapter(&(tgt_src->bplasma->maxm)));
      lib.client->AddField("bplasma_nrp", FusionFieldAdapter(&(tgt_src->bplasma->nrp)));
      lib.client->AddField("bplasma_dpsids", FusionFieldAdapter(tgt_src->bplasma->dpsids));
      lib.client->AddField("bplasma_t", FusionFieldAdapter(tgt_src->bplasma->t));
      lib.client->AddField("bplasma_mnum_r", FusionFieldAdapter(tgt_src->bplasma->mnum_r));
      lib.client->AddField("bplasma_mnum_i", FusionFieldAdapter(tgt_src->bplasma->mnum_i));
      lib.client->AddField("bplasma_val_r", FusionFieldAdapter(tgt_src->bplasma->val_r));
      lib.client->AddField("bplasma_val_i", FusionFieldAdapter(tgt_src->bplasma->val_i));

      lib.client->AddField("vplasma_maxm", FusionFieldAdapter(&(tgt_src->vplasma->maxm)));
      lib.client->AddField("vplasma_nrp", FusionFieldAdapter(&(tgt_src->vplasma->nrp)));
      lib.client->AddField("vplasma_dpsids", FusionFieldAdapter(tgt_src->vplasma->dpsids));
      lib.client->AddField("vplasma_t", FusionFieldAdapter(tgt_src->vplasma->t));
      lib.client->AddField("vplasma_mnum_r", FusionFieldAdapter(tgt_src->vplasma->mnum_r));
      lib.client->AddField("vplasma_mnum_i", FusionFieldAdapter(tgt_src->vplasma->mnum_i));
      lib.client->AddField("vplasma_val_r", FusionFieldAdapter(tgt_src->vplasma->val_r));
      lib.client->AddField("vplasma_val_i", FusionFieldAdapter(tgt_src->vplasma->val_i));

      lib.client->AddField("bplasma_maxm", FusionFieldAdapter(&(tgt_src->bplasma->maxm)));
      lib.client->AddField("bplasma_nrp", FusionFieldAdapter(&(tgt_src->bplasma->nrp)));
      lib.client->AddField("bplasma_dpsids", FusionFieldAdapter(tgt_src->bplasma->dpsids));
      lib.client->AddField("bplasma_t", FusionFieldAdapter(tgt_src->bplasma->t));
      lib.client->AddField("bplasma_mnum_r", FusionFieldAdapter(tgt_src->bplasma->mnum_r));
      lib.client->AddField("bplasma_mnum_i", FusionFieldAdapter(tgt_src->bplasma->mnum_i));
      lib.client->AddField("bplasma_val_r", FusionFieldAdapter(tgt_src->bplasma->val_r));
      lib.client->AddField("bplasma_val_i", FusionFieldAdapter(tgt_src->bplasma->val_i));
    }
    break;

    default:
      std::cerr << "Source type " << type << " unsupported." << std::endl;
      return FIO_UNSUPPORTED;
    }

    if(ierr != FIO_SUCCESS) {
      delete(*src);
      return ierr;
    }
    
    return FIO_SUCCESS;
  }
}