#include "fusion_io.h"
#include <stdio.h>
#include <iostream>

mars_source::mars_eqdata::mars_eqdata(const int nri, const int nci)
{
  cse = new double[nri];
  peq = new double[nri];
  t = new double[nri];
  ttp = new double[nri];
  ppeq= new double[nri];
  dpsids = new double[nri];
  req = new double*[nri];
  zeq = new double*[nri];
  rja = new double*[nri];
  g11l = new double*[nri];
  g22l = new double*[nri];
  g33l = new double*[nri];
  g12l = new double*[nri];
  rdcdz = new double*[nri];
  rdsdz = new double*[nri];
  rbz = new double*[nri];
  for(int i=0; i<nri; i++) {
    req[i] = new double[nci];
    zeq[i] = new double[nci];
    rja[i] = new double[nci];
    g11l[i] = new double[nci];
    g22l[i] = new double[nci];
    g33l[i] = new double[nci];
    g12l[i] = new double[nci];
    rdcdz[i] = new double[nci];
    rdsdz[i] = new double[nci];
    rbz[i] = new double[nci];
  }
  tp = new double[nri];

  nr = nri;
  nc = nci;
}

mars_source::mars_eqdata::~mars_eqdata()
{
  if(nr==0 || nc==0) return;
  
  delete[] cse;
  delete[] peq;
  delete[] t;
  delete[] ttp;
  delete[] dpsids;
  for(int i=0; i<nr; i++) {
    delete[] req[i];
    delete[] zeq[i];
    delete[] rja[i];
    delete[] g11l[i];
    delete[] g22l[i];
    delete[] g33l[i];
    delete[] g12l[i];
    delete[] rdcdz[i];
    delete[] rdsdz[i];
    delete[] rbz[i];
  }
  delete[] tp;
}

mars_source::mars_eigdata::mars_eigdata(const int m, const int n, 
					const bool ex)
{
  for(int i=0; i<3; i++) {
    if(ex) {
      dpsids[i] = new double[n];
      t[i] = new double[n];
    }
    mnum_r[i] = new double[m];
    mnum_i[i] = new double[m];
    val_r[i] = new double*[m];
    val_i[i] = new double*[m];
    for(int j=0; j<m; j++) {
      val_r[i][j] = new double[n];
      val_i[i][j] = new double[n];
    }
  }
  
  extra = ex;
  maxm = m;
  nrp = n;
}

mars_source::mars_eigdata::~mars_eigdata()
{
  if(maxm==0 || nrp==0) return;

  for(int i=0; i<3; i++) {
    if(extra) {
      delete[] dpsids[i];
      delete[] t[i];
    }
    delete[] mnum_r[i];
    delete[] mnum_i[i];
    for(int j=0; j<maxm; j++) {
      delete[] val_r[i][j];
      delete[] val_i[i][j];
    }
    delete[] val_r[i];
    delete[] val_i[i];
  }
}

mars_source::mars_source() 
  : psiiso(0), eqdata(0), eqdatam(0), xplasma(0), vplasma(0), bplasma(0)
{
}

mars_source::~mars_source()
{
  if(bplasma) delete(bplasma);
  if(xplasma) delete(xplasma);
  if(vplasma) delete(vplasma);

  close();
}

int mars_source::open(const char* )
{
  // Read plasma equilibrium
  std::cerr << "Reading OUTRMAR file" << std::endl;

  FILE* fid = fopen("OUTRMAR", "r");
  if(fid==0) return FIO_FILE_ERROR;
  fscanf(fid, "%i", &nrp1);
  fscanf(fid, "%i", &nchi);
  fscanf(fid, "%le", &aspect);
  fscanf(fid, "%le", &r0exp);
  fscanf(fid, "%le", &b0exp);

  std::cerr << "NRP1   = " << nrp1 << "\n";
  std::cerr << "NCHI   = " << nchi << "\n";
  std::cerr << "ASPECT = " << aspect << "\n";
  std::cerr << "R0EXP  = " << r0exp << "\n";
  std::cerr << "B0EXP  = " << b0exp << "\n";

  bfac = 4.e-7*M_PI/r0exp;
  vfac = 4.e-7*M_PI/b0exp;

  std::cerr << "bfac = " << bfac << " T\n";
  std::cerr << "vfac = " << vfac << " \n";


  std::cerr << "Allocating data structures for MARS data" << std::endl;
  psiiso = new double[2*nrp1];
  eqdata = new mars_eqdata(nrp1, nchi);
  eqdatam = new mars_eqdata(nrp1, nchi);

  std::cerr << "Reading data into data structures" << std::endl;
  for(int i=0; i<nrp1; i++) {
    fscanf(fid, "%le", &(eqdatam->cse[i]));
    fscanf(fid, "%le", &psiiso[2*i]);
    fscanf(fid, "%le", &(eqdatam->peq[i]));
    fscanf(fid, "%le", &(eqdatam->t[i]));
    fscanf(fid, "%le", &(eqdatam->ttp[i]));
    fscanf(fid, "%le", &(eqdatam->ppeq[i]));
    fscanf(fid, "%le", &(eqdatam->dpsids[i]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->req[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->zeq[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->rja[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->g11l[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->g22l[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->g33l[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->g12l[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->rdcdz[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->rdsdz[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdatam->rbz[i][j]));
    eqdatam->tp[i] = eqdatam->ttp[i] / eqdatam->t[i];

    fscanf(fid, "%le", &(eqdata->cse[i]));
    fscanf(fid, "%le", &psiiso[2*i+1]);
    fscanf(fid, "%le", &(eqdata->peq[i]));
    fscanf(fid, "%le", &(eqdata->t[i]));
    fscanf(fid, "%le", &(eqdata->ttp[i]));
    fscanf(fid, "%le", &(eqdata->ppeq[i]));
    fscanf(fid, "%le", &(eqdata->dpsids[i]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->req[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->zeq[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->rja[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->g11l[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->g22l[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->g33l[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->g12l[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->rdcdz[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->rdsdz[i][j]));
    for(int j=0; j<nchi; j++) fscanf(fid, "%le", &(eqdata->rbz[i][j]));
    eqdata->tp[i] = eqdata->ttp[i] / eqdata->t[i];
  }
  fclose(fid);

  // Read vacuum equilibrium


  return FIO_SUCCESS;
}

int mars_source::close()
{
  if(psiiso) { 
    delete(psiiso);    
    psiiso = 0;
  }
  if(eqdata) {
    delete(eqdata);    
    eqdata = 0;
  }
  if(eqdatam) {
    delete(eqdatam);  
    eqdatam = 0;
  }

  return FIO_SUCCESS;
}

int mars_source::load_eigdata(const char* filename, mars_eigdata** x, 
			      const bool ex)
{
  if(*x!=0) return FIO_SUCCESS;

  double dum;
  int m, n;
  FILE* fid = fopen(filename, "r");
  if(fid==0) return FIO_FILE_ERROR;

  fscanf(fid, "%le", &dum);
  m = dum;
  fscanf(fid, "%le", &dum);
  n = dum;
  fscanf(fid, "%le", &dum);
  fscanf(fid, "%le", &dum);
  fscanf(fid, "%le", &dum);
  fscanf(fid, "%le", &dum);

  *x = new mars_eigdata(m, n, ex);

  for(int i=0; i<(*x)->maxm; i++) {
    fscanf(fid, "%le %le", &((*x)->mnum_r[0][i]), &((*x)->mnum_i[0][i]));
    fscanf(fid, "%le %le", &((*x)->mnum_r[1][i]), &((*x)->mnum_i[1][i]));
    fscanf(fid, "%le %le", &((*x)->mnum_r[2][i]), &((*x)->mnum_i[2][i]));
  }

  if((*x)->extra) {
    for(int j=0; j<(*x)->nrp; j++) {
      for(int k=0; k<3; k++) fscanf(fid, "%le", &((*x)->dpsids[k][j]));
      for(int k=0; k<3; k++) fscanf(fid, "%le", &((*x)->t[k][j]));
    }
  }

  for(int i=0; i<(*x)->maxm; i++) {
    for(int j=0; j<(*x)->nrp; j++) {
      for(int k=0; k<3; k++) {
	fscanf(fid, "%le", &((*x)->val_r[k][i][j]));
	fscanf(fid, "%le", &((*x)->val_r[k][i][j]));
      }
    }
  }
  
  fclose(fid);

  return FIO_SUCCESS;
}

int mars_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();

  opt->add_option(FIO_TIMESLICE, 0);
  opt->add_option(FIO_LINEAR_SCALE, 1.);
  opt->add_option(FIO_PHASE, 0.);
  opt->add_option(FIO_PART, FIO_TOTAL);

  return FIO_SUCCESS;
}


int mars_source::get_field(const field_type t, fio_field** f,
			   const fio_option_list* opt)
{
  *f = 0;
  mars_field* mf;
  int result;

  result = FIO_SUCCESS;

  switch(t) {

  case(FIO_FLUID_VELOCITY):
    mf = new mars_fluid_velocity(this);
    break;
    
  case(FIO_MAGNETIC_FIELD):
    mf = new mars_magnetic_field(this);
    break;
    
  default:
    return FIO_UNSUPPORTED;
  };

  result = mf->load(opt);
  if(result == FIO_SUCCESS) {
    *f = mf;
  } else {
    delete(mf);
  }
  return result;
}

int mars_source::get_series(const series_type, fio_series**)
{
  return FIO_UNSUPPORTED;
}

int mars_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();

  return FIO_SUCCESS;
}

#ifdef PCMS_ENABLED
void mars_source::add_pcms_fields(PCMS_Library lib) {
  lib.client->AddField("psiiso", FusionIOFieldAdapter(psiiso, 2*nrp1));

  lib.client->AddField("eqdata_cse", FusionIOFieldAdapter(eqdata->cse, nrp1));
  lib.client->AddField("eqdata_peq", FusionIOFieldAdapter(eqdata->peq, nrp1));
  lib.client->AddField("eqdata_t", FusionIOFieldAdapter(eqdata->t, nrp1));
  lib.client->AddField("eqdata_ttp", FusionIOFieldAdapter(eqdata->ttp, nrp1));
  lib.client->AddField("eqdata_ppeq", FusionIOFieldAdapter(eqdata->ppeq, nrp1));
  lib.client->AddField("eqdata_dpsids", FusionIOFieldAdapter(eqdata->dpsids, nrp1));
  lib.client->AddField("eqdata_req", FusionIOFieldAdapter(eqdata->req, nrp1, nchi));
  lib.client->AddField("eqdata_zeq", FusionIOFieldAdapter(eqdata->zeq, nrp1, nchi));
  lib.client->AddField("eqdata_rja", FusionIOFieldAdapter(eqdata->rja, nrp1, nchi));
  lib.client->AddField("eqdata_g11l", FusionIOFieldAdapter(eqdata->g11l, nrp1, nchi));
  lib.client->AddField("eqdata_g22l", FusionIOFieldAdapter(eqdata->g22l, nrp1, nchi));
  lib.client->AddField("eqdata_g33l", FusionIOFieldAdapter(eqdata->g33l, nrp1, nchi));
  lib.client->AddField("eqdata_g12l", FusionIOFieldAdapter(eqdata->g12l, nrp1, nchi));
  lib.client->AddField("eqdata_rdcdz", FusionIOFieldAdapter(eqdata->rdcdz, nrp1, nchi));
  lib.client->AddField("eqdata_rdsdz", FusionIOFieldAdapter(eqdata->rdsdz, nrp1, nchi));
  lib.client->AddField("eqdata_rbz", FusionIOFieldAdapter(eqdata->rbz, nrp1, nchi));
  lib.client->AddField("eqdata_tp", FusionIOFieldAdapter(eqdata->tp, nrp1));

  lib.client->AddField("eqdatam_cse", FusionIOFieldAdapter(eqdatam->cse, nrp1));
  lib.client->AddField("eqdatam_peq", FusionIOFieldAdapter(eqdatam->peq, nrp1));
  lib.client->AddField("eqdatam_t", FusionIOFieldAdapter(eqdatam->t, nrp1));
  lib.client->AddField("eqdatam_ttp", FusionIOFieldAdapter(eqdatam->ttp, nrp1));
  lib.client->AddField("eqdatam_ppeq", FusionIOFieldAdapter(eqdatam->ppeq, nrp1));
  lib.client->AddField("eqdatam_dpsids", FusionIOFieldAdapter(eqdatam->dpsids, nrp1));
  lib.client->AddField("eqdatam_req", FusionIOFieldAdapter(eqdatam->req, nrp1, nchi));
  lib.client->AddField("eqdatam_zeq", FusionIOFieldAdapter(eqdatam->zeq, nrp1, nchi));
  lib.client->AddField("eqdatam_rja", FusionIOFieldAdapter(eqdatam->rja, nrp1, nchi));
  lib.client->AddField("eqdatam_g11l", FusionIOFieldAdapter(eqdatam->g11l, nrp1, nchi));
  lib.client->AddField("eqdatam_g22l", FusionIOFieldAdapter(eqdatam->g22l, nrp1, nchi));
  lib.client->AddField("eqdatam_g33l", FusionIOFieldAdapter(eqdatam->g33l, nrp1, nchi));
  lib.client->AddField("eqdatam_g12l", FusionIOFieldAdapter(eqdatam->g12l, nrp1, nchi));
  lib.client->AddField("eqdatam_rdcdz", FusionIOFieldAdapter(eqdatam->rdcdz, nrp1, nchi));
  lib.client->AddField("eqdatam_rdsdz", FusionIOFieldAdapter(eqdatam->rdsdz, nrp1, nchi));
  lib.client->AddField("eqdatam_rbz", FusionIOFieldAdapter(eqdatam->rbz, nrp1, nchi));
  lib.client->AddField("eqdatam_tp", FusionIOFieldAdapter(eqdatam->tp, nrp1));

  if(xplasma->extra) {
    lib.client->AddField("xplasma_dpsids", FusionIOFieldAdapter(xplasma->dpsids, 3, xplasma->nrp));
    lib.client->AddField("xplasma_t", FusionIOFieldAdapter(xplasma->t, 3, xplasma->nrp));
  }
  lib.client->AddField("xplasma_mnum_r", FusionIOFieldAdapter(xplasma->mnum_r, 3, xplasma->maxm));
  lib.client->AddField("xplasma_mnum_i", FusionIOFieldAdapter(xplasma->mnum_i, 3, xplasma->maxm));
  lib.client->AddField("xplasma_val_r", FusionIOFieldAdapter(xplasma->val_r, 3, xplasma->maxm, xplasma->nrp));
  lib.client->AddField("xplasma_val_i", FusionIOFieldAdapter(xplasma->val_i, 3, xplasma->maxm, xplasma->nrp));

  if(vplasma->extra) {
    lib.client->AddField("vplasma_dpsids", FusionIOFieldAdapter(vplasma->dpsids, 3, vplasma->nrp));
    lib.client->AddField("vplasma_t", FusionIOFieldAdapter(vplasma->t, 3, vplasma->nrp));
  }
  lib.client->AddField("vplasma_mnum_r", FusionIOFieldAdapter(vplasma->mnum_r, 3, vplasma->maxm));
  lib.client->AddField("vplasma_mnum_i", FusionIOFieldAdapter(vplasma->mnum_i, 3, vplasma->maxm));
  lib.client->AddField("vplasma_val_r", FusionIOFieldAdapter(vplasma->val_r, 3, vplasma->maxm, vplasma->nrp));
  lib.client->AddField("vplasma_val_i", FusionIOFieldAdapter(vplasma->val_i, 3, vplasma->maxm, vplasma->nrp));

  if(bplasma->extra) {
    lib.client->AddField("bplasma_dpsids", FusionIOFieldAdapter(bplasma->dpsids, 3, bplasma->nrp));
    lib.client->AddField("bplasma_t", FusionIOFieldAdapter(bplasma->t, 3, bplasma->nrp));
  }
  lib.client->AddField("bplasma_mnum_r", FusionIOFieldAdapter(bplasma->mnum_r, 3, bplasma->maxm));
  lib.client->AddField("bplasma_mnum_i", FusionIOFieldAdapter(bplasma->mnum_i, 3, bplasma->maxm));
  lib.client->AddField("bplasma_val_r", FusionIOFieldAdapter(bplasma->val_r, 3, bplasma->maxm, bplasma->nrp));
  lib.client->AddField("bplasma_val_i", FusionIOFieldAdapter(bplasma->val_i, 3, bplasma->maxm, bplasma->nrp));
}
#endif //PCMS_ENABLED
