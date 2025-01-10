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
using fusion_io::FieldAdapter;
void mars_source::add_pcms_fields(fusion_io::Library lib) {
  lib.client->AddField("nrp1", FieldAdapter(&nrp1));
  lib.client->AddField("nchi", FieldAdapter(&nchi));
  lib.client->AddField("aspect", FieldAdapter(&aspect));
  lib.client->AddField("r0exp", FieldAdapter(&r0exp));
  lib.client->AddField("b0exp", FieldAdapter(&b0exp));
  lib.client->AddField("bfac", FieldAdapter(&bfac));
  lib.client->AddField("vfac", FieldAdapter(&vfac));
  lib.client->AddField("psiiso", FieldAdapter(psiiso));

  lib.client->AddField("eqdata_cse", FieldAdapter(eqdata->cse));
  lib.client->AddField("eqdata_peq", FieldAdapter(eqdata->peq));
  lib.client->AddField("eqdata_t", FieldAdapter(eqdata->t));
  lib.client->AddField("eqdata_ttp", FieldAdapter(eqdata->ttp));
  lib.client->AddField("eqdata_ppeq", FieldAdapter(eqdata->ppeq));
  lib.client->AddField("eqdata_dpsids", FieldAdapter(eqdata->dpsids));
  lib.client->AddField("eqdata_req", FieldAdapter(eqdata->req));
  lib.client->AddField("eqdata_zeq", FieldAdapter(eqdata->zeq));
  lib.client->AddField("eqdata_rja", FieldAdapter(eqdata->rja));
  lib.client->AddField("eqdata_g11l", FieldAdapter(eqdata->g11l));
  lib.client->AddField("eqdata_g22l", FieldAdapter(eqdata->g22l));
  lib.client->AddField("eqdata_g33l", FieldAdapter(eqdata->g33l));
  lib.client->AddField("eqdata_g12l", FieldAdapter(eqdata->g12l));
  lib.client->AddField("eqdata_rdcdz", FieldAdapter(eqdata->rdcdz));
  lib.client->AddField("eqdata_rdsdz", FieldAdapter(eqdata->rdsdz));
  lib.client->AddField("eqdata_rbz", FieldAdapter(eqdata->rbz));
  lib.client->AddField("eqdata_tp", FieldAdapter(eqdata->tp));

  lib.client->AddField("eqdatam_cse", FieldAdapter(eqdatam->cse));
  lib.client->AddField("eqdatam_peq", FieldAdapter(eqdatam->peq));
  lib.client->AddField("eqdatam_t", FieldAdapter(eqdatam->t));
  lib.client->AddField("eqdatam_ttp", FieldAdapter(eqdatam->ttp));
  lib.client->AddField("eqdatam_ppeq", FieldAdapter(eqdatam->ppeq));
  lib.client->AddField("eqdatam_dpsids", FieldAdapter(eqdatam->dpsids));
  lib.client->AddField("eqdatam_req", FieldAdapter(eqdatam->req));
  lib.client->AddField("eqdatam_zeq", FieldAdapter(eqdatam->zeq));
  lib.client->AddField("eqdatam_rja", FieldAdapter(eqdatam->rja));
  lib.client->AddField("eqdatam_g11l", FieldAdapter(eqdatam->g11l));
  lib.client->AddField("eqdatam_g22l", FieldAdapter(eqdatam->g22l));
  lib.client->AddField("eqdatam_g33l", FieldAdapter(eqdatam->g33l));
  lib.client->AddField("eqdatam_g12l", FieldAdapter(eqdatam->g12l));
  lib.client->AddField("eqdatam_rdcdz", FieldAdapter(eqdatam->rdcdz));
  lib.client->AddField("eqdatam_rdsdz", FieldAdapter(eqdatam->rdsdz));
  lib.client->AddField("eqdatam_rbz", FieldAdapter(eqdatam->rbz));
  lib.client->AddField("eqdatam_tp", FieldAdapter(eqdatam->tp));

  lib.client->AddField("bplasma_maxm", FieldAdapter(&bplasma->maxm));
  lib.client->AddField("bplasma_nrp", FieldAdapter(&bplasma->nrp));
  lib.client->AddField("bplasma_dpsids", FieldAdapter(bplasma->dpsids));
  lib.client->AddField("bplasma_t", FieldAdapter(bplasma->t));
  lib.client->AddField("bplasma_mnum_r", FieldAdapter(bplasma->mnum_r));
  lib.client->AddField("bplasma_mnum_i", FieldAdapter(bplasma->mnum_i));
  lib.client->AddField("bplasma_val_r", FieldAdapter(bplasma->val_r));
  lib.client->AddField("bplasma_val_i", FieldAdapter(bplasma->val_i));

  lib.client->AddField("vplasma_maxm", FieldAdapter(&vplasma->maxm));
  lib.client->AddField("vplasma_nrp", FieldAdapter(&vplasma->nrp));
  lib.client->AddField("vplasma_dpsids", FieldAdapter(vplasma->dpsids));
  lib.client->AddField("vplasma_t", FieldAdapter(vplasma->t));
  lib.client->AddField("vplasma_mnum_r", FieldAdapter(vplasma->mnum_r));
  lib.client->AddField("vplasma_mnum_i", FieldAdapter(vplasma->mnum_i));
  lib.client->AddField("vplasma_val_r", FieldAdapter(vplasma->val_r));
  lib.client->AddField("vplasma_val_i", FieldAdapter(vplasma->val_i));

  lib.client->AddField("bplasma_maxm", FieldAdapter(&bplasma->maxm));
  lib.client->AddField("bplasma_nrp", FieldAdapter(&bplasma->nrp));
  lib.client->AddField("bplasma_dpsids", FieldAdapter(bplasma->dpsids));
  lib.client->AddField("bplasma_t", FieldAdapter(bplasma->t));
  lib.client->AddField("bplasma_mnum_r", FieldAdapter(bplasma->mnum_r));
  lib.client->AddField("bplasma_mnum_i", FieldAdapter(bplasma->mnum_i));
  lib.client->AddField("bplasma_val_r", FieldAdapter(bplasma->val_r));
  lib.client->AddField("bplasma_val_i", FieldAdapter(bplasma->val_i));
}
#endif //PCMS_ENABLED
