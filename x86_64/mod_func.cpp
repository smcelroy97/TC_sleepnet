#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern "C" void _Ih_reg(void);
extern "C" void _ampa_reg(void);
extern "C" void _ampa_D1_reg(void);
extern "C" void _ampa_D2_reg(void);
extern "C" void _ampa_NEURON_reg(void);
extern "C" void _cadecay_reg(void);
extern "C" void _gaba_A_reg(void);
extern "C" void _gaba_A_D2_reg(void);
extern "C" void _gaba_B_reg(void);
extern "C" void _hva_reg(void);
extern "C" void _iT_RE_reg(void);
extern "C" void _iT_TC_reg(void);
extern "C" void _kL_reg(void);
extern "C" void _kca_reg(void);
extern "C" void _kdr_reg(void);
extern "C" void _kdr_re_reg(void);
extern "C" void _kdr_tc_reg(void);
extern "C" void _km_reg(void);
extern "C" void _naf_reg(void);
extern "C" void _naf_re_reg(void);
extern "C" void _naf_tc_reg(void);
extern "C" void _nap_reg(void);
extern "C" void _nmda_D1_reg(void);
extern "C" void _no_diffusion_reg(void);
extern "C" void _xtra_reg(void);

extern "C" void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"mod/Ih.mod\"");
    fprintf(stderr, " \"mod/ampa.mod\"");
    fprintf(stderr, " \"mod/ampa_D1.mod\"");
    fprintf(stderr, " \"mod/ampa_D2.mod\"");
    fprintf(stderr, " \"mod/ampa_NEURON.mod\"");
    fprintf(stderr, " \"mod/cadecay.mod\"");
    fprintf(stderr, " \"mod/gaba_A.mod\"");
    fprintf(stderr, " \"mod/gaba_A_D2.mod\"");
    fprintf(stderr, " \"mod/gaba_B.mod\"");
    fprintf(stderr, " \"mod/hva.mod\"");
    fprintf(stderr, " \"mod/iT_RE.mod\"");
    fprintf(stderr, " \"mod/iT_TC.mod\"");
    fprintf(stderr, " \"mod/kL.mod\"");
    fprintf(stderr, " \"mod/kca.mod\"");
    fprintf(stderr, " \"mod/kdr.mod\"");
    fprintf(stderr, " \"mod/kdr_re.mod\"");
    fprintf(stderr, " \"mod/kdr_tc.mod\"");
    fprintf(stderr, " \"mod/km.mod\"");
    fprintf(stderr, " \"mod/naf.mod\"");
    fprintf(stderr, " \"mod/naf_re.mod\"");
    fprintf(stderr, " \"mod/naf_tc.mod\"");
    fprintf(stderr, " \"mod/nap.mod\"");
    fprintf(stderr, " \"mod/nmda_D1.mod\"");
    fprintf(stderr, " \"mod/no_diffusion.mod\"");
    fprintf(stderr, " \"mod/xtra.mod\"");
    fprintf(stderr, "\n");
  }
  _Ih_reg();
  _ampa_reg();
  _ampa_D1_reg();
  _ampa_D2_reg();
  _ampa_NEURON_reg();
  _cadecay_reg();
  _gaba_A_reg();
  _gaba_A_D2_reg();
  _gaba_B_reg();
  _hva_reg();
  _iT_RE_reg();
  _iT_TC_reg();
  _kL_reg();
  _kca_reg();
  _kdr_reg();
  _kdr_re_reg();
  _kdr_tc_reg();
  _km_reg();
  _naf_reg();
  _naf_re_reg();
  _naf_tc_reg();
  _nap_reg();
  _nmda_D1_reg();
  _no_diffusion_reg();
  _xtra_reg();
}
