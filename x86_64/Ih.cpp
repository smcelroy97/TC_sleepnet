/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#define _pval pval
// clang-format off
#include "md1redef.h"
#include "section_fwd.hpp"
#include "nrniv_mf.h"
#include "md2redef.h"
#include "nrnconf.h"
// clang-format on
#include "neuron/cache/mechanism_range.hpp"
static constexpr auto number_of_datum_variables = 5;
static constexpr auto number_of_floating_point_variables = 26;
namespace {
template <typename T>
using _nrn_mechanism_std_vector = std::vector<T>;
using _nrn_model_sorted_token = neuron::model_sorted_token;
using _nrn_mechanism_cache_range = neuron::cache::MechanismRange<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_mechanism_cache_instance = neuron::cache::MechanismInstance<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_non_owning_id_without_container = neuron::container::non_owning_identifier_without_container;
template <typename T>
using _nrn_mechanism_field = neuron::mechanism::field<T>;
template <typename... Args>
void _nrn_mechanism_register_data_fields(Args&&... args) {
  neuron::mechanism::register_data_fields(std::forward<Args>(args)...);
}
}
 
#if !NRNGPU
#undef exp
#define exp hoc_Exp
#if NRN_ENABLE_ARCH_INDEP_EXP_POW
#undef pow
#define pow hoc_pow
#endif
#endif
 
#define nrn_init _nrn_init__iar
#define _nrn_initial _nrn_initial__iar
#define nrn_cur _nrn_cur__iar
#define _nrn_current _nrn_current__iar
#define nrn_jacob _nrn_jacob__iar
#define nrn_state _nrn_state__iar
#define _net_receive _net_receive__iar 
#define activation activation__iar 
#define evaluate_fct evaluate_fct__iar 
#define ihkin ihkin__iar 
 
#define _threadargscomma_ _ml, _iml, _ppvar, _thread, _globals, _nt,
#define _threadargsprotocomma_ Memb_list* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt,
#define _internalthreadargsprotocomma_ _nrn_mechanism_cache_range* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt,
#define _threadargs_ _ml, _iml, _ppvar, _thread, _globals, _nt
#define _threadargsproto_ Memb_list* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt
#define _internalthreadargsproto_ _nrn_mechanism_cache_range* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 
#define t _nt->_t
#define dt _nt->_dt
#define ghbar _ml->template fpfield<0>(_iml)
#define ghbar_columnindex 0
#define fac_gh_TC _ml->template fpfield<1>(_iml)
#define fac_gh_TC_columnindex 1
#define h_inf _ml->template fpfield<2>(_iml)
#define h_inf_columnindex 2
#define tau_s _ml->template fpfield<3>(_iml)
#define tau_s_columnindex 3
#define m _ml->template fpfield<4>(_iml)
#define m_columnindex 4
#define c1 _ml->template fpfield<5>(_iml)
#define c1_columnindex 5
#define o1 _ml->template fpfield<6>(_iml)
#define o1_columnindex 6
#define o2 _ml->template fpfield<7>(_iml)
#define o2_columnindex 7
#define p0 _ml->template fpfield<8>(_iml)
#define p0_columnindex 8
#define p1 _ml->template fpfield<9>(_iml)
#define p1_columnindex 9
#define eh _ml->template fpfield<10>(_iml)
#define eh_columnindex 10
#define Dc1 _ml->template fpfield<11>(_iml)
#define Dc1_columnindex 11
#define Do1 _ml->template fpfield<12>(_iml)
#define Do1_columnindex 12
#define Do2 _ml->template fpfield<13>(_iml)
#define Do2_columnindex 13
#define Dp0 _ml->template fpfield<14>(_iml)
#define Dp0_columnindex 14
#define Dp1 _ml->template fpfield<15>(_iml)
#define Dp1_columnindex 15
#define cai _ml->template fpfield<16>(_iml)
#define cai_columnindex 16
#define ih _ml->template fpfield<17>(_iml)
#define ih_columnindex 17
#define gh _ml->template fpfield<18>(_iml)
#define gh_columnindex 18
#define alpha _ml->template fpfield<19>(_iml)
#define alpha_columnindex 19
#define beta _ml->template fpfield<20>(_iml)
#define beta_columnindex 20
#define k1ca _ml->template fpfield<21>(_iml)
#define k1ca_columnindex 21
#define k3p _ml->template fpfield<22>(_iml)
#define k3p_columnindex 22
#define tadj _ml->template fpfield<23>(_iml)
#define tadj_columnindex 23
#define v _ml->template fpfield<24>(_iml)
#define v_columnindex 24
#define _g _ml->template fpfield<25>(_iml)
#define _g_columnindex 25
#define _ion_eh *(_ml->dptr_field<0>(_iml))
#define _p_ion_eh static_cast<neuron::container::data_handle<double>>(_ppvar[0])
#define _ion_ih *(_ml->dptr_field<1>(_iml))
#define _p_ion_ih static_cast<neuron::container::data_handle<double>>(_ppvar[1])
#define _ion_dihdv *(_ml->dptr_field<2>(_iml))
#define _ion_cai *(_ml->dptr_field<3>(_iml))
#define _p_ion_cai static_cast<neuron::container::data_handle<double>>(_ppvar[3])
#define _ion_cao *(_ml->dptr_field<4>(_iml))
#define _p_ion_cao static_cast<neuron::container::data_handle<double>>(_ppvar[4])
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 static Prop* _extcall_prop;
 /* _prop_id kind of shadows _extcall_prop to allow validity checking. */
 static _nrn_non_owning_id_without_container _prop_id{};
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_activation(void);
 static void _hoc_evaluate_fct(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 static void _hoc_setdata();
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {"setdata_iar", _hoc_setdata},
 {"activation_iar", _hoc_activation},
 {"evaluate_fct_iar", _hoc_evaluate_fct},
 {0, 0}
};
 
/* Direct Python call wrappers to density mechanism functions.*/
 static double _npy_activation(Prop*);
 static double _npy_evaluate_fct(Prop*);
 
static NPyDirectMechFunc npy_direct_func_proc[] = {
 {"activation", _npy_activation},
 {"evaluate_fct", _npy_evaluate_fct},
 {0, 0}
};
 /* declare global and static user variables */
 #define gind 0
 #define _gth 0
#define Pc Pc_iar
 double Pc = 0.007;
#define cac cac_iar
 double cac = 0.0015;
#define ginc ginc_iar
 double ginc = 2;
#define k4 k4_iar
 double k4 = 0.001;
#define k2 k2_iar
 double k2 = 0.0004;
#define nexp nexp_iar
 double nexp = 1;
#define nca nca_iar
 double nca = 4;
#define taum taum_iar
 double taum = 20;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"cac_iar", "mM"},
 {"k2_iar", "1/ms"},
 {"k4_iar", "1/ms"},
 {"taum_iar", "ms"},
 {"ghbar_iar", "mho/cm2"},
 {"fac_gh_TC_iar", "mV"},
 {"tau_s_iar", "ms"},
 {0, 0}
};
 static double c10 = 0;
 static double delta_t = 0.01;
 static double o20 = 0;
 static double o10 = 0;
 static double p10 = 0;
 static double p00 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"cac_iar", &cac_iar},
 {"k2_iar", &k2_iar},
 {"Pc_iar", &Pc_iar},
 {"k4_iar", &k4_iar},
 {"nca_iar", &nca_iar},
 {"nexp_iar", &nexp_iar},
 {"ginc_iar", &ginc_iar},
 {"taum_iar", &taum_iar},
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 _prop_id = _nrn_get_prop_id(_prop);
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 static void nrn_alloc(Prop*);
static void nrn_init(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_state(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 static void nrn_cur(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_jacob(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(Prop*, int, neuron::container::data_handle<double>*, neuron::container::data_handle<double>*, double*, int);
static void _ode_spec(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void _ode_matsol(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
#define _cvode_ieq _ppvar[5].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"iar",
 "ghbar_iar",
 "fac_gh_TC_iar",
 0,
 "h_inf_iar",
 "tau_s_iar",
 "m_iar",
 0,
 "c1_iar",
 "o1_iar",
 "o2_iar",
 "p0_iar",
 "p1_iar",
 0,
 0};
 static Symbol* _h_sym;
 static Symbol* _ca_sym;
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     1.7e-05, /* ghbar */
     0, /* fac_gh_TC */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
   _ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 26);
 	/*initialize range parameters*/
 	ghbar = _parm_default[0]; /* 1.7e-05 */
 	fac_gh_TC = _parm_default[1]; /* 0 */
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 26);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_h_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0] = _nrn_mechanism_get_param_handle(prop_ion, 0); /* eh */
 	_ppvar[1] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ih */
 	_ppvar[2] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dihdv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3] = _nrn_mechanism_get_param_handle(prop_ion, 1); /* cai */
 	_ppvar[4] = _nrn_mechanism_get_param_handle(prop_ion, 2); /* cao */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 static void _thread_cleanup(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _Ih_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("h", 1.0);
 	ion_reg("ca", -10000.);
 	_h_sym = hoc_lookup("h_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 3);
  _extcall_thread.resize(2);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 hoc_register_parm_default(_mechtype, &_parm_default);
         hoc_register_npy_direct(_mechtype, npy_direct_func_proc);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"ghbar"} /* 0 */,
                                       _nrn_mechanism_field<double>{"fac_gh_TC"} /* 1 */,
                                       _nrn_mechanism_field<double>{"h_inf"} /* 2 */,
                                       _nrn_mechanism_field<double>{"tau_s"} /* 3 */,
                                       _nrn_mechanism_field<double>{"m"} /* 4 */,
                                       _nrn_mechanism_field<double>{"c1"} /* 5 */,
                                       _nrn_mechanism_field<double>{"o1"} /* 6 */,
                                       _nrn_mechanism_field<double>{"o2"} /* 7 */,
                                       _nrn_mechanism_field<double>{"p0"} /* 8 */,
                                       _nrn_mechanism_field<double>{"p1"} /* 9 */,
                                       _nrn_mechanism_field<double>{"eh"} /* 10 */,
                                       _nrn_mechanism_field<double>{"Dc1"} /* 11 */,
                                       _nrn_mechanism_field<double>{"Do1"} /* 12 */,
                                       _nrn_mechanism_field<double>{"Do2"} /* 13 */,
                                       _nrn_mechanism_field<double>{"Dp0"} /* 14 */,
                                       _nrn_mechanism_field<double>{"Dp1"} /* 15 */,
                                       _nrn_mechanism_field<double>{"cai"} /* 16 */,
                                       _nrn_mechanism_field<double>{"ih"} /* 17 */,
                                       _nrn_mechanism_field<double>{"gh"} /* 18 */,
                                       _nrn_mechanism_field<double>{"alpha"} /* 19 */,
                                       _nrn_mechanism_field<double>{"beta"} /* 20 */,
                                       _nrn_mechanism_field<double>{"k1ca"} /* 21 */,
                                       _nrn_mechanism_field<double>{"k3p"} /* 22 */,
                                       _nrn_mechanism_field<double>{"tadj"} /* 23 */,
                                       _nrn_mechanism_field<double>{"v"} /* 24 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 25 */,
                                       _nrn_mechanism_field<double*>{"_ion_eh", "h_ion"} /* 0 */,
                                       _nrn_mechanism_field<double*>{"_ion_ih", "h_ion"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"_ion_dihdv", "h_ion"} /* 2 */,
                                       _nrn_mechanism_field<double*>{"_ion_cai", "ca_ion"} /* 3 */,
                                       _nrn_mechanism_field<double*>{"_ion_cao", "ca_ion"} /* 4 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 5 */);
  hoc_register_prop_size(_mechtype, 26, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "h_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "h_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "h_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 iar /Users/scoot/TC_sleepnet/mod/Ih.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "anomalous rectifier channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int activation(_internalthreadargsprotocomma_ double, double);
static int evaluate_fct(_internalthreadargsprotocomma_ double, double);
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(static_cast<SparseObj*>(_so), _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[5], _dlist1[5]; static double *_temp1;
 static int ihkin (void* _so, double* _rhs, _internalthreadargsproto_);
 
static int ihkin (void* _so, double* _rhs, _internalthreadargsproto_)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=2;_i<5;_i++){
  	_RHS1(_i) = -_dt1*(_ml->data(_iml, _slist1[_i]) - _ml->data(_iml, _dlist1[_i]));
	_MATELM1(_i, _i) = _dt1;
      
} }
 evaluate_fct ( _threadargscomma_ v , cai ) ;
   /* ~ c1 <-> o1 ( alpha , beta )*/
 f_flux =  alpha * c1 ;
 b_flux =  beta * o1 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  alpha ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 3 ,2)  -= _term;
 _term =  beta ;
 _MATELM1( 2 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ p0 <-> p1 ( k1ca , k2 )*/
 f_flux =  k1ca * p0 ;
 b_flux =  k2 * p1 ;
 _RHS1( 4) -= (f_flux - b_flux);
 
 _term =  k1ca ;
 _MATELM1( 4 ,4)  += _term;
 _term =  k2 ;
 _MATELM1( 4 ,1)  -= _term;
 /*REACTION*/
  /* ~ o1 <-> o2 ( k3p , k4 )*/
 f_flux =  k3p * o1 ;
 b_flux =  k4 * o2 ;
 _RHS1( 3) -= (f_flux - b_flux);
 
 _term =  k3p ;
 _MATELM1( 3 ,3)  += _term;
 _term =  k4 ;
 _MATELM1( 3 ,0)  -= _term;
 /*REACTION*/
   /* p0 + p1 = 1.0 */
 _RHS1(1) =  1.0;
 _MATELM1(1, 1) = 1;
 _RHS1(1) -= p1 ;
 _MATELM1(1, 4) = 1;
 _RHS1(1) -= p0 ;
 /*CONSERVATION*/
  /* c1 + o1 + o2 = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= o2 ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= o1 ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= c1 ;
 /*CONSERVATION*/
   } return _reset;
 }
 
static int  evaluate_fct ( _internalthreadargsprotocomma_ double _lv , double _lcai ) {
   h_inf = 1.0 / ( 1.0 + exp ( ( _lv + 75.0 + fac_gh_TC ) / 5.5 ) ) ;
   tau_s = ( taum + 1000.0 / ( exp ( ( _lv + 71.5 + fac_gh_TC ) / 14.2 ) + exp ( - ( _lv + 89.0 + fac_gh_TC ) / 11.6 ) ) ) / tadj ;
   alpha = h_inf / tau_s ;
   beta = ( 1.0 - h_inf ) / tau_s ;
   k1ca = k2 * pow( ( _lcai / cac ) , nca ) ;
   k3p = k4 * pow( ( p1 / Pc ) , nexp ) ;
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for evaluate_fct_iar. Requires prior call to setdata_iar and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
_nt = nrn_threads;
 _r = 1.;
 evaluate_fct ( _threadargscomma_ *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static double _npy_evaluate_fct(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
_nt = nrn_threads;
 _r = 1.;
 evaluate_fct ( _threadargscomma_ *getarg(1) , *getarg(2) );
 return(_r);
}
 
static int  activation ( _internalthreadargsprotocomma_ double _lv , double _lcai ) {
   double _lcc ;
 evaluate_fct ( _threadargscomma_ _lv , _lcai ) ;
   _lcc = 1.0 / ( 1.0 + pow( ( cac / _lcai ) , nca ) ) ;
   m = 1.0 / ( 1.0 + beta / alpha + pow( ( _lcc / Pc ) , nexp ) ) ;
   m = ( 1.0 + ginc * pow( ( _lcc / Pc ) , nexp ) ) * m ;
    return 0; }
 
static void _hoc_activation(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for activation_iar. Requires prior call to setdata_iar and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
_nt = nrn_threads;
 _r = 1.;
 activation ( _threadargscomma_ *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static double _npy_activation(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
_nt = nrn_threads;
 _r = 1.;
 activation ( _threadargscomma_ *getarg(1) , *getarg(2) );
 return(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(_internalthreadargsproto_) {
  int _reset=0;
  {
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<5;_i++) _ml->data(_iml, _dlist1[_i]) = 0.0;}
 evaluate_fct ( _threadargscomma_ v , cai ) ;
 /* ~ c1 <-> o1 ( alpha , beta )*/
 f_flux =  alpha * c1 ;
 b_flux =  beta * o1 ;
 Dc1 -= (f_flux - b_flux);
 Do1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ p0 <-> p1 ( k1ca , k2 )*/
 f_flux =  k1ca * p0 ;
 b_flux =  k2 * p1 ;
 Dp0 -= (f_flux - b_flux);
 Dp1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ o1 <-> o2 ( k3p , k4 )*/
 f_flux =  k3p * o1 ;
 b_flux =  k4 * o2 ;
 Do1 -= (f_flux - b_flux);
 Do2 += (f_flux - b_flux);
 
 /*REACTION*/
   /* p0 + p1 = 1.0 */
 /*CONSERVATION*/
  /* c1 + o1 + o2 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, _internalthreadargsproto_) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<5;_i++){
  	_RHS1(_i) = _dt1*(_ml->data(_iml, _dlist1[_i]));
	_MATELM1(_i, _i) = _dt1;
      
} }
 evaluate_fct ( _threadargscomma_ v , cai ) ;
 /* ~ c1 <-> o1 ( alpha , beta )*/
 _term =  alpha ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 3 ,2)  -= _term;
 _term =  beta ;
 _MATELM1( 2 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ p0 <-> p1 ( k1ca , k2 )*/
 _term =  k1ca ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 1 ,4)  -= _term;
 _term =  k2 ;
 _MATELM1( 4 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ o1 <-> o2 ( k3p , k4 )*/
 _term =  k3p ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 0 ,3)  -= _term;
 _term =  k4 ;
 _MATELM1( 3 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
   /* p0 + p1 = 1.0 */
 /*CONSERVATION*/
  /* c1 + o1 + o2 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 5;}
 
static void _ode_spec(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
   Datum* _ppvar;
   size_t _iml;   _nrn_mechanism_cache_range* _ml;   Node* _nd{};
  double _v{};
  int _cntml;
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
  _ml = &_lmr;
  _cntml = _ml_arg->_nodecount;
  Datum *_thread{_ml_arg->_thread};
  double* _globals = nullptr;
  if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _ppvar = _ml_arg->_pdata[_iml];
    _nd = _ml_arg->_nodelist[_iml];
    v = NODEV(_nd);
  eh = _ion_eh;
  cai = _ion_cai;
     _ode_spec1 (_threadargs_);
  }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  Datum* _ppvar;
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 5; ++_i) {
    _pv[_i] = _nrn_mechanism_get_param_handle(_prop, _slist1[_i]);
    _pvdot[_i] = _nrn_mechanism_get_param_handle(_prop, _dlist1[_i]);
    _cvode_abstol(_atollist, _atol, _i);
  }
 }
 
static void _ode_matsol_instance1(_internalthreadargsproto_) {
 _cvode_sparse_thread(&(_thread[_cvspth1].literal_value<void*>()), 5, _dlist1, neuron::scopmath::row_view{_ml, _iml}, _ode_matsol1, _threadargs_);
 }
 
static void _ode_matsol(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
   Datum* _ppvar;
   size_t _iml;   _nrn_mechanism_cache_range* _ml;   Node* _nd{};
  double _v{};
  int _cntml;
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
  _ml = &_lmr;
  _cntml = _ml_arg->_nodecount;
  Datum *_thread{_ml_arg->_thread};
  double* _globals = nullptr;
  if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _ppvar = _ml_arg->_pdata[_iml];
    _nd = _ml_arg->_nodelist[_iml];
    v = NODEV(_nd);
  eh = _ion_eh;
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(static_cast<SparseObj*>(_thread[_cvspth1].get<void*>()));
   _nrn_destroy_sparseobj_thread(static_cast<SparseObj*>(_thread[_spth1].get<void*>()));
 }

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  c1 = c10;
  o2 = o20;
  o1 = o10;
  p1 = p10;
  p0 = p00;
 {
   tadj = 1.0 ;
   evaluate_fct ( _threadargscomma_ v , cai ) ;
   p1 = 0.0 ;
   o1 = 0.0 ;
   o2 = 0.0 ;
   c1 = 1.0 ;
   p0 = 1.0 ;
   }
 
}
}

static void nrn_init(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 v = _v;
  eh = _ion_eh;
  cai = _ion_cai;
 initmodel(_threadargs_);
 }
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   m = o1 + ginc * o2 ;
   ih = ghbar * m * ( v - eh ) ;
   }
 _current += ih;

} return _current;
}

static void nrn_cur(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_rhs = _nt->node_rhs_storage();
auto const _vec_sav_rhs = _nt->node_sav_rhs_storage();
auto const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
  eh = _ion_eh;
  cai = _ion_cai;
 auto const _g_local = _nrn_current(_threadargscomma_ _v + .001);
 	{ double _dih;
  _dih = ih;
 _rhs = _nrn_current(_threadargscomma_ _v);
  _ion_dihdv += (_dih - ih)/.001 ;
 	}
 _g = (_g_local - _rhs)/.001;
  _ion_ih += ih ;
	 _vec_rhs[_ni[_iml]] -= _rhs;
 
}
 
}

static void nrn_jacob(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_d = _nt->node_d_storage();
auto const _vec_sav_d = _nt->node_sav_d_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
  _vec_d[_ni[_iml]] += _g;
 
}
 
}

static void nrn_state(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
_ni = _ml_arg->_nodeindices;
size_t _cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (size_t _iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
   _v = _vec_v[_ni[_iml]];
 v=_v;
{
  eh = _ion_eh;
  cai = _ion_cai;
 {  sparse_thread(&(_thread[_spth1].literal_value<void*>()), 5, _slist1, _dlist1, neuron::scopmath::row_view{_ml, _iml}, &t, dt, ihkin, _linmat1, _threadargs_);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 5; ++_i) {
      _ml->data(_iml, _slist1[_i]) += dt*_ml->data(_iml, _dlist1[_i]);
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {o2_columnindex, 0};  _dlist1[0] = {Do2_columnindex, 0};
 _slist1[1] = {p1_columnindex, 0};  _dlist1[1] = {Dp1_columnindex, 0};
 _slist1[2] = {c1_columnindex, 0};  _dlist1[2] = {Dc1_columnindex, 0};
 _slist1[3] = {o1_columnindex, 0};  _dlist1[3] = {Do1_columnindex, 0};
 _slist1[4] = {p0_columnindex, 0};  _dlist1[4] = {Dp0_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/Users/scoot/TC_sleepnet/mod/Ih.mod";
    const char* nmodl_file_text = 
  "TITLE anomalous rectifier channel\n"
  "COMMENT\n"
  ":\n"
  ": Anomalous Rectifier Ih - cation (Na/K) channel in thalamocortical neurons\n"
  ":\n"
  ": Kinetic model of calcium-induced shift in the activation of Ih channels.\n"
  ": Model of Destexhe et al., Biophys J. 65: 1538-1552, 1993, based on the\n"
  ": voltage-clamp data on the calcium dependence of If in heart cells\n"
  ": (Harigawa & Irisawa, J. Physiol. 409: 121, 1989)\n"
  ":\n"
  ": The voltage-dependence is derived from Huguenard & McCormick, \n"
  ": J Neurophysiol. 68: 1373-1383, 1992, based on voltage-clamp data of \n"
  ": McCormick & Pape, J. Physiol. 431: 291, 1990. \n"
  ":\n"
  ": Modified model of the binding of calcium through a calcium-binding (CB)\n"
  ": protein, which in turn acts on Ih channels.  This model was described in\n"
  ": detail in the following reference:\n"
  ":    Destexhe, A., Bal, T., McCormick, D.A. and Sejnowski, T.J.  Ionic \n"
  ":    mechanisms underlying synchronized oscillations and propagating waves\n"
  ":    in a model of ferret thalamic slices. Journal of Neurophysiology 76:\n"
  ":    2049-2070, 1996.\n"
  ": See also http://www.cnl.salk.edu/~alain , http://cns.fmed.ulaval.ca\n"
  ":\n"
  ":   KINETIC MODEL:\n"
  ":\n"
  ":	  Normal voltage-dependent opening of Ih channels:\n"
  ":\n"
  ":		c1 (closed) <-> o1 (open)	; rate cst alpha(V),beta(V)\n"
  ":\n"
  ":	  Ca++ binding on CB protein\n"
  ":\n"
  ":		p0 (inactive) + nca Ca <-> p1 (active)	; rate cst k1,k2\n"
  ":\n"
  ":	  Binding of active CB protein on the open form (nexp binding sites) :\n"
  ":\n"
  ":		o1 (open) + nexp p1 <-> o2 (open)	; rate cst k3,k4\n"
  ":\n"
  ":\n"
  ":   PARAMETERS:\n"
  ":	It is more useful to reformulate the parameters k1,k2 into\n"
  ":	k2 and cac = (k2/k1)^(1/nca) = half activation calcium dependence, \n"
  ":	and idem for k3,k4 into k4 and Pc = (k4/k3)^(1/nexp) = half activation\n"
  ":	of Ih binding (this is like dealing with tau_m and m_inf instead of\n"
  ":	alpha and beta in Hodgkin-Huxley equations)\n"
  ":	- k2:	this rate constant is the inverse of the real time constant of \n"
  ":             	the binding of Ca to the CB protein\n"
  ":	- cac:	the half activation (affinity) of the CB protein;\n"
  ":		around 1 to 10 microM.  \n"
  ":	- k4:	this rate constant is the inverse of the real time constant of \n"
  ":             	the binding of the CB protein to Ih channels\n"
  ":		very low: it basically governs the interspindle period\n"
  ":	- Pc:	the half activation (affinity) of the Ih channels for the\n"
  ":		CB protein;\n"
  ":	- nca:	number of binding sites of calcium on CB protein; usually 4\n"
  ":	- nexp:	number of binding sites on Ih channels\n"
  ":       - ginc: augmentation of conductance associated with the Ca bound state\n"
  ":	  (about 2-3; see Harigawa & Hirisawa, 1989)\n"
  ":\n"
  ":\n"
  ":   IMPORTANT REMARKS:\n"
  ":       - This simple model for the binding of Ca++ on the open channel \n"
  ":	  suffies to account for the shift in the voltage-dependence of Ih\n"
  ":	  activation with calcium (see details in Destexhe et al, 1993).\n"
  ":	- It may be that calcium just binds to the Ih channel, preventing the \n"
  ":	  conformational change between open and closed; in this case one\n"
  ":	  should take into account binding on the closed state, which is \n"
  ":	  neglected here.\n"
  ":\n"
  ":   MODIFICATIONS\n"
  ":	- this file also contains a procedure (\"activation\") to estimate\n"
  ":	  the steady-state activation of the current; callable from outside\n"
  ":	- the time constant now contains a changeable minimal value (taum)\n"
  ":	- shift: new local variable to displace the voltage-dependence\n"
  ":	  (shift>0 -> depolarizing shift)\n"
  ":\n"
  ":\n"
  ": Alain Destexhe, Salk Institute and Laval University, 1995\n"
  ":\n"
  ": This file uses C++ code from Bazhenov 2002, found here:\n"
  ": https://modeldb.science/28189?tab=1\n"
  ": original code modified to remove dependence on temperature \n"
  ": (for computational efficiency)\n"
  ": also replaced 'shift' with 'fac_gh_TC'\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX iar\n"
  "	USEION h READ eh WRITE ih VALENCE 1\n"
  "	USEION ca READ cai\n"
  "    RANGE ghbar, h_inf, tau_s, m, fac_gh_TC\n"
  "	GLOBAL k2, cac, k4, Pc, nca, nexp, ginc, taum\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(molar)	= (1/liter)\n"
  "	(mM)	= (millimolar)\n"
  "	(mA) 	= (milliamp)\n"
  "	(mV) 	= (millivolt)\n"
  "	(msM)	= (ms mM)\n"
  "}\n"
  "\n"
  "\n"
  "PARAMETER {\n"
  "	eh	= -40	(mV) :Bazhenov C++ line 290\n"
  "	:celsius = 36	(degC)\n"
  "	ghbar	= 0.000017 (mho/cm2)    :Bazhenov C++ line 1506\n"
  "	cac	= 0.0015 (mM)		: half-activation of calcium dependence (Bazhenov 2002 neur271.c line 272)\n"
  "	k2	= 0.0004 (1/ms)		: inverse of time constant (Bazhenov 2002 neur271.c line 292)\n"
  "	Pc	= 0.007			: half-activation of CB protein dependence (Bazhenov 2002 neur271.c line 1497)\n"
  "	k4	= 0.001	(1/ms)		: backward binding on Ih (Bazhenov 2002 neur271.c line 274 & 1498)\n"
  "	nca	= 4			: number of binding sites of ca++ (Bazhenov 2002 neur271.c line 293)\n"
  "	nexp	= 1			: number of binding sites on Ih channels (Bazhenov 2002 neur271.c line 293)\n"
  "	ginc	= 2.0			: augmentation of conductance with Ca++ (Bazhenov 2002 neur271.c line 1493)\n"
  "	taum	= 20	(ms)		: min value of tau (Bazhenov 2002 neur271.c line 293)\n"
  "	fac_gh_TC	= 0	(mV)		: shift of Ih voltage-dependence\n"
  "}\n"
  "\n"
  "\n"
  "STATE {\n"
  "	c1	: closed state of channel\n"
  "	o1	: open state\n"
  "	o2	: CB-bound open state\n"
  "	p0	: resting CB\n"
  "	p1	: Ca++-bound CB\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	v	(mV)\n"
  "	cai	(mM)\n"
  "	ih	(mA/cm2)\n"
  "    gh	(mho/cm2)\n"
  "	h_inf\n"
  "	tau_s	(ms)\n"
  "	alpha	(1/ms)\n"
  "	beta	(1/ms)\n"
  "	k1ca	(1/ms)\n"
  "	k3p	(1/ms)\n"
  "	m\n"
  "	tadj\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE ihkin METHOD sparse\n"
  "\n"
  "	m = o1 + ginc * o2\n"
  "\n"
  "	ih = ghbar * m * (v - eh)\n"
  "}\n"
  "\n"
  "KINETIC ihkin {\n"
  ":\n"
  ":  Here k1ca and k3p are recalculated at each call to evaluate_fct\n"
  ":  because Ca or p1 have to be taken at some power and this does\n"
  ":  not work with the KINETIC block.\n"
  ":  So the kinetics is actually equivalent to\n"
  ":	c1 <-> o1\n"
  ":	p0 + nca Cai <-> p1\n"
  ":	o1 + nexp p1 <-> o2\n"
  "\n"
  "	evaluate_fct(v,cai)\n"
  "\n"
  "	~ c1 <-> o1		(alpha,beta)\n"
  "\n"
  "	~ p0 <-> p1		(k1ca,k2)\n"
  "\n"
  "	~ o1 <-> o2		(k3p,k4)\n"
  "\n"
  "	CONSERVE p0 + p1 = 1\n"
  "	CONSERVE c1 + o1 + o2 = 1\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "INITIAL {\n"
  ":\n"
  ":  Experiments of McCormick & Pape were at 36 deg.C\n"
  ":  Q10 is assumed equal to 3\n"
  ":\n"
  "        tadj = 1.0 :3.0 ^ ((celsius-36 (degC) )/10 (degC) ) : currents.h line 379\n"
  "\n"
  "	evaluate_fct(v,cai)\n"
  "\n"
  "	: commented-out intializations are what was used in 2002 Bazhenov C++ code (see 6/29/19 journal entry)\n"
  "	p1 = 0 :1/(1 + (cac/cai) ^ nca)\n"
  "	o1 = 0  :1/(1 + beta/alpha + (p1/Pc)^nexp )\n"
  "	o2 = 0 :((p1/Pc)^nexp) * o1\n"
  "	c1 = 1 :1-o1-o2\n"
  "	p0 = 1 :1-p1\n"
  "	\n"
  "}\n"
  "\n"
  "\n"
  "UNITSOFF\n"
  "PROCEDURE evaluate_fct(v (mV), cai (mM)) {\n"
  "	: Note: we replaced 'shift' with 'fac_gh_TC' and *added* fac_gh_TC (rather than subtracting shift), in order to replicate the approach in currents.cpp\n"
  "	h_inf = 1 / ( 1 + exp((v+75+fac_gh_TC)/5.5) )\n"
  "\n"
  "	tau_s = (taum + 1000 / ( exp((v+71.5+fac_gh_TC)/14.2) + exp(-(v+89+fac_gh_TC)/11.6) ) ) / tadj\n"
  "\n"
  "	alpha = h_inf / tau_s\n"
  "	beta  = ( 1 - h_inf ) / tau_s\n"
  "\n"
  "	k1ca = k2 * (cai/cac)^nca : currents.cpp line 178\n"
  "\n"
  "	k3p = k4 * (p1/Pc)^nexp : currents.cpp line 179\n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  ":\n"
  ":  procedure for evaluating the activation curve of Ih\n"
  ":\n"
  "PROCEDURE activation(v (mV), cai (mM)) { LOCAL cc\n"
  "\n"
  "	evaluate_fct(v,cai)\n"
  "\n"
  "	cc = 1 / (1 + (cac/cai)^nca ) 		: equil conc of CB-protein (Bazhenov 2002 neur271.c line 281)\n"
  "\n"
  "	m = 1 / ( 1 + beta/alpha + (cc/Pc)^nexp ) :Bazhenov 2002 neur271.c line 282\n"
  "\n"
  "	m = ( 1 + ginc * (cc/Pc)^nexp ) * m :Bazhenov 2002 neur271.c line 283\n"
  "}\n"
  "\n"
  "UNITSON\n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
