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
static constexpr auto number_of_datum_variables = 2;
static constexpr auto number_of_floating_point_variables = 19;
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
 
#define nrn_init _nrn_init__GABA_B
#define _nrn_initial _nrn_initial__GABA_B
#define nrn_cur _nrn_cur__GABA_B
#define _nrn_current _nrn_current__GABA_B
#define nrn_jacob _nrn_jacob__GABA_B
#define nrn_state _nrn_state__GABA_B
#define _net_receive _net_receive__GABA_B 
#define bindkin bindkin__GABA_B 
 
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
#define gmax _ml->template fpfield<0>(_iml)
#define gmax_columnindex 0
#define Erev _ml->template fpfield<1>(_iml)
#define Erev_columnindex 1
#define i _ml->template fpfield<2>(_iml)
#define i_columnindex 2
#define g _ml->template fpfield<3>(_iml)
#define g_columnindex 3
#define R _ml->template fpfield<4>(_iml)
#define R_columnindex 4
#define Ron _ml->template fpfield<5>(_iml)
#define Ron_columnindex 5
#define Roff _ml->template fpfield<6>(_iml)
#define Roff_columnindex 6
#define G _ml->template fpfield<7>(_iml)
#define G_columnindex 7
#define Gn _ml->template fpfield<8>(_iml)
#define Gn_columnindex 8
#define synon _ml->template fpfield<9>(_iml)
#define synon_columnindex 9
#define Rinf _ml->template fpfield<10>(_iml)
#define Rinf_columnindex 10
#define Rtau _ml->template fpfield<11>(_iml)
#define Rtau_columnindex 11
#define Beta _ml->template fpfield<12>(_iml)
#define Beta_columnindex 12
#define DRon _ml->template fpfield<13>(_iml)
#define DRon_columnindex 13
#define DRoff _ml->template fpfield<14>(_iml)
#define DRoff_columnindex 14
#define DG _ml->template fpfield<15>(_iml)
#define DG_columnindex 15
#define v _ml->template fpfield<16>(_iml)
#define v_columnindex 16
#define _g _ml->template fpfield<17>(_iml)
#define _g_columnindex 17
#define _tsav _ml->template fpfield<18>(_iml)
#define _tsav_columnindex 18
#define _nd_area *_ml->dptr_field<0>(_iml)
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 static void _hoc_setdata(void*);
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {0, 0}
};
 static Member_func _member_func[] = {
 {"loc", _hoc_loc_pnt},
 {"has_loc", _hoc_has_loc},
 {"get_loc", _hoc_get_loc_pnt},
 {0, 0}
};
 /* declare global and static user variables */
 #define gind 0
 #define _gth 0
#define Cdur Cdur_GABA_B
 double Cdur = 0.3;
#define Cmax Cmax_GABA_B
 double Cmax = 0.5;
#define KD KD_GABA_B
 double KD = 100;
#define K4 K4_GABA_B
 double K4 = 0.033;
#define K3 K3_GABA_B
 double K3 = 0.098;
#define K2 K2_GABA_B
 double K2 = 0.0013;
#define K1 K1_GABA_B
 double K1 = 0.52;
#define deadtime deadtime_GABA_B
 double deadtime = 1;
#define n n_GABA_B
 double n = 4;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"Cmax_GABA_B", "mM"},
 {"Cdur_GABA_B", "ms"},
 {"deadtime_GABA_B", "ms"},
 {"K1_GABA_B", "/ms"},
 {"K2_GABA_B", "/ms"},
 {"K3_GABA_B", "/ms"},
 {"K4_GABA_B", "/ms"},
 {"gmax", "umho"},
 {"Erev", "mV"},
 {"i", "nA"},
 {"g", "umho"},
 {0, 0}
};
 static double G0 = 0;
 static double Roff0 = 0;
 static double Ron0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"Cmax_GABA_B", &Cmax_GABA_B},
 {"Cdur_GABA_B", &Cdur_GABA_B},
 {"deadtime_GABA_B", &deadtime_GABA_B},
 {"K1_GABA_B", &K1_GABA_B},
 {"K2_GABA_B", &K2_GABA_B},
 {"K3_GABA_B", &K3_GABA_B},
 {"K4_GABA_B", &K4_GABA_B},
 {"KD_GABA_B", &KD_GABA_B},
 {"n_GABA_B", &n_GABA_B},
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 static void nrn_alloc(Prop*);
static void nrn_init(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_state(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 static void nrn_cur(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_jacob(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(Prop*, int, neuron::container::data_handle<double>*, neuron::container::data_handle<double>*, double*, int);
static void _ode_spec(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void _ode_matsol(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
#define _cvode_ieq _ppvar[3].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GABA_B",
 "gmax",
 "Erev",
 0,
 "i",
 "g",
 "R",
 0,
 "Ron",
 "Roff",
 "G",
 0,
 0};
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     0.0001, /* gmax */
     -95, /* Erev */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
  if (nrn_point_prop_) {
    _nrn_mechanism_access_alloc_seq(_prop) = _nrn_mechanism_access_alloc_seq(nrn_point_prop_);
    _ppvar = _nrn_mechanism_access_dparam(nrn_point_prop_);
  } else {
   _ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 19);
 	/*initialize range parameters*/
 	gmax = _parm_default[0]; /* 0.0001 */
 	Erev = _parm_default[1]; /* -95 */
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 19);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 
#define _tqitem &(_ppvar[2])
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _gaba_B_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 hoc_register_parm_default(_mechtype, &_parm_default);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"gmax"} /* 0 */,
                                       _nrn_mechanism_field<double>{"Erev"} /* 1 */,
                                       _nrn_mechanism_field<double>{"i"} /* 2 */,
                                       _nrn_mechanism_field<double>{"g"} /* 3 */,
                                       _nrn_mechanism_field<double>{"R"} /* 4 */,
                                       _nrn_mechanism_field<double>{"Ron"} /* 5 */,
                                       _nrn_mechanism_field<double>{"Roff"} /* 6 */,
                                       _nrn_mechanism_field<double>{"G"} /* 7 */,
                                       _nrn_mechanism_field<double>{"Gn"} /* 8 */,
                                       _nrn_mechanism_field<double>{"synon"} /* 9 */,
                                       _nrn_mechanism_field<double>{"Rinf"} /* 10 */,
                                       _nrn_mechanism_field<double>{"Rtau"} /* 11 */,
                                       _nrn_mechanism_field<double>{"Beta"} /* 12 */,
                                       _nrn_mechanism_field<double>{"DRon"} /* 13 */,
                                       _nrn_mechanism_field<double>{"DRoff"} /* 14 */,
                                       _nrn_mechanism_field<double>{"DG"} /* 15 */,
                                       _nrn_mechanism_field<double>{"v"} /* 16 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 17 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 18 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<void*>{"_tqitem", "netsend"} /* 2 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 3 */);
  hoc_register_prop_size(_mechtype, 19, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 4;
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GABA_B /Users/scoot/TC_sleepnet/mod/gaba_B.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[3], _dlist1[3];
 static int bindkin(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   DRon = synon * K1 * Cmax - ( K1 * Cmax + K2 ) * Ron ;
   DRoff = - K2 * Roff ;
   R = Ron + Roff ;
   DG = K3 * R - K4 * G ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 DRon = DRon  / (1. - dt*( ( - ( ( K1 * Cmax + K2 ) )*( 1.0 ) ) )) ;
 DRoff = DRoff  / (1. - dt*( ( - K2 )*( 1.0 ) )) ;
 R = Ron + Roff ;
 DG = DG  / (1. - dt*( ( - ( K4 )*( 1.0 ) ) )) ;
  return 0;
}
 /*END CVODE*/
 static int bindkin (_internalthreadargsproto_) { {
    Ron = Ron + (1. - exp(dt*(( - ( ( K1 * Cmax + K2 ) )*( 1.0 ) ))))*(- ( ( ( synon )*( K1 ) )*( Cmax ) ) / ( ( - ( ( ( K1 )*( Cmax ) + K2 ) )*( 1.0 ) ) ) - Ron) ;
    Roff = Roff + (1. - exp(dt*(( - K2 )*( 1.0 ))))*(- ( 0.0 ) / ( ( - K2 )*( 1.0 ) ) - Roff) ;
   R = Ron + Roff ;
    G = G + (1. - exp(dt*(( - ( K4 )*( 1.0 ) ))))*(- ( ( K3 )*( R ) ) / ( ( - ( K4 )*( 1.0 ) ) ) - G) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  Prop* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _nrn_mechanism_cache_instance _ml_real{_pnt->_prop};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
   _thread = nullptr; double* _globals = nullptr; _nt = (NrnThread*)_pnt->_vnt;   _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  if (_tsav > t){ hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = nullptr;}
 {
   if ( _lflag  == 0.0 ) {
     if ( ( t - _args[3] ) > ( Cdur + deadtime ) ) {
       synon = synon + _args[0] ;
       _args[1] = _args[1] * exp ( - Beta * ( t - _args[2] ) ) ;
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron + _args[1]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - ( ( K1 * Cmax + K2 ) )*( 1.0 ) ) ) ) )*( - ( ( ( synon )*( K1 ) )*( Cmax ) ) / ( ( - ( ( ( K1 )*( Cmax ) + K2 ) )*( 1.0 ) ) ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron + _args[1] ;
         }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff - _args[1]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - K2 )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - K2 )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff - _args[1] ;
         }
 _args[2] = t ;
       _args[3] = t ;
       net_send ( _tqitem, _args, _pnt, t +  Cdur , 1.0 ) ;
       }
     }
   if ( _lflag  == 1.0 ) {
     synon = synon - _args[0] ;
     _args[1] = _args[0] * Rinf + ( _args[1] - _args[0] * Rinf ) * exp ( - ( t - _args[2] ) / Rtau ) ;
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron - _args[1]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - ( ( K1 * Cmax + K2 ) )*( 1.0 ) ) ) ) )*( - ( ( ( synon )*( K1 ) )*( Cmax ) ) / ( ( - ( ( ( K1 )*( Cmax ) + K2 ) )*( 1.0 ) ) ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron - _args[1] ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff + _args[1]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - K2 )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - K2 )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff + _args[1] ;
       }
 _args[2] = t ;
     }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
     _nrn_mechanism_cache_instance _ml_real{_pnt->_prop};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
  Datum* _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  Datum* _thread = nullptr;
  double* _globals = nullptr;
  NrnThread* _nt = (NrnThread*)_pnt->_vnt;
 _args[1] = 0.0 ;
   _args[2] = 0.0 ;
   _args[3] = - 10.0 * ( Cdur + deadtime ) ;
   }
 
static int _ode_count(int _type){ return 3;}
 
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
     _ode_spec1 (_threadargs_);
 }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  Datum* _ppvar;
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 3; ++_i) {
    _pv[_i] = _nrn_mechanism_get_param_handle(_prop, _slist1[_i]);
    _pvdot[_i] = _nrn_mechanism_get_param_handle(_prop, _dlist1[_i]);
    _cvode_abstol(_atollist, _atol, _i);
  }
 }
 
static void _ode_matsol_instance1(_internalthreadargsproto_) {
 _ode_matsol1 (_threadargs_);
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
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  G = G0;
  Roff = Roff0;
  Ron = Ron0;
 {
   R = 0.0 ;
   G = 0.0 ;
   synon = 0.0 ;
   Rinf = K1 * Cmax / ( K1 * Cmax + K2 ) ;
   Rtau = 1.0 / ( K1 * Cmax + K2 ) ;
   Beta = K2 ;
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
 _tsav = -1e20;
   _v = _vec_v[_ni[_iml]];
 v = _v;
 initmodel(_threadargs_);
}
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   Gn = G * G * G * G ;
   g = gmax * Gn / ( Gn + KD ) ;
   i = g * ( v - Erev ) ;
   }
 _current += i;

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
 auto const _g_local = _nrn_current(_threadargscomma_ _v + .001);
 	{ _rhs = _nrn_current(_threadargscomma_ _v);
 	}
 _g = (_g_local - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
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
 {   bindkin(_threadargs_);
  }}}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {Ron_columnindex, 0};  _dlist1[0] = {DRon_columnindex, 0};
 _slist1[1] = {Roff_columnindex, 0};  _dlist1[1] = {DRoff_columnindex, 0};
 _slist1[2] = {G_columnindex, 0};  _dlist1[2] = {DG_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/Users/scoot/TC_sleepnet/mod/gaba_B.mod";
    const char* nmodl_file_text = 
  ": $Id: gabab.mod,v 1.9 2004/06/17 16:04:05 billl Exp $\n"
  "\n"
  "COMMENT\n"
  "This is a modified version of the mod file found here:\n"
  "https://senselab.med.yale.edu/modeldb/showmodel.cshtml?model=150538&file=%2fxietal2013%2fgababsyn.mod#tabs-2\n"
  "Modifications by Christian Fink\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "	Kinetic model of GABA-B receptors\n"
  "	=================================\n"
  "\n"
  "  MODEL OF SECOND-ORDER G-PROTEIN TRANSDUCTION AND FAST K+ OPENING\n"
  "  WITH COOPERATIVITY OF G-PROTEIN BINDING TO K+ CHANNEL\n"
  "\n"
  "  PULSE OF TRANSMITTER\n"
  "\n"
  "  SIMPLE KINETICS WITH NO DESENSITIZATION\n"
  "\n"
  "	Features:\n"
  "\n"
  "  	  - peak at 100 ms; time course fit to Tom Otis' PSC\n"
  "	  - SUMMATION (psc is much stronger with bursts)\n"
  "\n"
  "\n"
  "	Approximations:\n"
  "\n"
  "	  - single binding site on receptor	\n"
  "	  - model of alpha G-protein activation (direct) of K+ channel\n"
  "	  - G-protein dynamics is second-order; simplified as follows:\n"
  "		- saturating receptor\n"
  "		- no desensitization\n"
  "		- Michaelis-Menten of receptor for G-protein production\n"
  "		- \"resting\" G-protein is in excess\n"
  "		- Quasi-stat of intermediate enzymatic forms\n"
  "	  - binding on K+ channel is fast\n"
  "\n"
  "\n"
  "	Kinetic Equations:\n"
  "\n"
  "	  dR/dt = K1 * T * (1-R-D) - K2 * R\n"
  "\n"
  "	  dG/dt = K3 * R - K4 * G\n"
  "\n"
  "	  R : activated receptor\n"
  "	  T : transmitter\n"
  "	  G : activated G-protein\n"
  "	  K1,K2,K3,K4 = kinetic rate cst\n"
  "\n"
  "  n activated G-protein bind to a K+ channel:\n"
  "\n"
  "	n G + C <-> O		(Alpha,Beta)\n"
  "\n"
  "  If the binding is fast, the fraction of open channels is given by:\n"
  "\n"
  "	O = G^n / ( G^n + KD )\n"
  "\n"
  "  where KD = Beta / Alpha is the dissociation constant\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "  Parameters estimated from patch clamp recordings of GABAB PSP's in\n"
  "  rat hippocampal slices (Otis et al, J. Physiol. 463: 391-407, 1993).\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "  PULSE MECHANISM\n"
  "\n"
  "  Kinetic synapse with release mechanism as a pulse.  \n"
  "\n"
  "  Warning: for this mechanism to be equivalent to the model with diffusion \n"
  "  of transmitter, small pulses must be used...\n"
  "\n"
  "  For a detailed model of GABAB:\n"
  "\n"
  "  Destexhe, A. and Sejnowski, T.J.  G-protein activation kinetics and\n"
  "  spill-over of GABA may account for differences between inhibitory responses\n"
  "  in the hippocampus and thalamus.  Proc. Natl. Acad. Sci. USA  92:\n"
  "  9515-9519, 1995.\n"
  "\n"
  "  For a review of models of synaptic currents:\n"
  "\n"
  "  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of \n"
  "  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; \n"
  "  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1996.\n"
  "\n"
  "  This simplified model was introduced in:\n"
  "\n"
  "  Destexhe, A., Bal, T., McCormick, D.A. and Sejnowski, T.J.\n"
  "  Ionic mechanisms underlying synchronized oscillations and propagating\n"
  "  waves in a model of ferret thalamic slices. Journal of Neurophysiology\n"
  "  76: 2049-2070, 1996.  \n"
  "\n"
  "  See also http://www.cnl.salk.edu/~alain\n"
  "\n"
  "\n"
  "\n"
  "  Alain Destexhe, Salk Institute and Laval University, 1995\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS GABA_B\n"
  "	RANGE R, G, g, gmax, Erev\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	GLOBAL Cmax, Cdur, deadtime\n"
  "	GLOBAL K1, K2, K3, K4, KD\n"
  "}\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "	(mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    gmax = 0.0001 (umho)\n"
  "	Cmax	= 0.5	(mM)		: max transmitter concentration\n"
  "	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)\n"
  "	deadtime=1.0 (ms) : minimum time between release events\n"
  ":\n"
  ":	From Kfit with long pulse (5ms 0.5mM)\n"
  ":   see CellSyn.h for 'Kx' values\n"
  "	K1	= 0.52	(/ms mM)	: forward binding rate to receptor (currents.h line 665)\n"
  "	K2	= 0.0013 (/ms)		: backward (unbinding) rate of receptor (currents.h line 666)\n"
  "	K3	= 0.098 (/ms)		: rate of G-protein production (currents.h line 667)\n"
  "	K4	= 0.033 (/ms)		: rate of G-protein decay (currents.h line 668)\n"
  "	KD	= 100			: dissociation constant of K+ channel  (currents.cpp line 348)\n"
  "	n	= 4			: nb of binding sites of G-protein on K+ (currents.cpp line 348)\n"
  "	Erev	= -95	(mV)		: reversal potential (E_K) (currents.cpp line 346)\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)		: postsynaptic voltage\n"
  "	i 		(nA)		: current = g*(v - Erev)\n"
  "	g 		(umho)		: conductance\n"
  "	Gn\n"
  "	R				: fraction of activated receptor\n"
  "	synon\n"
  "	Rinf\n"
  "	Rtau (ms)\n"
  "	Beta (/ms)\n"
  "}\n"
  "\n"
  "\n"
  "STATE {\n"
  "	Ron Roff\n"
  "	G				: fraction of activated G-protein\n"
  "}\n"
  "\n"
  "\n"
  "INITIAL {\n"
  "	R = 0\n"
  "	G = 0\n"
  "	synon = 0\n"
  "	Rinf = K1*Cmax/(K1*Cmax + K2)\n"
  "	Rtau = 1/(K1*Cmax + K2)\n"
  "	Beta = K2\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE bindkin METHOD cnexp\n"
  "	Gn = G*G*G*G : ^n = 4\n"
  "	g = gmax * Gn / (Gn+KD)\n"
  "	i = g*(v - Erev)\n"
  "}\n"
  "\n"
  "\n"
  "DERIVATIVE bindkin {\n"
  "	Ron' = synon*K1*Cmax - (K1*Cmax + K2)*Ron\n"
  "	Roff' = -K2*Roff\n"
  "	R = Ron + Roff\n"
  "	G' = K3 * R - K4 * G\n"
  "}\n"
  "\n"
  ": following supports both saturation from single input and\n"
  ": summation from multiple inputs\n"
  ": Note: automatic initialization of all reference args to 0\n"
  ": except first\n"
  "\n"
  ":Note: we replaced the NET_RECEIVE block of the original file\n"
  ":with modifications simliar to those from Section 10.1.7 of The NEURON Book. \n"
  ":Differences from 10.1.7 are that a deadtime is included, and presynaptic spikes that occur\n"
  ":during an 'on' state or deadtime are simply ignored (though for Cdur of 0.3ms and deadtime of 1.0 ms,\n"
  ":this should never happen)\n"
  ":also, netcon weight is assumed to be 1.0; changes in synaptic strength should be implemented by changing \n"
  ":gmax, since the nonlinearity in the equations results in a nonlinear relationship between 'weight' and 'gmax' \n"
  "\n"
  "NET_RECEIVE(weight, r0, t0 (ms), lastspike (ms)) {\n"
  "	INITIAL{\n"
  "		r0 = 0\n"
  "		t0 = 0 (ms)\n"
  "		lastspike = -10*(Cdur + deadtime) :this statement must be here, and not in other INITIAL block, in order to prevent segmentation fault\n"
  "	}\n"
  "	:flag is an implicit argument of NET_RECEIVE, normally 0\n"
  "	if (flag == 0){ :flag==0 implies a spike is received\n"
  "		:a spike arrived; ignore it if we are already within either a spike state, or deadtime\n"
  "		if( (t-lastspike)>(Cdur + deadtime) ){\n"
  "			synon = synon + weight\n"
  "			r0 = r0*exp(-Beta*(t-t0)) :r0 at start of onset state\n"
  "			Ron = Ron + r0\n"
  "			Roff = Roff - r0\n"
  "			t0 = t\n"
  "			lastspike = t\n"
  "			:come again in Cdur with flag = 1\n"
  "			net_send(Cdur, 1)\n"
  "		}\n"
  "	}\n"
  "	if (flag == 1) {\n"
  "		: \"turn off transmitter\"\n"
  "		: i.e. this synapse enters the offset state\n"
  "		synon = synon - weight\n"
  "		: r0 at start of offset state\n"
  "		r0 = weight*Rinf + (r0-weight*Rinf)*exp(-(t-t0)/Rtau)\n"
  "		Ron = Ron - r0\n"
  "		Roff = Roff + r0\n"
  "		t0 = t\n"
  "	}\n"
  "\n"
  "COMMENT\n"
  "	:Note: this is the original NET_RECEIVE block from https://senselab.med.yale.edu/modeldb/showmodel.cshtml?model=150538&file=%2fxietal2013%2fgababsyn.mod#tabs-2\n"
  "	if (flag == 1) { : at end of Cdur pulse so turn off\n"
  "		r0 = weight*(Rinf + (r0 - Rinf)*exp(-(t - t0)/Rtau))\n"
  "		t0 = t\n"
  "		synon = synon - weight\n"
  "		Ron = Ron-r0\n"
  "		Roff = Roff+r0\n"
  "        }else{ : at beginning of Cdur pulse so turn on\n"
  "		r0 = weight*r0*exp(-Beta*(t - t0)) :CF: the factor of weight here does not make sense; it should only be applied to the exponential rise\n"
  "		t0 = t\n"
  "		synon = synon + weight C\n"
  "		Ron = Ron+r0\n"
  "		Roff = Roff-r0\n"
  "		:come again in Cdur\n"
  "		net_send(Cdur, 1)\n"
  "        }\n"
  "ENDCOMMENT\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
