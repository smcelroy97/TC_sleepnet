/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
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
#include <vector>
using std::size_t;
static auto& std_cerr_stream = std::cerr;
static constexpr auto number_of_datum_variables = 2;
static constexpr auto number_of_floating_point_variables = 12;
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
 
#define nrn_init _nrn_init__NMDA_D1
#define _nrn_initial _nrn_initial__NMDA_D1
#define nrn_cur _nrn_cur__NMDA_D1
#define _nrn_current _nrn_current__NMDA_D1
#define nrn_jacob _nrn_jacob__NMDA_D1
#define nrn_state _nrn_state__NMDA_D1
#define _net_receive _net_receive__NMDA_D1 
#define release release__NMDA_D1 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _internalthreadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
#define _internalthreadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gmax _ml->template fpfield<0>(_iml)
#define gmax_columnindex 0
#define Erev _ml->template fpfield<1>(_iml)
#define Erev_columnindex 1
#define i _ml->template fpfield<2>(_iml)
#define i_columnindex 2
#define g _ml->template fpfield<3>(_iml)
#define g_columnindex 3
#define Ron _ml->template fpfield<4>(_iml)
#define Ron_columnindex 4
#define Roff _ml->template fpfield<5>(_iml)
#define Roff_columnindex 5
#define fn _ml->template fpfield<6>(_iml)
#define fn_columnindex 6
#define synon _ml->template fpfield<7>(_iml)
#define synon_columnindex 7
#define DRon _ml->template fpfield<8>(_iml)
#define DRon_columnindex 8
#define DRoff _ml->template fpfield<9>(_iml)
#define DRoff_columnindex 9
#define _g _ml->template fpfield<10>(_iml)
#define _g_columnindex 10
#define _tsav _ml->template fpfield<11>(_iml)
#define _tsav_columnindex 11
#define _nd_area *_ml->dptr_field<0>(_iml)
 static _nrn_mechanism_cache_instance _ml_real{nullptr};
static _nrn_mechanism_cache_range *_ml{&_ml_real};
static size_t _iml{0};
static Datum *_ppvar;
 static int hoc_nrnpointerindex =  -1;
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
#define Alpha Alpha_NMDA_D1
 double Alpha = 1;
#define Beta Beta_NMDA_D1
 double Beta = 0.0067;
#define Cmax Cmax_NMDA_D1
 double Cmax = 0.5;
#define Cdur Cdur_NMDA_D1
 double Cdur = 0.3;
#define Rtau Rtau_NMDA_D1
 double Rtau = 0;
#define Rinf Rinf_NMDA_D1
 double Rinf = 0;
#define Tr Tr_NMDA_D1
 double Tr = 750;
#define Use Use_NMDA_D1
 double Use = 0;
#define deadtime deadtime_NMDA_D1
 double deadtime = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"Cdur_NMDA_D1", "ms"},
 {"deadtime_NMDA_D1", "ms"},
 {"Cmax_NMDA_D1", "mM"},
 {"Alpha_NMDA_D1", "/ms"},
 {"Beta_NMDA_D1", "/ms"},
 {"Tr_NMDA_D1", "ms"},
 {"Rtau_NMDA_D1", "ms"},
 {"gmax", "uS"},
 {"Erev", "mV"},
 {"i", "nA"},
 {"g", "umho"},
 {0, 0}
};
 static double Roff0 = 0;
 static double Ron0 = 0;
 static double delta_t = 0.01;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"Cdur_NMDA_D1", &Cdur_NMDA_D1},
 {"deadtime_NMDA_D1", &deadtime_NMDA_D1},
 {"Cmax_NMDA_D1", &Cmax_NMDA_D1},
 {"Alpha_NMDA_D1", &Alpha_NMDA_D1},
 {"Beta_NMDA_D1", &Beta_NMDA_D1},
 {"Use_NMDA_D1", &Use_NMDA_D1},
 {"Tr_NMDA_D1", &Tr_NMDA_D1},
 {"Rtau_NMDA_D1", &Rtau_NMDA_D1},
 {"Rinf_NMDA_D1", &Rinf_NMDA_D1},
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
_ppvar = _nrn_mechanism_access_dparam(_prop);
 Node * _node = _nrn_mechanism_access_node(_prop);
v = _nrn_mechanism_access_voltage(_node);
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
"NMDA_D1",
 "gmax",
 "Erev",
 0,
 "i",
 "g",
 0,
 "Ron",
 "Roff",
 0,
 0};
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     0.0001, /* gmax */
     0, /* Erev */
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
    assert(_nrn_mechanism_get_num_vars(_prop) == 12);
 	/*initialize range parameters*/
 	gmax = _parm_default[0]; /* 0.0001 */
 	Erev = _parm_default[1]; /* 0 */
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 12);
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

 extern "C" void _nmda_D1_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
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
                                       _nrn_mechanism_field<double>{"Ron"} /* 4 */,
                                       _nrn_mechanism_field<double>{"Roff"} /* 5 */,
                                       _nrn_mechanism_field<double>{"fn"} /* 6 */,
                                       _nrn_mechanism_field<double>{"synon"} /* 7 */,
                                       _nrn_mechanism_field<double>{"DRon"} /* 8 */,
                                       _nrn_mechanism_field<double>{"DRoff"} /* 9 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 10 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 11 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<void*>{"_tqitem", "netsend"} /* 2 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 3 */);
  hoc_register_prop_size(_mechtype, 12, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 5;
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 NMDA_D1 /Users/scoot/TC_sleepnet/mod/nmda_D1.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "nmda_D1.mod";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[2], _dlist1[2];
 static int release(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DRon = ( synon * Rinf - Ron ) / Rtau ;
   DRoff = - Beta * Roff ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DRon = DRon  / (1. - dt*( ( ( ( - 1.0 ) ) ) / Rtau )) ;
 DRoff = DRoff  / (1. - dt*( ( - Beta )*( 1.0 ) )) ;
  return 0;
}
 /*END CVODE*/
 static int release () {_reset=0;
 {
    Ron = Ron + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / Rtau)))*(- ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - Ron) ;
    Roff = Roff + (1. - exp(dt*(( - Beta )*( 1.0 ))))*(- ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - Roff) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{   neuron::legacy::set_globals_from_prop(_pnt->_prop, _ml_real, _ml, _iml);
    _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  if (_tsav > t){ hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = nullptr;}
 {
   if ( _lflag  == 0.0 ) {
     if ( ( t - _args[3] ) > ( Cdur + deadtime ) ) {
       _args[4] = 1.0 - ( 1.0 - _args[4] * ( 1.0 - Use ) ) * exp ( - ( t - _args[3] ) / Tr ) ;
       synon = synon + _args[4] * _args[0] ;
       _args[1] = _args[1] * exp ( - Beta * ( t - _args[2] ) ) ;
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron + _args[1]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / Rtau ) ) )*( - ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron + _args[1] ;
         }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff - _args[1]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - Beta )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - __primary );
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
     synon = synon - _args[4] * _args[0] ;
     _args[1] = _args[4] * _args[0] * Rinf + ( _args[1] - _args[4] * _args[0] * Rinf ) * exp ( - ( t - _args[2] ) / Rtau ) ;
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron - _args[1]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / Rtau ) ) )*( - ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron - _args[1] ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff + _args[1]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - Beta )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff + _args[1] ;
       }
 _args[2] = t ;
     }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
       _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
 _args[1] = 0.0 ;
   _args[2] = 0.0 ;
   _args[3] = - 100.0 * Rtau ;
   _args[4] = 1.0 ;
   }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
      Node* _nd{};
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
     _ode_spec1 ();
 }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 2; ++_i) {
    _pv[_i] = _nrn_mechanism_get_param_handle(_prop, _slist1[_i]);
    _pvdot[_i] = _nrn_mechanism_get_param_handle(_prop, _dlist1[_i]);
    _cvode_abstol(_atollist, _atol, _i);
  }
 }
 
static void _ode_matsol_instance1(_internalthreadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
      Node* _nd{};
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

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  Roff = Roff0;
  Ron = Ron0;
 {
   synon = 0.0 ;
   Rtau = 1.0 / ( ( Alpha * Cmax ) + Beta ) ;
   Rinf = Cmax * Alpha / ( Cmax * Alpha + Beta ) ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
Node *_nd; double _v; int* _ni; int _cntml;
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _tsav = -1e20;
   _v = _vec_v[_ni[_iml]];
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gmax * ( Ron + Roff ) ;
   fn = 1.0 / ( 1.0 + exp ( - ( v + 25.0 ) / 12.5 ) ) ;
   i = g * fn * ( v - Erev ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_rhs = _nt->node_rhs_storage();
auto const _vec_sav_rhs = _nt->node_sav_rhs_storage();
auto const _vec_v = _nt->node_voltage_storage();
Node *_nd; int* _ni; double _rhs, _v; int _cntml;
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 auto const _g_local = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g_local - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
	 _vec_rhs[_ni[_iml]] -= _rhs;
 
}}

static void nrn_jacob(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_d = _nt->node_d_storage();
auto const _vec_sav_d = _nt->node_sav_d_storage();
auto* const _ml = &_lmr;
Node *_nd; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
  _vec_d[_ni[_iml]] += _g;
 
}}

static void nrn_state(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _cntml;
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
   _v = _vec_v[_ni[_iml]];
 v=_v;
{
 { error =  release();
 if(error){
  std_cerr_stream << "at line 63 in file nmda_D1.mod:\nBREAKPOINT { : would be good to get this in terms of gmax\n";
  std_cerr_stream << _ml << ' ' << _iml << '\n';
  abort_run(error);
}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {Ron_columnindex, 0};  _dlist1[0] = {DRon_columnindex, 0};
 _slist1[1] = {Roff_columnindex, 0};  _dlist1[1] = {DRoff_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/Users/scoot/TC_sleepnet/mod/nmda_D1.mod";
    const char* nmodl_file_text = 
  "TITLE nmda_D1.mod\n"
  "\n"
  "COMMENT\n"
  "This model is a revised version of ampa_D1.mod, which added synaptic depression to ampa.mod, \n"
  "as in the model used in Krishnan 2016 (eLife)\n"
  "First-order synaptic dynamics originally proposed in \"An efficient method for computing synaptic\n"
  "conductances based on a kinetic model of receptor binding\" (Destexhe et. al., 1994).\n"
  "This is an updated version of a mod file originally by Alain Destexhe, ModelDB #18198.\n"
  "This updated version is based primarily on Section 10.1.7 from The NEURON Book, \n"
  "revised to include the parameters Cdur and gmax, as well as a \"deadtime\"\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON{\n"
  "	POINT_PROCESS NMDA_D1\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	GLOBAL deadtime, Cdur, Alpha, Beta,  Rinf, Rtau, Use :values of global variables are the same within a mechanism, but not across mechanisms (e.g., AMPA's deadtime may have a different value than AMPA_D1's deadtime)\n"
  "	RANGE g, gmax, Erev\n"
  "}\n"
  "\n"
  "UNITS{\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "	(mM) = (milli/liter)\n"
  "	(uS) = (micromho)\n"
  "}\n"
  "\n"
  "PARAMETER{ :see line 447 of currents.cpp (from Giri Krishnan) for parameter values\n"
  "	gmax   = 0.0001 (uS) :max conductance of *one* synapse (so in BREAKPOINT, g can be greater than this if there are multiple incoming connections)\n"
  "	Cdur   = 0.3  (ms) :transmitter duration (rising phase)\n"
  "	deadtime=1.0 (ms)  : minimum time between release events\n"
  "	Cmax   = 0.5	(mM)		: max transmitter concentration\n"
  "	Alpha  = 1.0  (/ms mM):forward (binding) rate\n"
  "	Beta   = 0.0067 (/ms):backward (dissociation) rate\n"
  "	Erev   = 0    (mV) :equilibrium potential\n"
  "	Use    = 0.0  :determines how quickly synaptic resources are depleted (the larger this value is, the greater the short-term depressive effect)\n"
  "	Tr     = 750 (ms) :time constant for short-term depression (see Krishnan's currents.h)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v    (mV)   : postsynaptic voltage\n"
  "	i    (nA)   : current = g*(v-Erev)\n"
  "	g    (umho) : conductance\n"
  "	fn			: for implementing NMDA receptor post-synaptic voltage dependence\n"
  "	Rtau (ms)   : time constant of channel building\n"
  "	Rinf        :fraction of open channels if xmtr is present \"forever\"\n"
  "	synon       :sum of weights of all synapses in the \"onset\" state (where weight is assumed to be a unitless factor which scales gmax)\n"
  "}\n"
  "\n"
  "STATE { Ron Roff }  :initialized to 0 by default\n"
  ": Ron and Roff are the total conductances of all synapses\n"
  ": that are in the \"onset\" (transmitter pulse ON)\n"
  ": and \"offset\" (transmitter pulse OFF) states, respectively\n"
  ":declared without units, so units are specified in BREAKPOINT block\n"
  "\n"
  "INITIAL {\n"
  "	:Ron and Roff default to being initialized to zero\n"
  "	synon = 0\n"
  "	Rtau = 1 / ((Alpha * Cmax) + Beta)\n"
  "	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)\n"
  "}\n"
  "\n"
  "BREAKPOINT { : would be good to get this in terms of gmax\n"
  "	SOLVE release METHOD cnexp\n"
  "	g = gmax * (Ron + Roff) :max value is gmax*synon*Rinf\n"
  "	fn = 1/(1+exp(-(v + 25 (mV))/12.5 (mV))) :see \"Synaptic Currents\" section of Bazhenov 2002 paper (also in Krishnan 2016 currents.cpp)\n"
  "	i = g*fn*(v - Erev)\n"
  "}\n"
  "\n"
  "DERIVATIVE release {\n"
  "	Ron'  = (synon*Rinf - Ron)/Rtau\n"
  "	Roff' = -Beta*Roff\n"
  "}\n"
  "\n"
  ":weight is assumed to be a unitless factor which scales gmax.\n"
  ":short-term depression achieved using the factor 'E,' which is proportion \n"
  ":of available presynaptic resources. \n"
  "\n"
  "NET_RECEIVE(weight, r0, t0 (ms), lastspike (ms), E) {\n"
  "	INITIAL{\n"
  "		r0 = 0\n"
  "		t0 = 0 (ms) :this value doesn't really matter\n"
  "		lastspike = -100*Rtau :initialize to large neg value so that do not get depression on first spike\n"
  "		E  = 1 :synapse starts out at full strength\n"
  "	}\n"
  "	:flag is an implicit argument of NET_RECEIVE, normally 0\n"
  "	if (flag == 0){ :flag==0 implies a spike is received \n"
  "		:a spike arrived; ignore it if we are already within either a spike state, or deadtime\n"
  "		if( (t-lastspike)>(Cdur + deadtime) ){\n"
  "			:for 'E,' see \"Synaptic Currents\" section of Bazhenov 2002 paper (and around line 500 of Krishnan 2016 currents.cpp)\n"
  "			E = 1 - (1 - E*(1-Use)) * exp(-(t-lastspike)/Tr) :note that we are assuming the last spike in this stream occurred at t0-Cdur\n"
  "			synon = synon + E*weight :weight is scaled by 'E' to implement synaptic depression\n"
  "			r0 = r0*exp(-Beta*(t-t0)) :r0 at start of onset state\n"
  "			Ron = Ron + r0\n"
  "			Roff = Roff - r0\n"
  "			t0 = t :update time of most recent state change\n"
  "			lastspike = t :update most recent spike time\n"
  "			:come again in Cdur with flag = 1\n"
  "			net_send(Cdur, 1)\n"
  "		}\n"
  "	}\n"
  "	if (flag == 1) {\n"
  "		: \"turn off transmitter\"\n"
  "		: i.e. this synapse enters the offset state\n"
  "		synon = synon - E*weight :note we need to include 'E' in order to undo the addition to synon in the above block\n"
  "		: r0 at start of offset state\n"
  "		r0 = E*weight*Rinf + (r0-E*weight*Rinf)*exp(-(t-t0)/Rtau)\n"
  "		Ron = Ron - r0\n"
  "		Roff = Roff + r0\n"
  "		t0 = t :update time of most recent state change\n"
  "	}\n"
  "}\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
