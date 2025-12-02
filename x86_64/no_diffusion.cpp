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
static constexpr auto number_of_datum_variables = 8;
static constexpr auto number_of_floating_point_variables = 13;
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
 
#define nrn_init _nrn_init__no_voxel
#define _nrn_initial _nrn_initial__no_voxel
#define nrn_cur _nrn_cur__no_voxel
#define _nrn_current _nrn_current__no_voxel
#define nrn_jacob _nrn_jacob__no_voxel
#define nrn_state _nrn_state__no_voxel
#define _net_receive _net_receive__no_voxel 
#define diffusion diffusion__no_voxel 
 
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
#define conc0 _ml->template fpfield<0>(_iml)
#define conc0_columnindex 0
#define dx_pos _ml->template fpfield<1>(_iml)
#define dx_pos_columnindex 1
#define dx_neg _ml->template fpfield<2>(_iml)
#define dx_neg_columnindex 2
#define dy_pos _ml->template fpfield<3>(_iml)
#define dy_pos_columnindex 3
#define dy_neg _ml->template fpfield<4>(_iml)
#define dy_neg_columnindex 4
#define dz_pos _ml->template fpfield<5>(_iml)
#define dz_pos_columnindex 5
#define dz_neg _ml->template fpfield<6>(_iml)
#define dz_neg_columnindex 6
#define lam _ml->template fpfield<7>(_iml)
#define lam_columnindex 7
#define F _ml->template fpfield<8>(_iml)
#define F_columnindex 8
#define conc _ml->template fpfield<9>(_iml)
#define conc_columnindex 9
#define Dconc _ml->template fpfield<10>(_iml)
#define Dconc_columnindex 10
#define v _ml->template fpfield<11>(_iml)
#define v_columnindex 11
#define _g _ml->template fpfield<12>(_iml)
#define _g_columnindex 12
#define _nd_area *_ml->dptr_field<0>(_iml)
#define conc_xp	*_ppvar[2].get<double*>()
#define _p_conc_xp _ppvar[2].literal_value<void*>()
#define conc_xn	*_ppvar[3].get<double*>()
#define _p_conc_xn _ppvar[3].literal_value<void*>()
#define conc_yp	*_ppvar[4].get<double*>()
#define _p_conc_yp _ppvar[4].literal_value<void*>()
#define conc_yn	*_ppvar[5].get<double*>()
#define _p_conc_yn _ppvar[5].literal_value<void*>()
#define conc_zp	*_ppvar[6].get<double*>()
#define _p_conc_zp _ppvar[6].literal_value<void*>()
#define conc_zn	*_ppvar[7].get<double*>()
#define _p_conc_zn _ppvar[7].literal_value<void*>()
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  2;
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
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"conc0", "nM"},
 {"dx_pos", "1/ms"},
 {"dx_neg", "1/ms"},
 {"dy_pos", "1/ms"},
 {"dy_neg", "1/ms"},
 {"dz_pos", "1/ms"},
 {"dz_neg", "1/ms"},
 {"lam", "1/ms"},
 {"F", "nM/ms"},
 {"conc", "nM"},
 {"conc_xp", "nM"},
 {"conc_xn", "nM"},
 {"conc_yp", "nM"},
 {"conc_yn", "nM"},
 {"conc_zp", "nM"},
 {"conc_zn", "nM"},
 {0, 0}
};
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
 
#define _cvode_ieq _ppvar[8].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"no_voxel",
 "conc0",
 "dx_pos",
 "dx_neg",
 "dy_pos",
 "dy_neg",
 "dz_pos",
 "dz_neg",
 "lam",
 "F",
 0,
 0,
 "conc",
 0,
 "conc_xp",
 "conc_xn",
 "conc_yp",
 "conc_yn",
 "conc_zp",
 "conc_zn",
 0};
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     0, /* conc0 */
     0, /* dx_pos */
     0, /* dx_neg */
     0, /* dy_pos */
     0, /* dy_neg */
     0, /* dz_pos */
     0, /* dz_neg */
     0, /* lam */
     0, /* F */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
  if (nrn_point_prop_) {
    _nrn_mechanism_access_alloc_seq(_prop) = _nrn_mechanism_access_alloc_seq(nrn_point_prop_);
    _ppvar = _nrn_mechanism_access_dparam(nrn_point_prop_);
  } else {
   _ppvar = nrn_prop_datum_alloc(_mechtype, 9, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 13);
 	/*initialize range parameters*/
 	conc0 = _parm_default[0]; /* 0 */
 	dx_pos = _parm_default[1]; /* 0 */
 	dx_neg = _parm_default[2]; /* 0 */
 	dy_pos = _parm_default[3]; /* 0 */
 	dy_neg = _parm_default[4]; /* 0 */
 	dz_pos = _parm_default[5]; /* 0 */
 	dz_neg = _parm_default[6]; /* 0 */
 	lam = _parm_default[7]; /* 0 */
 	F = _parm_default[8]; /* 0 */
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 13);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _no_diffusion_reg() {
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
                                       _nrn_mechanism_field<double>{"conc0"} /* 0 */,
                                       _nrn_mechanism_field<double>{"dx_pos"} /* 1 */,
                                       _nrn_mechanism_field<double>{"dx_neg"} /* 2 */,
                                       _nrn_mechanism_field<double>{"dy_pos"} /* 3 */,
                                       _nrn_mechanism_field<double>{"dy_neg"} /* 4 */,
                                       _nrn_mechanism_field<double>{"dz_pos"} /* 5 */,
                                       _nrn_mechanism_field<double>{"dz_neg"} /* 6 */,
                                       _nrn_mechanism_field<double>{"lam"} /* 7 */,
                                       _nrn_mechanism_field<double>{"F"} /* 8 */,
                                       _nrn_mechanism_field<double>{"conc"} /* 9 */,
                                       _nrn_mechanism_field<double>{"Dconc"} /* 10 */,
                                       _nrn_mechanism_field<double>{"v"} /* 11 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 12 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"conc_xp", "pointer"} /* 2 */,
                                       _nrn_mechanism_field<double*>{"conc_xn", "pointer"} /* 3 */,
                                       _nrn_mechanism_field<double*>{"conc_yp", "pointer"} /* 4 */,
                                       _nrn_mechanism_field<double*>{"conc_yn", "pointer"} /* 5 */,
                                       _nrn_mechanism_field<double*>{"conc_zp", "pointer"} /* 6 */,
                                       _nrn_mechanism_field<double*>{"conc_zn", "pointer"} /* 7 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 8 */);
  hoc_register_prop_size(_mechtype, 13, 9);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "pointer");
  hoc_register_dparam_semantics(_mechtype, 3, "pointer");
  hoc_register_dparam_semantics(_mechtype, 4, "pointer");
  hoc_register_dparam_semantics(_mechtype, 5, "pointer");
  hoc_register_dparam_semantics(_mechtype, 6, "pointer");
  hoc_register_dparam_semantics(_mechtype, 7, "pointer");
  hoc_register_dparam_semantics(_mechtype, 8, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 no_voxel /Users/scoot/TC_sleepnet/mod/no_diffusion.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "Voxel of Nitric Oxide Diffusion";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[1], _dlist1[1];
 static int diffusion(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   double _ldCdiff ;
 _ldCdiff = 0.0 ;
   _ldCdiff = _ldCdiff + dx_pos * ( conc_xp - conc ) ;
   _ldCdiff = _ldCdiff + dx_neg * ( conc_xn - conc ) ;
   _ldCdiff = _ldCdiff + dy_pos * ( conc_yp - conc ) ;
   _ldCdiff = _ldCdiff + dy_neg * ( conc_yn - conc ) ;
   _ldCdiff = _ldCdiff + dz_pos * ( conc_zp - conc ) ;
   _ldCdiff = _ldCdiff + dz_neg * ( conc_zn - conc ) ;
   Dconc = _ldCdiff - lam * conc + F ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 double _ldCdiff ;
 _ldCdiff = 0.0 ;
 _ldCdiff = _ldCdiff + dx_pos * ( conc_xp - conc ) ;
 _ldCdiff = _ldCdiff + dx_neg * ( conc_xn - conc ) ;
 _ldCdiff = _ldCdiff + dy_pos * ( conc_yp - conc ) ;
 _ldCdiff = _ldCdiff + dy_neg * ( conc_yn - conc ) ;
 _ldCdiff = _ldCdiff + dz_pos * ( conc_zp - conc ) ;
 _ldCdiff = _ldCdiff + dz_neg * ( conc_zn - conc ) ;
 Dconc = Dconc  / (1. - dt*( ( - ( lam )*( 1.0 ) ) )) ;
  return 0;
}
 /*END CVODE*/
 static int diffusion (_internalthreadargsproto_) { {
   double _ldCdiff ;
 _ldCdiff = 0.0 ;
   _ldCdiff = _ldCdiff + dx_pos * ( conc_xp - conc ) ;
   _ldCdiff = _ldCdiff + dx_neg * ( conc_xn - conc ) ;
   _ldCdiff = _ldCdiff + dy_pos * ( conc_yp - conc ) ;
   _ldCdiff = _ldCdiff + dy_neg * ( conc_yn - conc ) ;
   _ldCdiff = _ldCdiff + dz_pos * ( conc_zp - conc ) ;
   _ldCdiff = _ldCdiff + dz_neg * ( conc_zn - conc ) ;
    conc = conc + (1. - exp(dt*(( - ( lam )*( 1.0 ) ))))*(- ( _ldCdiff + F ) / ( ( - ( lam )*( 1.0 ) ) ) - conc) ;
   }
  return 0;
}
 
static int _ode_count(int _type){ return 1;}
 
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
  for (int _i=0; _i < 1; ++_i) {
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
  conc = conc0;
 {
   conc = conc0 ;
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
 initmodel(_threadargs_);
}
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{
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
 {   diffusion(_threadargs_);
  } {
   }
}}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {conc_columnindex, 0};  _dlist1[0] = {Dconc_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/Users/scoot/TC_sleepnet/mod/no_diffusion.mod";
    const char* nmodl_file_text = 
  "\n"
  "TITLE Voxel of Nitric Oxide Diffusion\n"
  "\n"
  "COMMENT Equations from\n"
  "        Pablo Fernandez Lopez, Patricio Garcia Baez, and Carmen Paz Suarez Araujo. Nitric Oxide Di\n"
  "\n"
  "\n"
  "usion and Multi-compartmental\n"
  "        Systems: Modeling and Implications DOI: 10.1007/978-3-319-26555-1 59\n"
  "\n"
  "        Instituto Universitario de Ciencias y Tecnologias Cibernticas,\n"
  "        Universidad de Las Palmas de Gran Canaria\n"
  "\n"
  "        Adapted to be used in NEURON and NetPyNE by Scott McElroy\n"
  "        SUNY Downstate scott.mcelroy@downstate.edu\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    THREADSAFE\n"
  "    POINT_PROCESS no_voxel\n"
  "    RANGE conc, conc0, dx_pos, dx_neg, dy_pos, dy_neg, dz_pos, dz_neg, lam, F\n"
  "    POINTER conc_xp, conc_xn, conc_yp, conc_yn, conc_zp, conc_zn\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (molar) = (1/liter)\n"
  "    (nM) = (nanomolar)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    conc0  = 0 (nM)\n"
  "    dx_pos = 0 (1/ms)\n"
  "    dx_neg = 0 (1/ms)\n"
  "    dy_pos = 0 (1/ms)\n"
  "    dy_neg = 0 (1/ms)\n"
  "    dz_pos = 0 (1/ms)\n"
  "    dz_neg = 0 (1/ms)\n"
  "    lam    = 0 (1/ms)\n"
  "    F      = 0 (nM/ms)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    conc_xp (nM)\n"
  "    conc_xn (nM)\n"
  "    conc_yp (nM)\n"
  "    conc_yn (nM)\n"
  "    conc_zp (nM)\n"
  "    conc_zn (nM)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    conc (nM)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    conc = conc0\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE diffusion METHOD cnexp\n"
  "}\n"
  "\n"
  "DERIVATIVE diffusion {\n"
  "    LOCAL dCdiff\n"
  "    dCdiff = 0\n"
  "    dCdiff = dCdiff + dx_pos*(conc_xp - conc)\n"
  "    dCdiff = dCdiff + dx_neg*(conc_xn - conc)\n"
  "    dCdiff = dCdiff + dy_pos*(conc_yp - conc)\n"
  "    dCdiff = dCdiff + dy_neg*(conc_yn - conc)\n"
  "    dCdiff = dCdiff + dz_pos*(conc_zp - conc)\n"
  "    dCdiff = dCdiff + dz_neg*(conc_zn - conc)\n"
  "\n"
  "    conc' = dCdiff - lam*conc + F\n"
  "}\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
