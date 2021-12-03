#ifndef COSMOLOGY_H
#define COSMOLOGY_H

# include "cosmology.h"
# include "gauleg.h"
# include <cmath>
# include <string>
# include <iostream>
# include <fstream>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_spline.h>
# include <gsl/gsl_matrix.h>
# include <gsl/gsl_odeiv.h>
# include <gsl/gsl_errno.h>

struct cosmo
{
    double  Om0,Omk,w0,wa,Omb,hval,th,s8,nspec,ximax,cfac;
};

struct gf_par{
    double Omega0,OmegaL,w0,wa;
};

class cosmology;

double dTime(double,void*);
double dChi (double,void*);
double findmvir(double, void*);
double E_sq(gf_par&, double&);
double dE_sqdz(gf_par&,double&);
void getall(gf_par&,double&,double&,double&,double&);
double d2lnE_sqdz2(gf_par&,double&);
int gf_func(double, const double[], double [], void*);
int gf_jac(double, const double[], double*, double [], void*);
double findzmax(double x, void *params);

class cosmology
{
    protected:
    /// Variables
    double Omega0,Omegal,Omegak,w0,wa,Omegab,h,theta,sigma8,ns,d2norm,t0,rho_crit_0,facmtor;
    double xiNLzetamax,cfactor;

#ifndef SWIG
    const static int N9_16=100;
#else
    const int N9_16=100;
#endif
    double x9_16[N9_16],w9_16[N9_16];

#ifndef SWIG
    const static int N0_2p=57;
#else
    const int N0_2p=57;
#endif
    double x0_2p[N0_2p],w0_2p[N0_2p];

#ifndef SWIG
    const static bool verbose=false;
#else
    const bool verbose=false;
#endif

    private:
    /// Some constants
#ifndef SWIG
    const static double kmpspMpctoGyr=977.813952;
    const static double gee=4.2994e-9;
    const static double c=299792.458;
    const static double e=2.71828183;

#else
    const double kmpspMpctoGyr=977.813952;
    const double gee=4.2994e-9;
    const double c=299792.458;
    const double e=2.71828183;

#endif

#ifndef SWIG
    const static int Nsigma=100;
    const static int Npower=1000;
#else
    const int Nsigma=100;
    const int Npower=1000;
#endif

    /// Options for various functions
    int opt_mf,opt_b,opt_c,opt_ps_L;

    /// Eisenstein and Hu 98 Power spectrum variables
    bool bool_initPS;
    double zeq, keq, zd, Rd, Req, sd, ksilk, alphac, betacinv, alphab, betanode, betab;

    /// Numerical interpolation units for comoving distance
    bool bool_init_Chi;
    gsl_interp_accel *Chi_acc;
    gsl_spline *Chi_spline;

    /// Numerical interpolation units for variance
    bool bool_init_varM_TH_spline;
    gsl_interp_accel *varM_TH_num_acc;
    gsl_spline *varM_TH_num_spline;

    /// Numerical interpolation units for growth factor and f3
    bool bool_init_GF;
    gsl_interp_accel *GF_acc;
    gsl_spline *GF_spline;

    /// Numerical interpolation units for linear power spectra
    bool bool_init_PSL0;
    gsl_interp_accel *PSL0_acc;
    gsl_spline *PSL0_spline;
    double PSL0_dlow, PSL0_dhigh;
    double PSL0_xlow,PSL0_ylow,PSL0_xhigh,PSL0_yhigh;
#ifndef SWIG
    const static double kmin=-5.0;
    const static double kmax=8.0;
#else
    const double kmin=-5.0;
    const double kmax=8.0;
#endif

    /// Numerical interpolation for ukofm
#ifndef SWIG
    const static double cmin=-1.0;
    const static double cmax=3.5;
    const static double krsmax= 8.0;
    const static double krsmin=-6.0;
    const static int Nuk=100;
#else
    const double cmin=-1.0;
    const double cmax=3.5;
    const double krsmax= 8.0;
    const double krsmin=-6.0;
    const int Nuk=100;
#endif
    double uk_krs[Nuk];
    double uk_c[Nuk];
    double ukrsc[Nuk][Nuk];
    bool bool_inituk;
    gsl_interp_accel *uk_c_acc;
    gsl_interp_accel *uk_krs_acc;

    /// Private functions
    void initialize();          //Initializations
    //double Lookback(double);    //Lookback time(z) units 1/H0
    //double Time(double);        //Time(z) units 1/H0
    //double Eofz(double);        //Eofz(z) 
    double Omega(double);       //Omega(z) 
    double Omegaw(double);       //Omegaw(z) 
    double Delta_crit(double);  //Delta_crit(z), Bryan and Norman '98
    double growthfactor(double);  //Growth factor dummy

    void init_growthfactor();   //Speed up for growth factor calculations
    void init_Chi();   //Speed up for growth factor calculations
    double initPS_EH();         // Initialize Power spectrum variables

    void fixd2norm();                // Fix the power spectrum normalization
    double Delta2_EH(double,double); // \Delta^2(k) Power spectrum, k should be in units of h Mpc^{-1}
    double Pk_EH(double,double);     // \P(k) Power spectrum, k should be in units of h Mpc^{-1}
        double TCold_EH(double);         // EH 98 Transfer function for CDM k in Mpc^{-1} 
        double Tb_EH(double);            // EH 98 Transfer function for baryons k in Mpc^{-1}
        double T0master_EH(double,double,double);  // EH 98 Transfer function master for CDM k in Mpc^{-1}
        double TF_EH(double);            // EH 98 Density-wtd Transfer function for CDM+Baryons k in Mpc^{-1}
        double var_TH(double,double);    // \sigma^2(R) R in h^{-1} Mpc assuming a top hat filter in real space
        double var_G(double,double);     // \sigma^2(R) R in h^{-1} Mpc assuming a top hat filter in real space
    void init_varM_TH_spline();          // Initialise the numerical calculation of variance
    double varM_TH(double, double);      // Variance of the density field smoothed with a Tophat filter of mass scale M
    double MF_WA(double,double);         // Warren et al . mass function now obsolete
    double MF_ST(double,double);         // Sheth Tormen mass function now obsolete
    double MF_BH(double,double);         // Bhattacharya mass function experimental
    double MF_TI09(double,double);       // Tinker et al. 2009 mass function also obsolete
    double MF_TI09_350(double,double);   // Tinker et al. 2009 SO(350) mass function, obsolete
    double bias_TWZZ(double, double);    // Tinker et al. 2006 bias function, obsolete

    // Tinker et al. 2010 mass function and bias function variables
    double alpTink;                      // Tinker et al. 2010, mass function normalization
    bool init_Tink;                      // Normalization initialise
    void initTinker(double);             // Initialise the normalization
    double MF_TI10(double,double);       // Tinker et al. 2010 mass function
    double bias_TI10(double, double);    // Tinker et al. 2010 bias
    
    double getzcoll(double);            // Redshift of collapse for a halo of mass M
    double getmstar();                  // M* defined such that sigma(M*)=1.686
    double getMvir(double, double);     // Get Mvir from M200 and redshift
    double getc200(double, double);     // Get c200
    double c_MAC(double,double);        // Concentration a'la Maccio 

    //Friends
    friend double dTime(double,void*);
    friend double dChi (double,void*);
    //friend double df3 (double,void*);
    friend double findmvir(double, void*);
    friend double E_sq(gf_par&, double&);
    friend double dE_sqdz(gf_par&,double&);
    friend void getall(gf_par&,double&,double&,double&,double&);
    friend double d2lnE_sqdz2(gf_par&,double&);
    friend int gf_func(double, const double[], double [], void*);
    friend int gf_jac(double, const double[], double*, double [], void*);
    friend class hod;
    friend double findzmax(double x, void *params);

    public:
    
        // Basic cosmology
        cosmology(); //Constructor
        ~cosmology(); //Destructor
        cosmology(double,double,double,double,double,double,double,double,double,double,double); //Constructor
        cosmology(cosmo); //Constructor
	void cosmo_free();

        double Chiofz_num(double);    //Comoving distance h^{-1} Mpc
        double Chiofz(double);    //Comoving distance h^{-1} Mpc
        double Dlofz(double);     //Luminosity distance h^{-1} Mpc
        double Daofz(double);     //Angular diameter distance h^{-1} Mpc
        double Daofzlh(double,double);     //Angular diameter distance h^{-1} Mpc

        double growthfactor_num(double); // Interpolating routine        
        double dlnDdln1pz(double z);

	void set_optmf(int); // Set the mass function option

        // Power spectrum calculations
        double Delta2_L(double,double); // \Delta^2(k) Power spectrum, k should be in units of h Mpc^{-1} wrapper
        double Pk_L(double,double);     // \P(k) Power spectrum, k should be in units of h Mpc^{-1} wrapper

	void init_powerspectra_L(); //Speed up for power spectra calculations
        double Delta2_L_num(double,double);     // Numerical \Delta^2(k) Power spectrum, k should be in units of h Mpc^{-1} wrapper

	// Mass and bias function wrappers
        double nofm(double,double);
        double bias(double,double);
	
	// Variance related functions
	double varM_TH_num(double,double);       
	double varM_TH_num_deriv(double,double);

	// Group catalog related functions
        double Nplus(double, double);
        double getM(double,double);

        //NFW profile related functions
        void modelNFWhalo(double,double,double&,double&,double&,double&,double&); // Radii in physical units
        void modelNFWhalo_com(double,double,double&,double&,double&,double&,double&); // Mvir, Rvir, cvir, R200, c200
        void modelNFWhalo_com(double,double,double&,double&,double&); //Mvir, Rvir, cvir
        double conc(double,double); //Concentration parameter wrapper
        double ukofm(double,double,double); // Fourier transform of NFW profile
        double uskofm(double,double,double,double); // Fourier transform of NFW profile for satellites
	double ukinterp(double,double); // Interpolation routine to find ukofm quickly

	//Temporarily public
	void ukinit();
	void ukinit2();
    double Eofz(double);        //Eofz(z) 

	/// Set new z
	void setnew_z(double);
	double z_glob;
	double gf_glob;

	/// Access to private variables
	double getOmb();
	double geth();
	double getns();
    double get_cfac();
    double set_cfac(double);

    /// SDSS survey specific functions
    double getzmax(double xL);
    double getLmin(double z, double L1);

    double Time(double);        //Time(z) units 1/H0
    double Lookback(double);    //Lookback time(z) units 1/H0

    /// New functionality added renew cosmology
    void renew(cosmo p);
};

/// This is from haloes.cpp
/// Passing cosmology object, Np and z for calculating M from N_{+}
struct np_params
{
    cosmology *cptr;
    double *z;
    double *Np;
};

///Passing cosmology object
struct c_params
{
    cosmology *cptr;
};


///Passing cosmology object and variance
struct coll_params
{
    cosmology *cptr;
    double *sig;
};

///Passing cosmology object, M200 and z
struct mvir_params
{
    cosmology *cptr;
    double *m200;
    double *z;
};

///Passing cosmology object, cvir, Omega(z), Deltacrit(z)
struct c200_params
{
    cosmology *cptr;
    double *cvir;
    double *omegaz;
    double *dcz;
};

/// This is for powerspectrum.cpp
/// Passing cosmology object, R and z for variance
struct cvar_params
{
    cosmology *cptr;
    double *R;
    double *z;
    bool *psinit;
};


/// Passing cosmology object and z for calculating k_{\sigma}
struct ksig_params
{
    cosmology *cptr;
    double *z;
};

///Passing cosmology object
struct z_params
{
    cosmology *cptr;
    double *mag;
};

#endif
