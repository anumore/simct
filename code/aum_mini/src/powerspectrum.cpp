/// Power spectrum routines part of cosmology class
#include "cosmology.h"

/// Linear Powerspectrum wrappers
double cosmology::Delta2_L(double k,double z)
{
    double result;

    if(opt_ps_L==1)
    {
        result=Delta2_EH(k, z);
    } else
    {
        std::cout<<"Linear power spectrum option not supported yet."<<std::endl;
        exit(0);
    }

    return result;

}

double cosmology::Pk_L(double k,double z)
{
    double result;

    if(opt_ps_L==1)
    {
        result=Pk_EH(k, z);
    } else
    {
        std::cout<<"Linear power spectrum option not supported yet."<<std::endl;
        exit(0);
    }

    return result;

}

/// Initialize power spectrum variables Eisenstein and Hu 1998
double cosmology::initPS_EH()
{

    double O0h2=Omega0*h*h;
    double Obh2=Omegab*h*h;

    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<" O0h2="<<O0h2<<" "<<" Obh2="<<Obh2<<std::endl;
    }

    /// Equations refer to Eisenstein and Hu '98.
    // Eq. 2
    zeq = 2.50*10000*O0h2*pow(theta,-4.);
    //std::cout<<"# "<<" zeq="<<zeq<<" "<<O0h2<<" "<<theta<<std::endl;
    //std::cout<<"# "<<" zeq="<<zeq<<std::endl;

    // Eq. 3
    //keq = sqrt(2*O0h2*pow(100./c,2.)*zeq);
    keq = 0.0746*Omega0*h*h*pow(theta,-2.);
    //std::cout<<"# "<<" keq="<<keq<<std::endl;//" "<<7.46e-2*Omega0*h*h*pow(theta,-2.)<<std::endl;

    //Eq. 4, define this to be 1+zdrag, needs a fix in all releases.
    double b1=0.313*pow(O0h2,-0.419)*(1.+0.607*pow(O0h2,0.674));
    double b2=0.238*pow(O0h2,0.223);
    zd = 1+1291*pow(O0h2,0.251)*(1.+b1*pow(Obh2,b2))/(1.+0.659*pow(O0h2,0.828));
    //std::cout<<"# "<<" zd="<<zd<<std::endl;

    //Eq. 5
    Rd = 31.5*Obh2*pow(theta,-4.)*(1000./zd);
    Req = 31.5*Obh2*pow(theta,-4.)*(1000./zeq);
    //std::cout<<"# "<<" Rd="<<Rd<<std::endl;
    //std::cout<<"# "<<" Req="<<Req<<std::endl;

    //Eq.6
    sd = 2./(3.*keq)*sqrt(6./Req)*log( ( sqrt(1.+Rd) + sqrt(Rd+Req) )/(1.+sqrt(Req))  );
    //std::cout<<"# "<<" sd="<<sd<<std::endl;

    //Eq. 7
    //double factor=pow(10.4*O0h2,-0.95);
    ksilk = 1.6*pow(Obh2,0.52)*pow(O0h2,0.73)*(1.+pow(10.4*O0h2,-0.95));
    //std::cout<<"# "<<" ksilk="<<ksilk<<std::endl;

    //Eq. 11
    double a1 = pow(46.9*O0h2,0.670)*(1.+pow(32.1*O0h2,-0.532));
    double a2 = pow(12.0*O0h2,0.424)*(1.+pow(45.0*O0h2,-0.582));
    alphac = pow(a1,-Omegab/Omega0)*pow(a2,-pow(Omegab/Omega0,3.));
    //std::cout<<"# "<<" alphac="<<alphac<<std::endl;

    //Eq. 12
    double Omegac=Omega0-Omegab;
    double bb1 = 0.944/(1.+pow(458*O0h2,-0.708));
    double bb2 = pow(0.395*O0h2,-0.0266);
    betacinv=1.+bb1*(pow(Omegac/Omega0,bb2)-1.);
    //std::cout<<"# "<<" betac="<<1./betacinv<<std::endl;

    //Eq. 14
    double y = zeq/zd;
    double Gy= y*( -6.*sqrt(1.+y) + (2.+3.*y)*log( (sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)  )  );
    alphab=2.07*keq*sd*pow(1.+Rd,-0.75)*Gy;
    //std::cout<<"# "<<" alphab="<<alphab<<" "<<Gy<<" "<<log( (sqrt(1.+y)+1.)/(sqrt(1.+y)-1.))<<std::endl;
    //std::cout<<"# "<<-6.*sqrt(1.+y)<<" "<<(2.+3.*y)*log( (sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)  )<<std::endl;

    //Eq. 23
    betanode=8.41*pow(O0h2,0.435);
    //std::cout<<"# "<<" betanode="<<betanode<<std::endl;

    //Eq. 24
    betab = 0.5 + Omegab/Omega0 + (3. - 2.*Omegab/Omega0)*sqrt(1.+pow(17.2*O0h2,2.));
    //std::cout<<"# "<<" betab="<<betab<<std::endl;

    bool_initPS=true;

    if(verbose){
    std::cout<<"# "<<"Eisenstein and Hu Power spectrum variables initialised."<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }

    return 1;

}

/// T0master for Eisenstein and Hu power spectrum
double cosmology::T0master_EH(double q, double alpc, double betc)
{
    // Eq. 20
    // The value q = k/(keq*13.41) is passed instead of k to avoid repititive
    // calculation;
    double CT0 = 14.2/alpc+386./(1+69.9*pow(q,1.08));

    //std::cout<<" Check master "<<CT0<<" "<<alpc<<" "<<betc<<" "<<q<<std::endl;
    
    // Eq. 21 This is also called TF_pressureless
    return log( e + 1.8*betc*q )/( log(e+1.8*betc*q)+ CT0*q*q);

}

/// Transfer function a'la Eisenstein and Hu, k should be supplied in Mpc^{-1}
double cosmology::TF_EH(double kpassed)
{

    //double k = h*kpassed;
    double k = kpassed;
    if (!bool_initPS) 
    {
	if(verbose){
        std::cout<<"# "<<" Power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = initPS_EH(); 

	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }

    double q=k/(keq*13.41);
    // Eq. 21 first part
    double befj0= T0master_EH(q,1.,1.)/(1.+pow(k*sd/5.2,2.)) + alphab/(1.+pow(betab/(k*sd),3.))*exp(-pow(k/ksilk,1.4) );
    // Eq. 22
    double stilde=sd/pow((1.+pow(betanode/(k*sd),3.)),1./3.);
    // Eq. 21 remaining part
    double Tb=befj0*gsl_sf_bessel_j0(k*stilde);

    // Eq. 18
    double f= 1./(1.+pow(k*sd/5.4,4.));
    // Eq. 17
    double TC=f*T0master_EH(q,1.,1./betacinv)+(1-f)*T0master_EH(q,alphac,1./betacinv);
    double TF=Omegab/Omega0*Tb + (Omega0-Omegab)/Omega0*TC;
    return TF;

    //return Omegab/Omega0*Tb_EH(k) + (Omega0-Omegab)/Omega0*TCold_EH(k);
}

/// Delta^2(k) a'la Eisenstein and Hu, k should be in Mpc^{-1}
/// kpassed should be supplied as h Mpc^{-1}.
/// Result is dimensionless as it should be
double cosmology::Delta2_EH(double kpassed, double z)
{
    double k = h*kpassed;

    if (!bool_initPS) 
    {
	if(verbose){
        std::cout<<"# "<<" Power spectrum varibles were not initialised. Now initialising. "<<std::endl;
	}
        double init = initPS_EH(); 
	if(verbose){
        std::cout<<"# "<<init<<std::endl;
	}
    }
    //Normalize the scale to 100./c; If d2norm is 1, then the absolute value as
    //such has no meaning. d2norm should be fixed using sigma8. Check routine,
    //fix_d2norm()
    //if(d2norm==1.){std::cout<<"Did you normalise the power spectrum?";}
    double res=d2norm*pow(c*kpassed/(100.*h),3.+ns)*pow(TF_EH(k),2.);
    if(z!=z_glob){
        res=res*pow(growthfactor_num(z),2.);
    }else{
        res=res*pow(gf_glob,2.);
    }

    return res;
}

/// P(k) a'la Eisenstein and Hu, k should be supplied in Mpc^{-1}
/// kpassed should be supplied as h Mpc^{-1}.
/// Result is in the units of (h^{-1} Mpc)^3
double cosmology::Pk_EH(double kpassed, double z)
{
    double k = kpassed;

    double res=Delta2_EH(k,z);

    return res*(2*M_PI*M_PI)/(pow(k,3.));
}

/// Integrand for var_TH: 1/k \Delta^2(k) (3.*j1(kR)/(kR))^2  k in hMpc^{-1}
double dvar_TH(double x, void * params)
{
    cvar_params c1 = *(cvar_params *) params;
//    std::cout<<1<<std::endl;
    cosmology *c2;
    double *R;
    double *z;
    bool *psi;
    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    psi=c1.psinit;
    double arg=(*R)*x;
    double d2;
    if(*psi){
	d2=(*c2).Delta2_L_num(x,*z);
    }else{
	d2=(*c2).Delta2_L(x,*z);
    }
    double bes= pow(3.*gsl_sf_bessel_j1(arg)/(arg) ,2.);
    //std::cout<<1./x*d2*bes<<"Is this quick"<<std::endl;
    return 1./x*d2*bes;
}

/// Calculate the variance of the power spectrum on a scale R in h^{-1} Mpc
double cosmology::var_TH(double R, double z)
{
    ///Int_{0}^{\infty} dk/k Delta^2(k) (3j1(kR)/kR)^2
    double result1, result2, error;

    gsl_function F;
    F.function = &(dvar_TH);

    //Initialize parameter block
    cvar_params p;
    p.cptr = this;
    p.R=&R;
    p.z=&z;
    p.psinit=&bool_init_PSL0;

    F.params = &p;

    // Do it in two steps as suggested by Frank. First upto 1/Rf and then from
    // 1/Rf to infty.
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F, 0., 1./R, 0, 1e-6, 1000, w, &result1, &error); 
    gsl_integration_workspace_free (w);
    
    gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagiu (&F, 1./R, 0, 1e-6, 1000, v, &result2, &error); 
    gsl_integration_workspace_free (v);
    
    //gsl_integration_qagi (&F, 0, 1e-4, 1000, w, &result, &error); 

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);


    return result1+result2;
    
}

/// Calculate the variance of the power spectrum on a mass scale M in h^{-1} Msun
double cosmology::varM_TH(double M, double z)
{
    ///Calculate R of top hat for M
    double r=pow(M*facmtor,1./3.);
    //std::cout<<M<<" "<<r<<std::endl;
    return var_TH(r,z);
    
}

/// Initialize the spline for numerically calculating the variance
/// This is valid for redshift z=0.0. Should be appropriately multiplied by the
/// growth factor to extend to redshift z.
void cosmology::init_varM_TH_spline()
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for numerical calculation of variance was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance
    double mmin=5.0;
    double mmax=18.0;

    double xx[Nsigma],yy[Nsigma];
    for (int i=0;i<Nsigma;i++)
    {
        double m=mmin+i/(Nsigma-1.)*(mmax-mmin);
        xx[i]=m;
        m=pow(10.,m);
        yy[i]=log(varM_TH(m,0.0))/log(10.);
        //std::cout<<xx[i]<<" "<<yy[i]<<std::endl;
    }
    varM_TH_num_acc = gsl_interp_accel_alloc ();
    varM_TH_num_spline = gsl_spline_alloc (gsl_interp_cspline, Nsigma);
    gsl_spline_init (varM_TH_num_spline,xx,yy,Nsigma);

    bool_init_varM_TH_spline=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of variance is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

/// Calculate the variance of the power spectrum on a mass scale M in h^{-1} Msun
/// numerically by interpolating
double cosmology::varM_TH_num(double M, double z)
{
    //
    if(!bool_init_varM_TH_spline)
    {
        //std::cout<<"# I called it: varM_TH_spline"<<std::endl;
        init_varM_TH_spline();
    }
    double logm=log(M)/log(10.);
    double result;
    if(logm>5.0&&logm<18.0)
    {
        result=gsl_spline_eval (varM_TH_num_spline,logm, varM_TH_num_acc);
    } else
    {
        std::cout<<"# "<<"varM_TH_num_spline: Interpolation not possible M not within range "<<M<<std::endl;
        //exit(0);
        result=log10(varM_TH(M,0.0));
    }
    result = pow(10.,result);
    if(z!=z_glob){
        result=result*pow(growthfactor_num(z),2.);
    }else{
        result=result*pow(gf_glob,2.);
    }
    return result;

}

/// Calculate dlog\sigma/dlogM for a mass M in h^{-1} Msun
/// numerically using the spline
double cosmology::varM_TH_num_deriv(double M, double z)
{
    //
    if(!bool_init_varM_TH_spline)
    {
        //std::cout<<"# I called it: varM_TH_num_deriv"<<std::endl;
        init_varM_TH_spline();
    }
    double logm=log(M)/log(10.);
    double result;
    if(logm>5&&logm<18.0)
    {
        result=gsl_spline_eval_deriv (varM_TH_num_spline,logm, varM_TH_num_acc);
    } else
    {
        std::cout<<"# "<<"varM_TH_num_spline: Interpolation for derivative not possible, M not within range "<<M<<std::endl;
        exit(0);
        result=varM_TH(M,z);
    }
    result = pow(10.,result);
    if(z!=z_glob){
        result=result*pow(growthfactor_num(z),2.);
    }else{
        result=result*pow(gf_glob,2.);
    }
    return result;

}

/// Integrand for var_G: 1/k \Delta^2(k) exp(-(kR)^2)  k in hMpc^{-1}
double dvar_G(double x, void * params)
{
    cvar_params c1 = *(cvar_params *) params;
//    std::cout<<1<<std::endl;
    cosmology *c2;
    double *R;
    double *z;
    c2=c1.cptr;
    R=c1.R;
    z=c1.z;
    double d2=(*c2).Delta2_L(x,*z);
    double res=1./x*d2*exp(-pow(x*(*R),2));
    return res;
}

/// Calculate the variance of the power spectrum on a scale R in h^{-1} Mpc
/// using a Gaussian filter. 
double cosmology::var_G(double R, double z)
{
    //Smith et al. 2003 Eq. 54
    //Int_0^{\infty} dk/k Delta^2(k) exp(-(kR)^2)
    double result1, result2, error;

    gsl_function F;
    F.function = &(dvar_G);

    //Initialize parameter block
    cvar_params p;
    p.cptr = this;
    p.R=&R;
    p.z=&z;

    F.params = &p;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qags (&F, 0., 1./R, 0, 1e-6, 1000, w, &result1, &error); 
    gsl_integration_workspace_free (w);

    gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);
    gsl_integration_qagiu (&F, 1./R, 0, 1e-6, 1000, v, &result2, &error); 
    gsl_integration_workspace_free (v);

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);


    return result1+result2;
    
}

/// Fix the normalization of the power spectrum using the value of \f$\sigma_8\$f supplied by the user.
void cosmology::fixd2norm()
{
    double sig82=var_TH(8.,0.0);
    d2norm=pow(sigma8,2)/sig82;
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The power spectrum has been normalized. Checking:";
    }
    if(fabs(var_TH(8.,0.0)-pow(sigma8,2))>1.0e-4)
    {
        std::cout<<"# "<<"Problem in fixing the power spectrum normalization. Quitting"<<sigma8<<" "<<var_TH(8.,0.0)<<" "<<bool_initPS<<std::endl;
        exit(0);
    }
    if(verbose){
    std::cout<<"# "<<"Check successful."<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
}

/// Power spectra numerical implementation. Initializing spline. 
void cosmology::init_powerspectra_L()
{
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"The spline for linear power spectra was not initialized. Initializing..."<<std::endl;
    }
    //Define the limits for which to calculate the variance

    double xx[Npower],yy[Npower];
    for (int i=0;i<Npower;i++)
    {
        double k=kmin+i/(Npower-1.)*(kmax-kmin);
        xx[i]=k;
	k=pow(10.0,k);
	/// Improvement possible here
        yy[i]=log10(Delta2_L(k,0.0));
    }
    PSL0_acc = gsl_interp_accel_alloc ();
    PSL0_spline = gsl_spline_alloc (gsl_interp_cspline, Npower);
    gsl_spline_init (PSL0_spline,xx,yy,Npower);
    PSL0_dlow=(yy[2]-yy[0])/(xx[2]-xx[0]);
    PSL0_xlow=xx[1]; PSL0_ylow=yy[1];
    PSL0_dhigh=(yy[Npower-3]-yy[Npower-1])/(xx[Npower-3]-xx[Npower-1]);
    PSL0_xhigh=xx[Npower-2]; PSL0_yhigh=yy[Npower-2];

    bool_init_PSL0=true;
    if(verbose){
    std::cout<<"# "<<"The spline for numerical calculation of linear power spectra is now initialized"<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
       
}

/// Calculate Delta^2L numerically with k in h Mpc^{-1}
double cosmology::Delta2_L_num(double k, double z)
{
    //
    if(!bool_init_PSL0)
    {
        init_powerspectra_L();
    }
    double logk=log10(k);
    double result;
    if(logk>kmin&&logk<kmax)
    {
        result=gsl_spline_eval (PSL0_spline,logk, PSL0_acc);
    } else if (logk<=kmin)
    {
        //std::cout<<"# "<<"PSL_spline: Interpolation not possible k not within range "<<k<<std::endl;
        //exit(0);
        //result=log10(Delta2_L(k,z));
        //result=1.0e-30;
	//Extrapolate now
	result=PSL0_ylow+PSL0_dlow*(logk-PSL0_xlow);
    } else if (logk>=kmax)
    {
	result=PSL0_yhigh+PSL0_dhigh*(logk-PSL0_xhigh);
    }

    result = pow(10.,result);
    if(z!=z_glob){
        result=result*pow(growthfactor_num(z),2.);
    }else{
        result=result*pow(gf_glob,2.);
    }
    return result;

}

/// Set new redshift
/// Linear power spectrum is not required to be reinited
void cosmology::setnew_z(double z)
{
    init_Tink=false;
    z_glob=z;
    gf_glob=growthfactor_num(z);
    if(verbose){
        std::cout<<"# NEW REDSHIFT SET"<<std::endl;
    }

}
