#include<stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286 
#define G 4*PI*PI
#define M 1
#define DAY 0.00273972602
#define HOUR 0.00011415525

double dj1a[3];
double de1a[3];
double L1,L2,koct,kphi0;
double mu1;
double a1,a2;
double m0,m1,m2;
double hk1,hk2;
double I0,I1;
double R0,R1;
double mu;
double tf0,tf1;
double tv0,tv1;
double klove0,klove1;
double k0,k1;
double tf;
double Q0,Q1;
double c=63197.639580028;

double dotp(double *v1,double *v2)
{
	int i;
	double sum=0;
	for(i=0;i<3;i++)
	{
		sum+=(v1[i]*v2[i]);
	}
	return sum;
}

int divideVec(double *v1,double scalar,double *v1bs)
{
	v1bs[0]=v1[0]/scalar;
	v1bs[1]=v1[1]/scalar;
	v1bs[2]=v1[2]/scalar;
	return 1;
}

int crossp(double *v1,double *v2, double *cp)
{
	cp[0]=(v1[1]*v2[2])-(v1[2]*v2[1]);
	cp[1]=(v1[2]*v2[0])-(v1[0]*v2[2]);
	cp[2]=(v1[0]*v2[1])-(v1[1]*v2[0]);
	return 1;
}

int updatePreFactors(double *e1v,double *h1v)
{	
	double e1=sqrt(dotp(e1v,e1v));
	double h1=sqrt(dotp(h1v,h1v));

	a1=pow(h1*(m0+m1)/(m0*m1),2)/(G*(m0+m1)*(1-e1*e1));	

	L1=m0*m1*sqrt(G*(m0+m1)*a1)/(m0+m1);
	L2=(m0+m1)*m2*sqrt(G*(m0+m1+m2)*a2)/(m0+m1+m2);
	mu1=m0*m1/(m0+m1);
	kphi0=G*m2*a1*a1/(a2*a2*a2);
	//koct=((m0-m1)/(m0+m1))*(a1/a2);
	koct=0;
	//printf("a1: %e, mu1:%e\n",a1,mu1);
	//tf0=(tv0/9)*pow(a1/R0,8)*((m0*m0)/((m0+m1)*m1))*pow(1+2*k0,-2);
	//tf1=(tv1/9)*pow(a1/R1,8)*((m1*m1)/((m1+m0)*m0))*pow(1+2*k1,-2);

	tf0 = (1./6)*(Q0/klove0)*(m0/m1)*sqrt(1/(G*(m0+m1)))*pow(1./R0,5)*pow(a1,6.5);
	tf1 = (1./6)*(Q1/klove1)*(m1/m0)*sqrt(1/(G*(m1+m0)))*pow(1./R1,5)*pow(a1,6.5);

	//printf("Q0,k0,m0,m1,sqrt,R0,a1 : %e,%e,%e,%e,%e,%e,%e\n",Q0,k0,m0,m1,sqrt(1/(G*(m0+m1))),R0,a1);
	
	//printf("Q0:%e\n",Q0);
	//printf("k0:%e\n",k0);
	//printf("m0:%e\n",m0);
	//printf("m1:%e\n",m1);
	//printf("R0:%e\n",R0);
	//printf("a1:%e\n",a1);

	//printf("tf0 : %e\n",tf0);
	hk1=m0*m1*sqrt(G*(m0+m1)*a1)/(m0+m1);
	hk2=(m0+m1)*m2*sqrt(G*(m0+m1+m2)*a2)/(m0+m1+m2);
}

double dedh_td(double *dedt, double *dhdt,double *dsOm0dt, double *dsOm1dt, double *e1vc, double *h1vc,double *sOm0v,double *sOm1v)
{
	double e1=sqrt(dotp(e1vc,e1vc));
	double h1=sqrt(dotp(h1vc,h1vc));
	double ehat[3];
	double hhat[3];
	double qhat[3];
	double sOm0vh,sOm0ve,sOm0vq;
	double sOm1vh,sOm1ve,sOm1vq;
	double l1d;
	
	divideVec(e1vc,e1,ehat);
	divideVec(h1vc,h1,hhat);
	crossp(hhat,ehat,qhat);
	
	sOm0vh=dotp(sOm0v,hhat);
	sOm0ve=dotp(sOm0v,ehat);
	sOm0vq=dotp(sOm0v,qhat);

	sOm1vh=dotp(sOm1v,hhat);
	sOm1ve=dotp(sOm1v,ehat);
	sOm1vq=dotp(sOm1v,qhat);

	//sOm0ve=10;sOm0vq=10;sOm0vh=10;
	//printf("sOmh:%e ,sOme:%e ,sOmq:%e \n",sOm0vh,sOm0ve,sOm0vq);
	l1d=sqrt(G*(m0+m1)/(a1*a1*a1));
	//printf("l1d: %e\n",l1d);
		
	double V0,W0,X0,Y0,Z0;
	double V1,W1,X1,Y1,Z1;
	V1=0;W1=0;X1=0;Y1=0;Z1=0;

	V0=(9/tf0)*((1+(15.0*pow(e1,2)/4.0)+(15.0*pow(e1,4)/8.0)+(5.0*pow(e1,6)/64))/pow(1-(e1*e1),6.5));
	V0+=-(9/tf0)*(11*sOm0vh/(18*l1d))*((1+(3*pow(e1,2)/2)+(1*pow(e1,4)/8))/pow(1-(e1*e1),5));
	
	V1=(9/tf1)*((1+(15.0*pow(e1,2)/4.0)+(15.0*pow(e1,4)/8.0)+(5*pow(e1,6)/64))/pow(1-(e1*e1),6.5));
	V1+=-(9/tf1)*(11*sOm1vh/(18*l1d))*((1+(3*pow(e1,2)/2)+(1*pow(e1,4)/8))/pow(1-(e1*e1),5));

	W0=(1/tf0)*((1+(15.0*pow(e1,2)/2.0)+(45.0*pow(e1,4)/8.0)+(5*pow(e1,6)/16))/pow(1-(e1*e1),6.5));
	W0+=-(1/tf0)*(sOm0vh/(l1d))*((1+(3*pow(e1,2))+(3*pow(e1,4)/8))/pow(1-(e1*e1),5));

	W1=(1/tf1)*((1+(15.0*pow(e1,2)/2.0)+(45.0*pow(e1,4)/8.0)+(5*pow(e1,6)/16))/pow(1-(e1*e1),6.5));
	W1+=-(1/tf1)*(sOm1vh/(l1d))*((1+(3*pow(e1,2))+(3*pow(e1,4)/8))/pow(1-(e1*e1),5));

	X0=-m1*k0*pow(R0,5)*sOm0vh*sOm0ve/(mu1*l1d*pow(a1,5)*pow(1-e1*e1,2));
	X0+=-(sOm0vq*(1+(9.0*pow(e1,2)/2)+(5.0*pow(e1,4)/8)))/(2*l1d*tf0*pow(1-e1*e1,5));
	
	X1=-m0*k1*pow(R1,5)*sOm1vh*sOm1ve/(mu1*l1d*pow(a1,5)*pow(1-e1*e1,2));
	X1+=-(sOm1vq*(1+(9.0*pow(e1,2)/2)+(5.0*pow(e1,4)/8)))/(2*l1d*tf1*pow(1-e1*e1,5));
	
	Y0=-m1*k0*pow(R0,5)*sOm0vh*sOm0vq/(mu1*l1d*pow(a1,5)*pow(1-e1*e1,2));
	Y0+=(sOm0ve*(1+(3.0*pow(e1,2)/2)+(1.0*pow(e1,4)/8)))/(2*l1d*tf0*pow(1-e1*e1,5));

	Y1=-m0*k1*pow(R1,5)*sOm1vh*sOm1vq/(mu1*l1d*pow(a1,5)*pow(1-e1*e1,2));
	Y1+=(sOm1ve*(1+(3.0*pow(e1,2)/2)+(1.0*pow(e1,4)/8)))/(2*l1d*tf1*pow(1-e1*e1,5));

	Z0=((m1*k0*pow(R0,5))/(mu1*l1d*pow(a1,5)))*((2*pow(sOm0vh,2)-pow(sOm0vq,2)-pow(sOm0ve,2))/(2*pow(1-e1*e1,2)));
	Z0+=((m1*k0*pow(R0,5))/(mu1*l1d*pow(a1,5)))*((15*G*m1)/(a1*a1*a1))*((1+(3.0*pow(e1,2)/2)+(pow(e1,4)/8))/pow(1-e1*e1,5));

	Z1=((m0*k1*pow(R1,5))/(mu1*l1d*pow(a1,5)))*((2*pow(sOm1vh,2)-pow(sOm1vq,2)-pow(sOm1ve,2))/(2*pow(1-e1*e1,2)));
	Z1+=((m0*k1*pow(R1,5))/(mu1*l1d*pow(a1,5)))*((15*G*m0)/(a1*a1*a1))*((1+(3.0*pow(e1,2)/2)+(pow(e1,4)/8))/pow(1-e1*e1,5));

	dedt[0]=e1*((Z0+Z1)*qhat[0]-(Y0+Y1)*hhat[0]-(V0+V1)*ehat[0]);
	dedt[1]=e1*((Z0+Z1)*qhat[1]-(Y0+Y1)*hhat[1]-(V0+V1)*ehat[1]);
	dedt[2]=e1*((Z0+Z1)*qhat[2]-(Y0+Y1)*hhat[2]-(V0+V1)*ehat[2]);
	
	dhdt[0]=h1*((Y0+Y1)*ehat[0]-(X0+X1)*qhat[0]-(W0+W1)*hhat[0]);
	dhdt[1]=h1*((Y0+Y1)*ehat[1]-(X0+X1)*qhat[1]-(W0+W1)*hhat[1]);
	dhdt[2]=h1*((Y0+Y1)*ehat[2]-(X0+X1)*qhat[2]-(W0+W1)*hhat[2]);
	
	dsOm0dt[0]=(h1/I0)*((-Y0)*ehat[0]+(X0)*qhat[0]+(W0)*hhat[0]);
	dsOm0dt[1]=(h1/I0)*((-Y0)*ehat[1]+(X0)*qhat[1]+(W0)*hhat[1]);
	dsOm0dt[2]=(h1/I0)*((-Y0)*ehat[2]+(X0)*qhat[2]+(W0)*hhat[2]);

	dsOm1dt[0]=(h1/I1)*((-Y1)*ehat[0]+(X1)*qhat[0]+(W1)*hhat[0]);
	dsOm1dt[1]=(h1/I1)*((-Y1)*ehat[1]+(X1)*qhat[1]+(W1)*hhat[1]);
	dsOm1dt[2]=(h1/I1)*((-Y1)*ehat[2]+(X1)*qhat[2]+(W1)*hhat[2]);
	
	//printf("V0,W0,X0,Y0,Z0: %e,%e,%e,%e,%e\n",V0,W0,X0,Y0,Z0);	
	//printf("de1dt,dh1dt: %e,%e\n",-e1*V0,-h1*W0);	
	//printf("dsOm0dt0,dsOm0dt1,dsOm0dt1: %e,%e,%e\n",dsOm0dt[0],dsOm0dt[1],dsOm0dt[2]);	

	//printf("V1,W1,X1,Y1,Z1: %e,%e,%e,%e,%e\n",V1,W1,X1,Y1,Z1);	

}



int de_gr(double *dedt,double *e1vc,double *j1vc)
{
	double e1=sqrt(dotp(e1vc,e1vc));
	double j1=sqrt(dotp(j1vc,j1vc));
	
	double ehat[3];
	double jhat[3];
	double qhat[3];

	divideVec(e1vc,e1,ehat);
	divideVec(j1vc,j1,jhat);
	crossp(jhat,ehat,qhat);
	double zgr=3*pow(G,1.5)*pow(m0+m1,1.5)/(pow(a1,2.5)*c*c*(1-e1*e1));

	dedt[0]=zgr*e1*qhat[0];
	dedt[1]=zgr*e1*qhat[1];
	dedt[2]=zgr*e1*qhat[2];
}

double func_dphide1(double *dphide1,double *e1vc,double *e2vc,double *j1vc,double *j2vc)
{
	double e1=sqrt(dotp(e1vc,e1vc));
	double e2=sqrt(dotp(e2vc,e2vc));

	double e1vi,j1vi,e2vi,j2vi;
	int index;
	for(index=0;index<3;index++)
	{	
		e1vi=e1vc[index];
		j1vi=j1vc[index];
		e2vi=e2vc[index];
		j2vi=j2vc[index];
		dphide1[index]=(kphi0*mu1*(-12*e1vi*(1 - pow(e2,2)) + 30*j2vi*dotp(e1vc,j2vc)))/
   (8.*pow(1 - pow(e2,2),2.5)) + 
  (15*koct*kphi0*mu1*(dotp(e1vc,e2vc)*
        (16*e1vi*(1 - pow(e2,2)) - 70*j2vi*dotp(e1vc,j2vc)) + 
       10*j2vi*dotp(e2vc,j1vc)*dotp(j1vc,j2vc) + 
       e2vi*((-1 + 8*pow(e1,2))*(1 - pow(e2,2)) - 
          35*pow(dotp(e1vc,j2vc),2) + 5*pow(dotp(j1vc,j2vc),2))))/
   (64.*pow(1 - pow(e2,2),3.5));
	}
}


double func_dphidj1(double *dphidj1,double *e1vc,double *e2vc,double *j1vc,double *j2vc)
{
	double e1=sqrt(dotp(e1vc,e1vc));
	double e2=sqrt(dotp(e2vc,e2vc));

	double e1vi,j1vi,e2vi,j2vi;
	int index;
	for(index=0;index<3;index++)
	{	
		e1vi=e1vc[index];
		j1vi=j1vc[index];
		e2vi=e2vc[index];
		j2vi=j2vc[index];
		dphidj1[index]=(-3*j2vi*kphi0*mu1*dotp(j1vc,j2vc))/(4.*pow(1 - pow(e2,2),2.5)) + 
  (15*koct*kphi0*mu1*(10*j2vi*dotp(e1vc,j2vc)*dotp(e2vc,j1vc) + 
       10*j2vi*dotp(e1vc,e2vc)*dotp(j1vc,j2vc) + 
       10*e2vi*dotp(e1vc,j2vc)*dotp(j1vc,j2vc)))/
   (64.*pow(1 - pow(e2,2),3.5));
	}
}

double func_dphide2(double *dphide2,double *e1vc,double *e2vc,double *j1vc,double *j2vc)
{
	double e1=sqrt(dotp(e1vc,e1vc));
	double e2=sqrt(dotp(e2vc,e2vc));

	double e1vi,j1vi,e2vi,j2vi;
	int index;
	for(index=0;index<3;index++)
	{	
		e1vi=e1vc[index];
		j1vi=j1vc[index];
		e2vi=e2vc[index];
		j2vi=j2vc[index];
		dphide2[index]=(-3*kphi0*mu1*(25*(-1 + pow(e2,2))*(8*e2vi - 7*e1vi*koct)*
       pow(dotp(e1vc,j2vc),2) + 
      50*koct*dotp(e1vc,j2vc)*((-1 + pow(e2,2))*j1vi - 
         7*e2vi*dotp(e2vc,j1vc))*dotp(j1vc,j2vc) + 
      25*e2vi*koct*dotp(e1vc,e2vc)*
       ((-1 + 8*pow(e1,2))*(-1 + pow(e2,2)) + 
         49*pow(dotp(e1vc,j2vc),2) - 7*pow(dotp(j1vc,j2vc),2)) + 
      (-1 + pow(e2,2))*((-1 + pow(e2,2))*
          (8*(-1 + 6*pow(e1,2))*e2vi + 5*(1 - 8*pow(e1,2))*e1vi*koct) - 
         5*(8*e2vi - 5*e1vi*koct)*pow(dotp(j1vc,j2vc),2))))/
  (64.*pow(1 - pow(e2,2),4.5));
	}
}

double func_dphidj2(double *dphidj2,double *e1vc,double *e2vc,double *j1vc,double *j2vc)
{
	double e1=sqrt(dotp(e1vc,e1vc));
	double e2=sqrt(dotp(e2vc,e2vc));

	double e1vi,j1vi,e2vi,j2vi;
	int index;
	for(index=0;index<3;index++)
	{	
		e1vi=e1vc[index];
		j1vi=j1vc[index];
		e2vi=e2vc[index];
		j2vi=j2vc[index];
		dphidj2[index]=(-3*kphi0*mu1*(5*dotp(e1vc,j2vc)*
       (8*e1vi*(-1 + pow(e2,2)) + 35*e1vi*koct*dotp(e1vc,e2vc) - 
         5*j1vi*koct*dotp(e2vc,j1vc)) - 
      (8*(-1 + pow(e2,2))*j1vi + 25*j1vi*koct*dotp(e1vc,e2vc) + 
         25*e1vi*koct*dotp(e2vc,j1vc))*dotp(j1vc,j2vc)))/
  (32.*pow(1 - pow(e2,2),3.5));
	}
}

int func (double t, const double org[], double f[],
      void *params)
{
  	(void)(t); /* avoid unused parameter warning */

	double h1v[3],e1v[3];
	h1v[0]=org[0];e1v[0]=org[3];
	h1v[1]=org[1];e1v[1]=org[4];
	h1v[2]=org[2];e1v[2]=org[5];

	//printf("h1v[0]:%e h1v[1]:%e h1v[2]:%e\n",h1v[0],h1v[1],h1v[2]);
	//printf("e1v[0]:%e e1v[1]:%e e1v[2]:%e\n",e1v[0],e1v[1],e1v[2]);

	double h2v[3],e2v[3];
	h2v[0]=org[6];e2v[0]=org[9];
	h2v[1]=org[7];e2v[1]=org[10];
	h2v[2]=org[8];e2v[2]=org[11];

	updatePreFactors(e1v,h1v);

	double sOm0v[3],sOm1v[3];
	sOm0v[0]=org[12];sOm1v[0]=org[15];
	sOm0v[1]=org[13];sOm1v[1]=org[16];
	sOm0v[2]=org[14];sOm1v[2]=org[17];

	double dpde1v[3];
	double dpdj1v[3];
	double dpde2v[3];
	double dpdj2v[3];

	double e1,e2;
	e1=sqrt(dotp(e1v,e1v));
	e2=sqrt(dotp(e2v,e2v));
	
	double h1,h2;
	h1=sqrt(dotp(h1v,h1v));
	h2=sqrt(dotp(h2v,h2v));
	
	double j1v[3],j2v[3];
	divideVec(h1v,hk1,j1v);
	divideVec(h2v,hk2,j2v);
		
	func_dphide1(dpde1v,e1v,e2v,j1v,j2v);
	func_dphidj1(dpdj1v,e1v,e2v,j1v,j2v);
	func_dphide2(dpde2v,e1v,e2v,j1v,j2v);
	func_dphidj2(dpdj2v,e1v,e2v,j1v,j2v);
			
	double jcrossde1[3];double ecrossdj1[3];
	double ecrossde1[3];double jcrossdj1[3];

	double jcrossde2[3];double ecrossdj2[3];
	double ecrossde2[3];double jcrossdj2[3];

	crossp(j1v,dpde1v,jcrossde1);
	crossp(j1v,dpdj1v,jcrossdj1);
	crossp(e1v,dpde1v,ecrossde1);
	crossp(e1v,dpdj1v,ecrossdj1);
	
	crossp(j2v,dpde2v,jcrossde2);
	crossp(j2v,dpdj2v,jcrossdj2);
	crossp(e2v,dpde2v,ecrossde2);
	crossp(e2v,dpdj2v,ecrossdj2);
	
	f[0]=-(hk1/L1)*(jcrossdj1[0]+ecrossde1[0]);
	f[1]=-(hk1/L1)*(jcrossdj1[1]+ecrossde1[1]);
	f[2]=-(hk1/L1)*(jcrossdj1[2]+ecrossde1[2]);

	f[3]=-(1/L1)*(jcrossde1[0]+ecrossdj1[0]);
	f[4]=-(1/L1)*(jcrossde1[1]+ecrossdj1[1]);
	f[5]=-(1/L1)*(jcrossde1[2]+ecrossdj1[2]);
	
	f[6]=-(hk2/L2)*(jcrossdj2[0]+ecrossde2[0]);
	f[7]=-(hk2/L2)*(jcrossdj2[1]+ecrossde2[1]);
	f[8]=-(hk2/L2)*(jcrossdj2[2]+ecrossde2[2]);

	f[9] =-(1/L2)*(jcrossde2[0]+ecrossdj2[0]);
	f[10]=-(1/L2)*(jcrossde2[1]+ecrossdj2[1]);
	f[11]=-(1/L2)*(jcrossde2[2]+ecrossdj2[2]);

	//f[0]=0; f[1]=0; f[2]=0; 
	//f[3]=0; f[4]=0; f[5]=0;	

	//f[6]=0; f[7]=0; f[8]=0; 
	//f[9]=0; f[10]=0; f[11]=0;	

	f[12]=0; f[13]=0; f[14]=0; 
	f[15]=0; f[16]=0; f[17]=0;	
	//exit(0);

	double dedt[3],dhdt[3],dsOm0dt[3],dsOm1dt[3];
	de_gr(dedt,e1v,j1v);
	f[3]+=dedt[0];
	f[4]+=dedt[1];
	f[5]+=dedt[2];

	dedh_td(dedt,dhdt,dsOm0dt,dsOm1dt,e1v,h1v,sOm0v,sOm1v);	
	f[0]+=dhdt[0];
	f[1]+=dhdt[1];
	f[2]+=dhdt[2];	
	
	f[3]+=dedt[0];
	f[4]+=dedt[1];
	f[5]+=dedt[2];

	f[12]+=dsOm0dt[0];
	f[13]+=dsOm0dt[1];
	f[14]+=dsOm0dt[2];
	
	//printf("SPIN:: %e\n",sqrt(dotp(sOm1v,sOm1v)));

	f[15]+=dsOm1dt[0];
	f[16]+=dsOm1dt[1];
	f[17]+=dsOm1dt[2];

	return GSL_SUCCESS;
}


int integrateGSL(double *h1v,double *e1v,double *h2v,double *e2v,double *sOm0v,double *sOm1v,double tmax,FILE *fp)
{

	gsl_odeiv2_system sys = {func, NULL, 18, NULL};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-8, 1e-10, 0);
	double org[18];	
	double porg[18];	
	//printf("####-> START SPIN:: %e\n",sqrt(dotp(sOm1v,sOm1v)));

	int status;
	org[0]=h1v[0]; 	org[3]=e1v[0]; 
	org[1]=h1v[1]; 	org[4]=e1v[1]; 
	org[2]=h1v[2]; 	org[5]=e1v[2]; 
	
	org[6]=h2v[0]; 	org[9] =e2v[0]; 
	org[7]=h2v[1]; 	org[10]=e2v[1]; 
	org[8]=h2v[2]; 	org[11]=e2v[2];
	
	org[12]=sOm0v[0]; org[15]=sOm1v[0]; 
	org[13]=sOm0v[1]; org[16]=sOm1v[1]; 
	org[14]=sOm0v[2]; org[17]=sOm1v[2];

	double t=0;

	double dt=pow(a1,1.5)/10;
	//0.3*pow(a2/a1,3)*(m0/m2)*pow(a1,1.5)/5;
	fprintf(fp,"time,a1,a2,e1,e2,i1,omega1,Omega1,s0,s1\n");
	double ti=t+dt;
	int count=0;
	double h1,h2,e1,e2;
	double somega,comega,inc;
	double Omega,omega,cinc,sinc;
	double sv0[3],sv1[3];
	double s0,s1;
	double pe1;
	e1=sqrt(dotp(e1v,e1v));
	int td=0;
	double tsckz,tsctd,tscgr;
	double tscmn;
	double pt;

	tsckz=pow(a2/a1,3)*((m0+m1)/m2)*pow(a1*a1*a1/(m0+m1),0.5);
	tsctd=tf0;
	tscgr=pow((3*pow(G,1.5)*pow(m0+m1,1.5)/(pow(a1,2.5)*c*c)),-1);
	tscmn=tsckz<tsctd?tsckz:tsctd;
	tscmn=tscgr<tscmn?tscgr:tscmn;	
	dt=tscmn/1000;
	
	while(t<tmax)
	{
		porg[0]=org[0]; porg[3]=org[3]; porg[9]=org[9];   porg[12]=org[12]; porg[15]=org[15];
		porg[1]=org[1]; porg[4]=org[4]; porg[10]=org[10]; porg[13]=org[13]; porg[16]=org[16];
		porg[2]=org[2]; porg[5]=org[5]; porg[11]=org[11]; porg[14]=org[14]; porg[17]=org[17];
						
		porg[6]=org[6]; 
		porg[7]=org[7]; 
		porg[8]=org[8];
		pt=t;
		pe1=e1;	
		status = gsl_odeiv2_driver_apply (d, &t, ti, org);
		if (status != GSL_SUCCESS)
		  break;
		
		h1v[0]=org[0];e1v[0]=org[3];
		h1v[1]=org[1];e1v[1]=org[4];
		h1v[2]=org[2];e1v[2]=org[5];
		
		h2v[0]=org[6];e2v[0]=org[9];
		h2v[1]=org[7];e2v[1]=org[10];
		h2v[2]=org[8];e2v[2]=org[11];	

		sv0[0]=org[12];sv1[0]=org[15];
		sv0[1]=org[13];sv1[1]=org[16];
		sv0[2]=org[14];sv1[2]=org[17];	

		s0=sqrt(dotp(sv0,sv0));
		s1=sqrt(dotp(sv1,sv1));
		h1=sqrt(dotp(h1v,h1v));
		h2=sqrt(dotp(h2v,h2v));
                e1=sqrt(dotp(e1v,e1v));
		e2=sqrt(dotp(e2v,e2v));

		if(fabs(pe1-e1)>0.001)
		{
			dt=dt/2;
			t=pt;
			ti=t+dt;
			org[0]=porg[0]; org[3]=porg[3]; org[9]=porg[9];   org[12]=porg[12]; org[15]=porg[15];
			org[1]=porg[1]; org[4]=porg[4]; org[10]=porg[10]; org[13]=porg[13]; org[16]=porg[16];
			org[2]=porg[2]; org[5]=porg[5]; org[11]=porg[11]; org[14]=porg[14]; org[17]=porg[17];
						
			org[6]=porg[6]; 
			org[7]=porg[7]; 
			org[8]=porg[8];
			gsl_odeiv2_driver_reset(d);
			continue;			
		}

		ti=t+dt;
		
		//printf("kz: %e, td: %e, tscgr:%e tscmn: %e\n",tsckz,tsctd,tscgr,tscmn);
		//printf("dt: %e\n",dt);
		//exit(0);	
		
		Omega=atan2(h1v[0],-h1v[1]);
		sinc=(h1v[0]*sin(Omega))-(h1v[1]*cos(Omega));
		cinc=h1v[2];
		inc=atan2(sinc,cinc);
		somega=((-e1v[0]*sin(Omega))+(e1v[1]*cos(Omega)))/cos(inc);
		comega=(e1v[0]*cos(Omega))+(e1v[1]*sin(Omega));
		omega=atan2(somega,comega);
		
		if(count%10==0)
		{
			fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",t,a1,a2,e1,e2,inc,omega,Omega,(2*PI/s0)/DAY,((2*PI)/s1)/DAY);
		}

		if(a1*(1-e1)<(0.0127))
		{
			td=1;
			break;
		}
		if(a1<0.1 && e1<0.01)
		{
			break;
		}
		
		count++;
	}
	//fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%d\n",t,a1,a2,e1,e2,inc,omega,Omega,td);
	gsl_odeiv2_driver_free (d);
	

}

int main(int argn,char *argv[])
{

	double e1,e2;
	double inc1,omega1,Omega1;
	double inc2,omega2,Omega2;

	char *fname=argv[1];	
	const gsl_rng_type * T;
	gsl_rng * r;
        gsl_rng_env_setup();
	T = gsl_rng_default; // Generator setup
    	r = gsl_rng_alloc (T);
   	unsigned long mySeed =getpid() + (time(0) << (CHAR_BIT * sizeof(pid_t)));
	gsl_rng_set(r, mySeed);

        a1=5.0;
	a2=1000.0;//exp(gsl_ran_flat(r,log(100),log(1500)));
	e1=0.1;//gsl_ran_rayleigh(r,0.01);
	inc1=85.6*PI/180;
	omega1=45.0*PI/180;//gsl_ran_flat(r,0,2*PI);
	Omega1=0.0;//gsl_ran_flat(r,0,2*PI);
	
	e2=0.5;//sqrt(gsl_ran_flat(r,0,0.81));
	inc2=0.0001;//acos(gsl_ran_flat(r,-1,1))+0.01;
	omega2=0.0001;//gsl_ran_flat(r,0,2*PI);
	Omega2=0.0001;//Omega1;

	m0=1.1;
	m1=7.8*0.0009543;
	m2=1.1;

	R0=0.00465047;
	R1=0.000467329598;
	
	I0=(0.08*m0*R0*R0);
	I1=(0.25*m1*R1*R1);

	klove0=0.028;
	klove1=0.50;

	k0=klove0/2;
	k1=klove1/2;
	
	tv0=50;
	tv1=0.1;

	Q0=1e6;
	Q1=3e5;
	
	L1=m0*m1*sqrt(G*(m0+m1)*a1)/(m0+m1);
	L2=(m0+m1)*m2*sqrt(G*(m0+m1+m2)*a2)/(m0+m1+m2);
	mu1=m0*m1/(m0+m1);
	kphi0=G*m2*a1*a1/(a2*a2*a2);
	koct=((m0-m1)/(m0+m1))*(a1/a2);
	
	tf0=(tv0/9)*pow(a1/R0,8)*((m0*m0)/((m0+m1)*m1))*pow(1+2*k0,-2);
	tf1=(tv1/9)*pow(a1/R1,8)*((m1*m1)/((m1+m0)*m0))*pow(1+2*k1,-2);

	hk1=m0*m1*sqrt(G*(m0+m1)*a1)/(m0+m1);
	hk2=(m0+m1)*m2*sqrt(G*(m0+m1+m2)*a2)/(m0+m1+m2);

	double h1mod=hk1*sqrt(1-(e1*e1));
	double h2mod=hk2*sqrt(1-(e2*e2));

	double h1v[3]={h1mod*sin(Omega1)*sin(inc1),-h1mod*cos(Omega1)*sin(inc1),h1mod*cos(inc1)};
	double e1v[3];
	e1v[0]=e1*((cos(Omega1)*cos(omega1))-(sin(Omega1)*sin(omega1)*cos(inc1)));
	e1v[1]=e1*((sin(Omega1)*cos(omega1))+(cos(Omega1)*sin(omega1)*cos(inc1)));
	e1v[2]=e1*(sin(omega1)*sin(inc1));


	double h2v[3]={h2mod*sin(Omega2)*sin(inc2),-h2mod*cos(Omega2)*sin(inc2),h2mod*cos(inc2)};
	double e2v[3];
	e2v[0]=e2*((cos(Omega2)*cos(omega2))-(sin(Omega2)*sin(omega2)*cos(inc2)));
	e2v[1]=e2*((sin(Omega2)*cos(omega2))+(cos(Omega2)*sin(omega2)*cos(inc2)));
	e2v[2]=e2*(sin(omega2)*sin(inc2));

	double sOm1v[3];
	double s0mod,s1mod;
	s0mod=2*PI/(20*DAY);
	s1mod=2*PI/(10*HOUR);

	
	double sOm0v[3];
	sOm0v[0]=s0mod*h1v[0]/h1mod;  sOm1v[0]=s1mod*h1v[0]/h1mod;;
	sOm0v[1]=s0mod*h1v[1]/h1mod;  sOm1v[1]=s1mod*h1v[1]/h1mod;;
	sOm0v[2]=s0mod*h1v[2]/h1mod;  sOm1v[2]=s1mod*h1v[2]/h1mod;;

	double dedt[3],dhdt[3],dsOm0dt[3],dsOm1dt[3];
	double tmax=gsl_ran_flat(r,0,1e10);
	
	double lhs=a2/a1;
	double muc=m2/(m0+m1);
	double rhs=2.8*pow(1+muc,0.4)*(pow(1+e2,0.4)/pow(1-e2,1.2))*(1-(0.3*inc2/PI));
	if(lhs>rhs)
	{
		FILE *fp=fopen(fname,"w");
		printf("a1:%e,e1:%e,a2:%e,e2:%e,tmax:%e\n",a1,e1,a2,e2,tmax);
	//exit(0);
		integrateGSL(h1v,e1v,h2v,e2v,sOm0v,sOm1v,tmax,fp);
		fclose(fp);
	}

	gsl_rng_free (r);	
}



