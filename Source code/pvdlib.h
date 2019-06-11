#ifndef PVDLIB_H
#define PVDLIB_H
#include <math.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>
#include <QObject>
#include <QDebug>

using namespace std;
using namespace boost::math;

class PvdLib : public QObject
{
    Q_OBJECT
public:
    explicit PvdLib(QObject *parent = nullptr);
    //math function
    double factorial (int n);
    double besselj(double u,double x);
    double bessely(double u,double x);

    //PVD no smear function - freestrain
    void setParametersFreeStrain(double re,double rw,double Cvr,double p0);
    double PvdNoSmearBessel(double x,double n);
    double rootFunctionNoSmearBiSection(double initial,double sizeStep,double n);
    void findRootNoSmear(vector<double> &alpha,double n);
    double freeStrain(double r,double time);
    void printInforNoSmear();

    //PVD smear zone function - freestrain
    void setParametersFreeStrainSmear(double re, double rw, double rs, double Cvr,double kh, double ks,double p0);
    double PvdSmearBessel(double x,double kh,double ks,double s,double n);
    double rootFunctionSmearBiSection(double initial,double sizeStep,double kh, double ks, double s, double n);
    void findRootSmear(vector<double> &alpha, double kh, double ks, double s, double n);
    double freeStrainSmear(double r,double time);
    void printInforSmear();
    
    //PVD equal strain - no smear
    void setParametersEqualStrain(double re,double rw,double Cvr,double p0);
    double equalStrain(double r,double time);

    //PVD no smear zone - equal strain
    void setParametersEqualStrainSmear(double re,double rw,double rs,double Cvr,double kh,double ks,double p0);
    double equalStrainSmear(double r,double time);

    //PVD smear zone + well resistance - equal strain
    void setParametersEqualStrainResistance(double re,double rw,double rs,double Cvr,double kh,double ks,double p0,double kw,double H);
    double equalStrainResistance(double r,double z,double time);

    //support function
    void pauseSystem();
    void setParameterMechanic(double ve,double gf,double kh);
    void printInformation();
    double getK();
    double getG();

signals:

public slots:

private:    
    double re=0.565;
    double rw=0.0264;
    double n=re/rw; //Not depend on r
    double dt=86400;
    double r=0.565;
    double p0=100;
    double Cvr=1e-7;

    double kh=1e-9;
    double ks=0.5*kh;
    double kw=1e-4;
    double H=1;
    double rs=0.102;
    double s=rs/rw;

    double ve=0.15;
    double gf=9.81;
    double Ke,Ge;

    vector<double> alpha;
    int analysisType=0;
    bool rootNoSmearCheck=false;
    bool rootSmearCheck=false;
};

#endif // PVDLIB_H
