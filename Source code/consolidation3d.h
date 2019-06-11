#ifndef CONSOLIDATION3D_H
#define CONSOLIDATION3D_H
#define EIGEN_USE_MKL_ALL
#include <QDebug>
#include <QString>
#include <QFile>
#include <QFileDialog>
#include <QElapsedTimer>
#include <QObject>
#include <iostream>
#include <fstream>

#include "getfile.h"
#include "writetofile.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/PardisoSupport>


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <QThread>
#include "pvdlib.h"
#include "importmeshdata.h"
#include "soilparameter.h"
#include "pvdparameters.h"
#include "analysistype.h"

using namespace std;
using namespace Eigen;

class Consolidation3D:public QObject
{
    Q_OBJECT

public:
    Consolidation3D();
    //constant
    const double PI=3.14159265359;
    const double machine_error=2.2204e-16;

    //Functions
    void importData();
    void initialData();
    void resetData();
    void applyVerticalLoad();
    void createBoundaryCondition();
    void assemblyGlobalMatrix();
    void solveDirect();    
    void exportResults();
    void calculateAnalytical();
    void calculateError();
    double getError();
    void WriteLogFile();
    void clearData();
    MatrixXd GetXDisp(){return U;}
    MatrixXd GetYDisp(){return V;}
    MatrixXd GetZDisp(){return W;}
    MatrixXd GetPore(){return PFull;}

    //Common functions
    void undrainedSolve();
    void drainedSolve();
    void calculateStress(int &eleNum);
    void createGauss4(); //For tet
    void createGauss8(); //For hexagon 8-nodes 2x2x2
    void createGauss9(); //For prism
    void createGauss27(); //For hexagon 20-nodes 3x3x3 + pyramid
    void createGauss13(); //For pyramid
    void pauseSystem();

    void compareMatrix(Ref<MatrixXi> matrixA, Ref<MatrixXi> matrixB, Ref<MatrixXi> matrixC);
    void calculateNodfat();
    void getPore();
    int countPoreBoundary(Ref<MatrixXd> matrixA,Ref<MatrixXi> nodfmt);
    void createNewPoreBoundary(Ref<MatrixXd>matrixOld, Ref<MatrixXd>matrixNew,Ref<MatrixXi> nodfmt);
    void createGlobalBoundary(Ref<MatrixXd>LocalMatrix, Ref<MatrixXd>GlobalMatrix);
    double calculationError();
    void setCd(double CdValue);

    //element Functions
    void tet10pMatrix(int &eleNum);
    void hex20pMatrix(int &eleNum);
    void pyra13pMatrix(int &eleNum);
    void prism15pMatrix(int &eleNum);
    void drain1DMatrix(int &eleNum);
    double error, error0, error1;

public slots:
    void GetMeshData(meshDataNames meshName);
    void GetSoilData(Soil soilData);
    void GetPVDData(PVD PVDData);
    void GetAnalysisType(SolutionType typeData);

signals:
    void sendCalculationInfor(int step,int numberOfStep);

private:
    QElapsedTimer timer;
    //---------------------
    double kv,kh,K,G,poissionRatio;
    double kw, Aw;
    double Ke, Ge, kve,khe,ve,Se;
    double voidRatio,re,Cd;
    double Ss;

    double H;
    int analysisType;

    QString fileName, folderName;
    QString coordFile, eleFile, fixXFile, fixYFile, fixZFile, fixHFile,forceYFile;
    int noe,non,numberOfStep,dof,doff,totalDof,step;
    double sigma1,sigma3,sigma0,timeInterval, timeInterval0;
    int eleNum;
    double p0;

    double A,B,Cf,Cs,gf;
    double tol, maxIt;      
    double Cd0, Cd1;

    double R,Rw,Rs,ks,ksRatio,n,Cvr,s;
    bool smearAna, noSmearAna;
    //---------------------
    MatrixXd coordinates,elements,fixx,fixy,fixz,fixh,forcey,gauss,pizo,analytical,poreModel;
    MatrixXi elementsIndex;
    MatrixXd gauss4, gauss8, gauss27, gauss13, gauss9;
    MatrixXd fixAll,fixhNew;
    MatrixXd DirichletAll, DirichletAllU, DirichletAll1;
    MatrixXd DirichletU, Dirichlet1, Dirichlet;

    MatrixXi tet10p,hex20p,hex8p,pyramid13p, prism15p, drain1Dp;
    MatrixXi nodfat,nodfmt;
    MatrixXd X0, XX, X, F;
    MatrixXd U,V,W,P,PFull;
    MatrixXd Perror, Panalytical;
    MatrixXd Pizometer;
    double errorNew;
    //---------------------
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trip_total;
    SparseMatrix<double,RowMajor> KK;

    //---------------------
    GetFile getfile;
    WriteToFile exportFile;

    //---------------------
    PvdLib analyTest;
    meshDataNames meshName;
    Soil soilData;
    PVD PVDData;
    SolutionType typeData;

};

#endif // CONSOLIDATION3D_H
