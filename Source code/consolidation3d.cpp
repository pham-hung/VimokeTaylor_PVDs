#include "consolidation3d.h"

Consolidation3D::Consolidation3D()
{
    cout<<"---------------------------------------"<<endl;
}

void Consolidation3D::importData()
{
    getfile.DoGetFile(coordFile);
    coordinates=MatrixXd::Zero(getfile.row,getfile.col);
    coordinates=getfile.data_file;

    getfile.DoGetFile(eleFile);
    elements=MatrixXd::Zero(getfile.row,getfile.col);
    elements=getfile.data_file;

    getfile.DoGetFile(fixXFile);
    fixx=MatrixXd::Zero(getfile.row,getfile.col);
    fixx=getfile.data_file;

    getfile.DoGetFile(fixYFile);
    fixy=MatrixXd::Zero(getfile.row,getfile.col);
    fixy=getfile.data_file;


    getfile.DoGetFile(fixZFile);
    fixz=MatrixXd::Zero(getfile.row,getfile.col);
    fixz=getfile.data_file;

    getfile.DoGetFile(fixHFile);
    fixh=MatrixXd::Zero(getfile.row,getfile.col);
    fixh=getfile.data_file;

    getfile.DoGetFile(forceYFile);
    forcey=MatrixXd::Zero(getfile.row,getfile.col);
    forcey=getfile.data_file;
}

void Consolidation3D::initialData()
{
    createGauss4();
    createGauss8();
    createGauss9();
    createGauss13();
    createGauss27();
    //---------------------------------
    step=0;
    A=1; B=1;
    Cs=0; Cf=1e-7;
    double n=voidRatio/(1+voidRatio);
    gf=9.81;
    Ss=n*Cf+(A-n)*Cs;
    tol=1e-4;
    maxIt=10000;
    //---------------------------------
    non=coordinates.rows();
    noe=elements.col(0).maxCoeff();
    elementsIndex=MatrixXi::Zero(noe,1);
    int ecount=0;
    for (int j=0;j<elements.rows();j++)
    {
        int eleIndex=elements(j,0);
        if (eleIndex!=0)
        {
            elementsIndex(ecount,0)=j;
            ecount=ecount+1;
        }
    }
    nodfat=MatrixXi::Zero(non,4);
    nodfmt=MatrixXi::Zero(non,2);
    calculateNodfat();
    //---------------------------------
    totalDof=0;
    for (int j=0;j<non;j++)
    {
        nodfmt(j,1)=nodfat.row(j).sum();
        totalDof=totalDof+nodfat.row(j).sum();
        if(j<(non-1)){nodfmt(j+1,0)=nodfmt(j,0)+nodfat.row(j).sum();}
    }
    cout<<"Total DOF:"<<totalDof<<endl;
    cout<<"Number of node:"<<non<<endl;
    cout<<"Number of element: "<<noe<<endl;

    //---------------------------------
    F=MatrixXd::Zero(totalDof,1); //Right hand side vector
    X=MatrixXd::Zero(totalDof,numberOfStep); //Solution vector
    XX=MatrixXd::Zero(totalDof,1);
    X0=MatrixXd::Zero(totalDof,1);  //Solution vector for initial step;
    U=MatrixXd::Zero(non,numberOfStep);
    V=MatrixXd::Zero(non,numberOfStep);
    W=MatrixXd::Zero(non,numberOfStep);
    PFull=MatrixXd::Zero(non,numberOfStep);
    P=MatrixXd::Zero(totalDof-3*non,numberOfStep);
    Pizometer=MatrixXd::Zero(numberOfStep,pizo.rows());
    //    for (int j=0;j<pizo.rows();j++)
    //    {
    //        pizo(j,0)=nodfmt(pizo(j,0)-1,0)+3;
    //    }
    poreModel=MatrixXd::Zero(numberOfStep,1);
    cout<<"Finish initialData"<<endl;
    //---------
    Perror=MatrixXd::Zero(non,numberOfStep);
    Panalytical=MatrixXd::Zero(non,numberOfStep);

}

void Consolidation3D::resetData()
{
    F=MatrixXd::Zero(totalDof,1); //Right hand side vector
    X=MatrixXd::Zero(totalDof,numberOfStep); //Solution vector
    XX=MatrixXd::Zero(totalDof,1);
    X0=MatrixXd::Zero(totalDof,1);  //Solution vector for initial step;
    U=MatrixXd::Zero(non,numberOfStep);
    V=MatrixXd::Zero(non,numberOfStep);
    W=MatrixXd::Zero(non,numberOfStep);
    PFull=MatrixXd::Zero(non,numberOfStep);
    P=MatrixXd::Zero(totalDof-3*non,numberOfStep);
    Perror=MatrixXd::Zero(non,numberOfStep);
}

void Consolidation3D::applyVerticalLoad()
{
    //Sigma 1 from forcey matrix
    //forcey matrix is vector load from load = 100 Pa
    double scale=sigma1/100;
    for (int j=0;j<forcey.rows();j++)
    {
        int node=forcey(j,0)-1; //mode index start from 0
        int equaNum=nodfmt(node,0)+1;
        F(equaNum,0)=scale*forcey(j,1);
    }

}

void Consolidation3D::createBoundaryCondition()
{
    int count=0;
    MatrixXd fixxNew=MatrixXd::Zero(fixx.rows(),fixx.cols());
    MatrixXd fixyNew=MatrixXd::Zero(fixy.rows(),fixy.cols());
    MatrixXd fixzNew=MatrixXd::Zero(fixz.rows(),fixz.cols());

    for (int j=0;j<fixx.rows();j++)
    {
        int index =fixx(j,0)-1;
        fixxNew(j,0)=nodfmt(index,0)+0;
        fixxNew(j,1)=fixx(j,1);
    }

    for (int j=0;j<fixy.rows();j++)
    {
        int index=fixy(j,0)-1;
        fixyNew(j,0)=nodfmt(index,0)+1;
        fixyNew(j,1)=fixy(j,1);
    }

    for (int j=0;j<fixz.rows();j++)
    {
        int index=fixz(j,0)-1;
        fixzNew(j,0)=nodfmt(index,0)+2;
        fixzNew(j,1)=fixz(j,1);
    }
    count=countPoreBoundary(fixh,nodfmt);
    fixhNew=MatrixXd::Zero(count,fixh.cols());
    createNewPoreBoundary(fixh,fixhNew,nodfmt);

    Dirichlet1=MatrixXd::Zero(fixx.rows()+fixy.rows()+fixz.rows()+fixhNew.rows(),fixx.cols()); //h1
    DirichletU=MatrixXd::Zero(fixx.rows()+fixy.rows()+fixz.rows(),fixx.cols()); //undrained first step
    Dirichlet=MatrixXd::Zero(1,1);

    Dirichlet1<<fixxNew,fixyNew,fixzNew,fixhNew;
    DirichletU<<fixxNew,fixyNew,fixzNew;

    DirichletAll=MatrixXd::Zero(totalDof,3);
    DirichletAllU=MatrixXd::Zero(totalDof,3);
    DirichletAll1=MatrixXd::Zero(totalDof,3);

    createGlobalBoundary(DirichletU,DirichletAllU);
    createGlobalBoundary(Dirichlet1,DirichletAll1);
}

void Consolidation3D::assemblyGlobalMatrix()
{
    //3 Threads + 1 main Thread for 4 cores;
    Ke=K;
    ve=poissionRatio;
    Ge=3*Ke*(1-2*ve)/2/(1+ve);
    Se=Ss;
    kve=kv;
    khe=re*kve;
    qDebug()<<"Cd: "<<Cd<<endl;

    //main thread
    trip_total.clear();
    trip_total.resize(0);
    trip_total.reserve(30*30*noe);

    for (int i=0;i<noe;i++)
    {
        eleNum=i;
        kve=kv;
        khe=re*kve;
        int eleIndex=elementsIndex(eleNum,0);
        int eleMat=elements(eleIndex,12);

        if(eleMat==2)
        {
            khe=Cd*khe;
        }

        int nodeCount=elements(eleIndex,11);
        if(nodeCount==20){hex20pMatrix(eleNum);}
        else if(nodeCount==10){tet10pMatrix(eleNum);}
        else if(nodeCount==13){pyra13pMatrix(eleNum);}
        else if(nodeCount==15){prism15pMatrix(eleNum);}
        else if(nodeCount==2){drain1DMatrix(eleNum);}
    }

    //Create global sparse matrix
    KK.setFromTriplets(trip_total.begin(),trip_total.end());

    //Apply boundary condition, global level
    for (int j=0;j<Dirichlet.rows();j++)
    {
        int jj=Dirichlet(j,0);
        F(jj,0)=Dirichlet(j,1);
        KK.coeffRef(jj,jj)=1;
    }
    KK.prune(0.0);
    KK.makeCompressed();
    trip_total.clear(); //Release memory
}

void Consolidation3D::solveDirect()
{   
    PardisoLU<SparseMatrix<double,RowMajor> > solver;
    solver.compute(KK);
    XX.setZero();
    XX=solver.solve(F);
}

void Consolidation3D::exportResults()
{
    fileName=folderName+"/"+"X-Displacement.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(U);

    fileName=folderName+"/"+"Y-Displacement..txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(V);

    fileName=folderName+"/"+"Z-Displacement..txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(W);

    fileName=folderName+"/"+"Pore-Pressure.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(PFull);

    fileName=folderName+"/"+"Pore-Analytical-Solution.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Panalytical);

    WriteLogFile();
}

void Consolidation3D::calculateAnalytical()
{
    double distance;
    if(analysisType==0)
    {
        analyTest.setParametersEqualStrain(R,Rw,Cvr,p0);
    }
    else if (analysisType==1)
    {
        analyTest.setParametersEqualStrainSmear(R,Rw,Rs,Cvr,kh,ks,p0);
    }
    else if(analysisType==2)
    {
        analyTest.setParametersEqualStrainResistance(R,Rw,Rs,Cvr,kh,ks,p0,kw,H);
    }
    else if (analysisType==3)
    {
        analyTest.setParametersFreeStrain(R,Rw,Cvr,p0);
    }
    else if(analysisType==4)
    {
        analyTest.setParametersFreeStrainSmear(R,Rw,Rs,Cvr,kh,ks,p0);
    }

    double minY=coordinates.col(2).minCoeff();
    double maxY=coordinates.col(2).maxCoeff();

    for (int i=0;i<numberOfStep;i++)
    {
        for (int j=0;j<non;j++)
        {
            if(nodfmt(j,1)==4) //calculate pore presure
            {
                double xCoord=coordinates(j,1);
                double zCoord=coordinates(j,3);
                double yCoord=coordinates(j,2);
                double z=maxY-yCoord;
                distance=sqrt(xCoord*xCoord+zCoord*zCoord);

                if(analysisType==0)
                {
                    Panalytical(j,i)=analyTest.equalStrain(distance,i*timeInterval0);
                }
                else if(analysisType==1)
                {
                    Panalytical(j,i)=analyTest.equalStrainSmear(distance,i*timeInterval0);
                }
                else if(analysisType==2)
                {
                    Panalytical(j,i)=analyTest.equalStrainResistance(distance,z,i*timeInterval0);
                }
                else if(analysisType==3)
                {
                    Panalytical(j,i)=analyTest.freeStrain(distance,i*timeInterval0);
                }
                else if(analysisType==4)
                {
                    Panalytical(j,i)=analyTest.freeStrainSmear(distance,i*timeInterval0);
                }
            }
        }
    }

    fileName=folderName+"/"+"Panalytical.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Panalytical);
}

void Consolidation3D::calculateError()
{
    int count=0;
    double errorTotal=0;
    double errorEle=0;
    for(int i=0;i<PFull.rows();i++)
    {
        for (int j=0;j<PFull.cols();j++)
        {
            if(nodfmt(i,1)==4)
            {
                double poreModel=PFull(i,j);
                double poreAnalytical=Panalytical(i,j);
                if(poreAnalytical==0)
                {
                    Perror(i,j)=0;
                }
                else
                {
                    errorEle=abs(100*(poreAnalytical-poreModel)/poreAnalytical);
                    Perror(i,j)=errorEle;
                    count=count+1;
                    errorTotal=errorTotal+errorEle;
                }
            }
        }
    }
    errorNew=errorTotal/count;
    cout<<"Average error: "<<errorNew<<endl;

}

double Consolidation3D::getError()
{
    return errorNew;
}

void Consolidation3D::WriteLogFile()
{
    fileName=folderName+"/"+"log.txt";
    ofstream file_;
    file_.open(fileName.toStdString());
    file_<<"The found Cd : "<<Cd<<'\n';
    file_<<"Minimum error: "<<errorNew<<'\n';
    file_<<"Minimum settlement: "<<V.col(numberOfStep-1).minCoeff()<<'\n';

    file_<<"Re, Rw, Rs: "<<R<<'\t'<<Rw<<'\t'<<Rs<<'\n';
    file_<<"Cvr,kv,kh/kv,ks/kh,p0,step,dt,Poisson's ratio: "<<Cvr<<'\t'<<kv<<'\t'<<re<<'\t'<<ksRatio<<'\t'<<p0<<'\t'<<numberOfStep<<'\t'<<timeInterval0<<'\t'<<poissionRatio<<'\n';
    file_<<"kw,Aw: "<<kw<<'\t'<<Aw<<'\n';
    file_<<"Analysis Type "<<analysisType<<'\n';
    file_<<"Analysis Type == 0: Equal strain, no smear-zone"<<'\n';
    file_<<"Analysis Type == 1: Equal strain, with smear-zone"<<'\n';
    file_<<"Analysis Type == 2: Equal strain, with well resistance"<<'\n';
    file_<<"Analysis Type == 3: Free strain, no smear-zone"<<'\n';
    file_<<"Analysis Type == 4: Free strain, with smear-zone"<<'\n';
    file_.close();
}

void Consolidation3D::clearData()
{
    F=MatrixXd::Zero(0,0); //Right hand side vector
    X=MatrixXd::Zero(0,0); //Solution vector
    XX=MatrixXd::Zero(0,0);
    X0=MatrixXd::Zero(0,1);  //Solution vector for initial step;
    U=MatrixXd::Zero(0,0);
    V=MatrixXd::Zero(0,0);
    W=MatrixXd::Zero(0,0);
    PFull=MatrixXd::Zero(0,0);
    P=MatrixXd::Zero(0,0);
    Perror=MatrixXd::Zero(0,0);
}

void Consolidation3D::undrainedSolve()
{
    sigma1=p0;
    timeInterval=0;
    step=0;
    cout<<"Pressure: "<<sigma1<<endl;
    Dirichlet.setZero();
    Dirichlet.resize(DirichletU.rows(),DirichletU.cols());
    Dirichlet=DirichletU;
    DirichletAll.setZero();
    DirichletAll=DirichletAllU;
    X0.setZero();
    trip_total.clear();
    trip_total.reserve(30*30*noe);
    KK.setZero();
    KK.resize(totalDof,totalDof);
    F.setZero();
    applyVerticalLoad();
    assemblyGlobalMatrix();
    solveDirect();
    X.col(step)=XX.col(0);

    int pore_count=0;
    for (int j=0;j<non;j++)
    {
        U(j,step)=XX(nodfmt(j,0)+0,0);
        V(j,step)=XX(nodfmt(j,0)+1,0);
        W(j,step)=XX(nodfmt(j,0)+2,0);
        if(nodfmt(j,1)==4)
        {
            pore_count=pore_count+1;
            P(pore_count-1,step)=XX(nodfmt(j,0)+3,0);
            PFull(j,step)=XX(nodfmt(j,0)+3,0);
        }
    }
    getPore();
}

void Consolidation3D::drainedSolve()
{
    for (int st=1;st<numberOfStep;st++)
    {
        step=st;
        sigma1=p0;
        timeInterval=timeInterval0;
        X0=X.col(step-1);

        Dirichlet.setZero();
        Dirichlet.resize(Dirichlet1.rows(),Dirichlet1.cols());
        Dirichlet=Dirichlet1;
        DirichletAll.setZero();
        DirichletAll=DirichletAll1;

        trip_total.clear();
        trip_total.reserve(30*30*noe);
        KK.setZero();
        KK.resize(totalDof,totalDof);
        F.setZero();
        applyVerticalLoad();
        assemblyGlobalMatrix();
        //solveIterative();
        solveDirect();
        X.col(step)=XX.col(0);

        int pore_count=0;
        for (int j=0;j<non;j++)
        {
            U(j,step)=XX(nodfmt(j,0)+0,0);
            V(j,step)=XX(nodfmt(j,0)+1,0);
            W(j,step)=XX(nodfmt(j,0)+2,0);
            if(nodfmt(j,1)==4)
            {
                pore_count=pore_count+1;
                P(pore_count-1,step)=XX(nodfmt(j,0)+3,0);
                PFull(j,step)=XX(nodfmt(j,0)+3,0);
            }
        }
        getPore();
        cout<<"Step: "<<step<<endl;
        cout<<"Settlement: "<<V.col(step).minCoeff()<<endl;
        emit sendCalculationInfor(step,numberOfStep);
    }
    //    poreModel.col(0)=Pizometer.col(0);
    if(pizo.rows()>=1)
    {
        fileName=folderName+"/"+"Pmodel.txt";
        exportFile.fileName=fileName.toStdString();
        exportFile.ToFile(poreModel);
    }
    WriteLogFile();
}

void Consolidation3D::createGauss4()
{
    gauss4=MatrixXd::Zero(4,4);
    gauss4(0,0)=0.585410;
    gauss4(0,1)=0.138197;
    gauss4(0,2)=0.138197;
    gauss4(0,3)=0.041667;

    gauss4(1,0)=0.138197;
    gauss4(1,1)=0.585410;
    gauss4(1,2)=0.138197;
    gauss4(1,3)=0.041667;

    gauss4(2,0)=0.138197;
    gauss4(2,1)=0.138197;
    gauss4(2,2)=0.585410;
    gauss4(2,3)=0.041667;

    gauss4(3,0)=0.138197;
    gauss4(3,1)=0.138197;
    gauss4(3,2)=0.138197;
    gauss4(3,3)=0.041667;

}

void Consolidation3D::createGauss8()
{
    //2x2x2, for brick 8 nodes
    gauss8=MatrixXd::Zero(8,4);
    double gaussValue=0.57735;
    gauss8(0,0)=-gaussValue;
    gauss8(0,1)=-gaussValue;
    gauss8(0,2)=-gaussValue;
    gauss8(0,3)=1;

    gauss8(1,0)=gaussValue;
    gauss8(1,1)=-gaussValue;
    gauss8(1,2)=-gaussValue;
    gauss8(1,3)=1;

    gauss8(2,0)=gaussValue;
    gauss8(2,1)=-gaussValue;
    gauss8(2,2)=gaussValue;
    gauss8(2,3)=1;

    gauss8(3,0)=-gaussValue;
    gauss8(3,1)=-gaussValue;
    gauss8(3,2)=gaussValue;
    gauss8(3,3)=1;

    gauss8(4,0)=-gaussValue;
    gauss8(4,1)=gaussValue;
    gauss8(4,2)=-gaussValue;
    gauss8(4,3)=1;

    gauss8(5,0)=gaussValue;
    gauss8(5,1)=gaussValue;
    gauss8(5,2)=-gaussValue;
    gauss8(5,3)=1;

    gauss8(6,0)=gaussValue;
    gauss8(6,1)=gaussValue;
    gauss8(6,2)=gaussValue;
    gauss8(6,3)=1;

    gauss8(7,0)=-gaussValue;
    gauss8(7,1)=gaussValue;
    gauss8(7,2)=gaussValue;
    gauss8(7,3)=1;
}

void Consolidation3D::createGauss9()
{
    gauss9=MatrixXd::Zero(9,4);
    double gauss1=sqrt(3.0/5.0);
    gauss9<<(1.0/6.0),(1.0/6.0),-gauss1,(1.0/6.0)*(5.0/9.0),
            (1.0/6.0),(2.0/3.0),-gauss1,(1.0/6.0)*(5.0/9.0),
            (2.0/3.0),(1.0/6.0),-gauss1,(1.0/6.0)*(5.0/9.0),
            (1.0/6.0),(1.0/6.0),0.0,(1.0/6.0)*(8.0/9.0),
            (1.0/6.0),(2.0/3.0),0.0,(1.0/6.0)*(8.0/9.0),
            (2.0/3.0),(1.0/6.0),0.0,(1.0/6.0)*(8.0/9.0),
            (1.0/6.0),(1.0/6.0),gauss1,(1.0/6.0)*(5.0/9.0),
            (1.0/6.0),(2.0/3.0),gauss1,(1.0/6.0)*(5.0/9.0),
            (2.0/3.0),(1.0/6.0),gauss1,(1.0/6.0)*(5.0/9.0);
}

void Consolidation3D::createGauss27()
{
    //3x3x3,for brick 20-nodes
    gauss27=MatrixXd::Zero(27,4);
    double g1=-sqrt(3.0/5.0);
    double g2=0;
    double g3=sqrt(3.0/5.0);
    double w1=5.0/9.0;
    double w2=8.0/9.0;
    double w3=5.0/9.0;

    gauss27<<g1,g1,g1,w1*w1*w1,
            g2,g1,g1,w2*w1*w1,
            g3,g1,g1,w3*w1*w1,
            g3,g2,g1,w3*w2*w1,
            g3,g3,g1,w3*w3*w1,
            g2,g3,g1,w2*w3*w1,
            g1,g3,g1,w1*w3*w1,
            g1,g2,g1,w1*w2*w1,
            g2,g2,g1,w2*w2*w1,
            g1,g1,g2,w1*w1*w2,
            g2,g1,g2,w2*w1*w2,
            g3,g1,g2,w3*w1*w2,
            g3,g2,g2,w3*w2*w2,
            g3,g3,g2,w3*w3*w2,
            g2,g3,g2,w2*w3*w2,
            g1,g3,g2,w1*w3*w2,
            g1,g2,g2,w1*w2*w2,
            g2,g2,g2,w2*w2*w2,
            g1,g1,g3,w1*w1*w3,
            g2,g1,g3,w2*w1*w3,
            g3,g1,g3,w3*w1*w3,
            g3,g2,g3,w3*w2*w3,
            g3,g3,g3,w3*w3*w3,
            g2,g3,g3,w2*w3*w3,
            g1,g3,g3,w1*w3*w3,
            g1,g2,g3,w1*w2*w3,
            g2,g2,g3,w2*w2*w3;
}

void Consolidation3D::createGauss13()
{
    //13 gauss points, for Pyramid elements
    gauss13=MatrixXd::Zero(13,4);
    double g1=(7.0/8.0)*sqrt(35.0/39.0);
    double g2=0.6106396;
    double g3=(1.0/56.0)*sqrt(37043.0/35.0);
    double g4=(-1.0/7.0);
    double g5=(-9.0/28.0);
    double g6=0.524394;
    double g7=(-127.0/153.0);
    double w1=0.515003;
    double w2=0.2571837;
    double w3=2.4740050;
    double w4=0.4195157;
    gauss13<<-g1,-g1,g4,w1,
            g1,-g1,g4,w1,
            g1,g1,g4,w1,
            -g1,g1,g4,w1,
            -g2,0,g5,w2,
            g2,0,g5,w2,
            0,-g2,g5,w2,
            0,g2,g5,w2,
            0,0,g6,w3,
            -g3,-g3,g7,w4,
            g3,-g3,g7,w4,
            g3,g3,g7,w4,
            -g3,g3,g7,w4;

}

void Consolidation3D::pauseSystem()
{
    do {
        cout << '\n' << "Press the Enter key to continue.";
    } while (cin.get() != '\n');
}

void Consolidation3D::compareMatrix(Ref<MatrixXi> matrixA, Ref<MatrixXi> matrixB, Ref<MatrixXi> matrixC)
{
    int cols=matrixA.cols();
    for (int i=0;i<cols;i++)
    {
        if (matrixA(0,i)==matrixB(0,i)){matrixC(0,i)=matrixA(0,i);}
        else if(matrixA(0,i)>matrixB(0,i)){matrixC(0,i)=matrixA(0,i);}
        else{matrixC(0,i)=matrixB(0,i);}
    }

}

void Consolidation3D::calculateNodfat()
{
    tet10p=MatrixXi::Zero(10,4);
    hex20p=MatrixXi::Zero(20,4);
    hex8p=MatrixXi::Zero(8,4);
    pyramid13p=MatrixXi::Zero(13,4);
    prism15p=MatrixXi::Zero(15,4);
    drain1Dp=MatrixXi::Zero(2,4);

    tet10p<<1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0;

    hex20p<<1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0;

    hex8p<<1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1;

    pyramid13p<<1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0;

    prism15p<<1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0,
            1,1,1,0;

    drain1Dp<<0,0,0,1,
            0,0,0,1;

    for (int i=0;i<noe;i++)
    {
        int eleIndex=elementsIndex(i,0);
        int eleNode=elements(eleIndex,11);
        int node;
        if(eleNode==20)
        {
            for (int j=0;j<20;j++)
            {
                MatrixXi matrixA= MatrixXi::Zero(1,4);
                MatrixXi matrixB= MatrixXi::Zero(1,4);
                MatrixXi matrixC= MatrixXi::Zero(1,4);
                if(j<10){node=elements(eleIndex,j+1)-1;}
                else{node=elements(eleIndex+1,j-10+1)-1;}
                matrixA=nodfat.row(node);
                matrixB=hex20p.row(j);
                compareMatrix(matrixA,matrixB,matrixC);
                nodfat.row(node)=matrixC;
            }
        }
        else if(eleNode==10)
        {
            for (int j=0;j<10;j++)
            {
                MatrixXi matrixA= MatrixXi::Zero(1,4);
                MatrixXi matrixB= MatrixXi::Zero(1,4);
                MatrixXi matrixC= MatrixXi::Zero(1,4);
                node=elements(eleIndex,j+1)-1;
                matrixA=nodfat.row(node);
                matrixB=tet10p.row(j);
                compareMatrix(matrixA,matrixB,matrixC);
                nodfat.row(node)=matrixC;
            }
        }
        else if(eleNode==13)
        {
            for (int j=0;j<13;j++)
            {
                MatrixXi matrixA= MatrixXi::Zero(1,4);
                MatrixXi matrixB= MatrixXi::Zero(1,4);
                MatrixXi matrixC= MatrixXi::Zero(1,4);
                if(j<5){node=elements(eleIndex,j+1)-1;}
                else{node=elements(eleIndex+1,j-5+1)-1;}
                matrixA=nodfat.row(node);
                matrixB=pyramid13p.row(j);
                compareMatrix(matrixA,matrixB,matrixC);
                nodfat.row(node)=matrixC;
            }
        }
        else if(eleNode==15)
        {
            for (int j=0;j<15;j++)
            {
                MatrixXi matrixA= MatrixXi::Zero(1,4);
                MatrixXi matrixB= MatrixXi::Zero(1,4);
                MatrixXi matrixC= MatrixXi::Zero(1,4);
                if(j<6){node=elements(eleIndex,j+1)-1;}
                else{node=elements(eleIndex+1,j-6+1)-1;}
                matrixA=nodfat.row(node);
                matrixB=prism15p.row(j);
                compareMatrix(matrixA,matrixB,matrixC);
                nodfat.row(node)=matrixC;
            }
        }
        else if(eleNode==2)
        {
            for (int j=0;j<2;j++)
            {
                MatrixXi matrixA= MatrixXi::Zero(1,4);
                MatrixXi matrixB= MatrixXi::Zero(1,4);
                MatrixXi matrixC= MatrixXi::Zero(1,4);
                node=elements(eleIndex,j+1)-1;
                matrixA=nodfat.row(node);
                matrixB=drain1Dp.row(j);
                compareMatrix(matrixA,matrixB,matrixC);
                nodfat.row(node)=matrixC;
            }
        }
    }
}

void Consolidation3D::getPore()
{
    if(pizo.rows()>=1)
    {
        int count=0;
        double avgPore=0;
        for (int j=0;j<pizo.rows();j++)
        {
            int index=pizo(j,0)-1;
            if(nodfmt(index,1)==4)
            {
                count=count++;
                avgPore=avgPore+X(nodfmt(index,0)+3,step);
            }
        }
        poreModel(step,0)=double(avgPore/count);
    }

    else
    {
        return;
    }
}

int Consolidation3D::countPoreBoundary(Ref<MatrixXd> matrixA, Ref<MatrixXi> nodfmt)
{
    int count=0;
    for (int j=0;j<matrixA.rows();j++)

    {
        int index=matrixA(j,0)-1;
        if(nodfmt(index,1)==4)
        {
            count++;
        }
    }
    return count;
}

void Consolidation3D::createNewPoreBoundary(Ref<MatrixXd> matrixOld, Ref<MatrixXd> matrixNew, Ref<MatrixXi> nodfmt)
{
    int count=0;
    for(int j=0;j<matrixOld.rows();j++)
    {
        int index=matrixOld(j,0)-1;
        if(nodfmt(index,1)==4)
        {
            matrixNew(count,0)=nodfmt(index,0)+3;
            matrixNew(count,1)=matrixOld(j,1);
            count++;
        }
    }
}

void Consolidation3D::createGlobalBoundary(Ref<MatrixXd> LocalMatrix, Ref<MatrixXd> GlobalMatrix)
{
    for (int j=0;j<LocalMatrix.rows();j++)
    {
        int equaNum=LocalMatrix(j,0);
        GlobalMatrix(equaNum,0)=equaNum;
        GlobalMatrix(equaNum,1)=1;
        GlobalMatrix(equaNum,2)=LocalMatrix(j,1);
    }
}

double Consolidation3D::calculationError()
{
    error=0;
    int countStep=0;
    for (int j=int(numberOfStep/3);j<numberOfStep;j++)
    {
        countStep=countStep+1;
        double modelResult=poreModel(j,0);
        double analyticalResult=analytical(j,0);
        cout<<"Model results, analytical result: "<<modelResult<<" "<<analyticalResult<<endl;
        error=error+100*(modelResult-analyticalResult)/analyticalResult;
    }
    error=error/countStep;
    return error;
    cout<<"Error: "<<error<<endl;
}

void Consolidation3D::setCd(double CdValue)
{
    this->Cd=CdValue;
}

void Consolidation3D::tet10pMatrix(int &eleNum)
{
    int ii=elementsIndex(eleNum,0); //shorten name
    //Initilize variable
    gauss=MatrixXd::Zero(4,4);
    gauss=gauss4;
    MatrixXi nodeIndex=MatrixXi::Zero(10,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(10,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(10,1); //X coordinates
    MatrixXd Z_coor=MatrixXd::Zero(10,1); //X coordinates
    MatrixXd u0l=MatrixXd::Zero(10,1); //initial x-displacment of elements
    MatrixXd v0l=MatrixXd::Zero(10,1); //initial y-displacment of elements
    MatrixXd w0l=MatrixXd::Zero(10,1); //initial z-displcament of elements
    MatrixXd p0l=MatrixXd::Zero(4,1); //initial pore-pressure of elements

    MatrixXi index_u=MatrixXi::Zero(10,1);
    MatrixXi index_v=MatrixXi::Zero(10,1);
    MatrixXi index_w=MatrixXi::Zero(10,1);
    MatrixXi index_p=MatrixXi::Zero(4,1);
    MatrixXi index_total=MatrixXi::Zero(34,1);

    //element matrices
    MatrixXd Al = MatrixXd::Zero(10,10);
    MatrixXd Bl = MatrixXd::Zero(10,10);
    MatrixXd Cl = MatrixXd::Zero(10,10);
    MatrixXd Dl = MatrixXd::Zero(10,4);

    MatrixXd El = MatrixXd::Zero(10,10);
    MatrixXd Fl = MatrixXd::Zero(10,10);
    MatrixXd Gl = MatrixXd::Zero(10,10);
    MatrixXd Hl = MatrixXd::Zero(10,4);

    MatrixXd Il = MatrixXd::Zero(10,10);
    MatrixXd Jl = MatrixXd::Zero(10,10);
    MatrixXd Ll = MatrixXd::Zero(10,10);
    MatrixXd Ml = MatrixXd::Zero(10,4);

    MatrixXd Nl = MatrixXd::Zero(4,10);
    MatrixXd Ol = MatrixXd::Zero(4,10);
    MatrixXd Pl = MatrixXd::Zero(4,10);
    MatrixXd Ql1 = MatrixXd::Zero(4,4);
    MatrixXd Ql2 = MatrixXd::Zero(4,4);
    MatrixXd Ql = MatrixXd::Zero(4,4);
    MatrixXd QQl = MatrixXd::Zero(4,1);
    MatrixXd Kl=MatrixXd::Zero(34,34);
    //Shape function matrices
    double g1,g2,g3,g4,wi; //gauss-point and weight point
    //For displacment field
    MatrixXd N= MatrixXd::Zero(1,10);
    MatrixXd dNL1= MatrixXd::Zero(1,10);
    MatrixXd dNL2= MatrixXd::Zero(1,10);
    MatrixXd dNL3= MatrixXd::Zero(1,10);
    MatrixXd jacobi=MatrixXd::Zero(3,3);
    MatrixXd jacobiLeft=MatrixXd::Zero(3,10);
    MatrixXd jacobiRight=MatrixXd::Zero(10,3);

    MatrixXd dN=MatrixXd::Zero(3,10);
    MatrixXd dNx=MatrixXd::Zero(1,10);
    MatrixXd dNy=MatrixXd::Zero(1,10);
    MatrixXd dNz=MatrixXd::Zero(1,10);

    //for pore pressure field
    MatrixXd Npore= MatrixXd::Zero(1,4);
    MatrixXd dNporeL1= MatrixXd::Zero(1,4);
    MatrixXd dNporeL2= MatrixXd::Zero(1,4);
    MatrixXd dNporeL3= MatrixXd::Zero(1,4);
    MatrixXd jacobipore=MatrixXd::Zero(3,3);
    MatrixXd jacobiporeLeft=MatrixXd::Zero(3,4);
    MatrixXd jacobiporeRight=MatrixXd::Zero(4,3);

    MatrixXd dNpore=MatrixXd::Zero(3,4);
    MatrixXd dNporex=MatrixXd::Zero(1,4);
    MatrixXd dNporey=MatrixXd::Zero(1,4);
    MatrixXd dNporez=MatrixXd::Zero(1,4);

    //Set all matrix to zero
    u0l.setZero(); v0l.setZero();w0l.setZero();p0l.setZero();
    Al.setZero(); Bl.setZero(); Cl.setZero(); Dl.setZero();
    El.setZero(); Fl.setZero(); Gl.setZero(); Hl.setZero();
    Il.setZero(); Jl.setZero(); Ll.setZero(); Ml.setZero();
    Nl.setZero(); Ol.setZero(); Pl.setZero(); Ql.setZero();    Ql1.setZero(); Ql2.setZero();
    Kl.setZero();

    //get node coordinates
    for (int j=0;j<10;j++)
    {
        nodeIndex(j,0)=elements(ii,j+1)-1;
        index_u(j,0)=nodfmt(nodeIndex(j,0),0)+0;
        index_v(j,0)=nodfmt(nodeIndex(j,0),0)+1;
        index_w(j,0)=nodfmt(nodeIndex(j,0),0)+2;

        if(j<4){index_p(j,0)=nodfmt(nodeIndex(j,0),0)+3;}

        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);
        Z_coor(j,0)=coordinates(nodeIndex(j,0),3);

        u0l(j,0)=X0(index_u(j,0),0);
        v0l(j,0)=X0(index_v(j,0),0);
        w0l(j,0)=X0(index_w(j,0),0);

        if(j<4){p0l(j,0)=X0(index_p(j,0),0);}
    }
    index_total<<index_u,index_v,index_w,index_p;

    //Calculate volume of elements
    MatrixXd volumematrix=MatrixXd::Zero(4,4);
    volumematrix<<1,1,1,1,
            X_coor(0,0),X_coor(1,0),X_coor(2,0),X_coor(3,0),
            Y_coor(0,0),Y_coor(1,0),Y_coor(2,0),Y_coor(3,0),
            Z_coor(0,0),Z_coor(1,0),Z_coor(2,0),Z_coor(3,0);

    double vol;
    vol=volumematrix.determinant();
    vol=vol/6;

    if (vol<0)
    {
        vol=-vol;
    }
    double detj=6*vol;

    //Loop over gauss point
    for (int jj=0;jj<gauss.rows();jj++)
    {
        g1=gauss(jj,0); g2=gauss(jj,1); g3=gauss(jj,2); g4=1-g1-g2-g3; wi=gauss(jj,3);
        N<<(2*g1-1)*g1,(2*g2-1)*g2,(2*g3-1)*g3,(2*g4-1)*g4,4*g1*g2,4*g2*g3,4*g1*g3,4*g1*g4,4*g2*g4,4*g3*g4;
        dNL1<<4*g1-1,0,0,1-4*g4,4*g2,0,4*g3,4*(g4-g1),-4*g2,-4*g3;
        dNL2<<0,4*g2-1,0,1-4*g4,4*g1,4*g3,0,-4*g1,4*(g4-g2),-4*g3;
        dNL3<<0,0,4*g3-1,1-4*g4,0,4*g2,4*g1,-4*g1,-4*g2,4*(g4-g3);
        jacobiLeft<<dNL1,dNL2,dNL3;
        jacobiRight<<X_coor,Y_coor,Z_coor;
        jacobi=jacobiLeft*jacobiRight;
        dN=jacobi.inverse()*jacobiLeft;
        dNx=dN.row(0);
        dNy=dN.row(1);
        dNz=dN.row(2);

        Al=Al+wi*detj*((Ke+4*Ge/3)*dNx.transpose()*dNx+Ge*dNy.transpose()*dNy+Ge*dNz.transpose()*dNz);
        Bl=Bl+wi*detj*((Ke-2*Ge/3)*dNx.transpose()*dNy+Ge*dNy.transpose()*dNx);
        Cl=Cl+wi*detj*((Ke-2*Ge/3)*dNx.transpose()*dNz+Ge*dNz.transpose()*dNx);

        Fl=Fl+wi*detj*((Ke+4*Ge/3)*dNy.transpose()*dNy+Ge*dNx.transpose()*dNx+Ge*dNz.transpose()*dNz);
        El=El+wi*detj*((Ke-2*Ge/3)*dNy.transpose()*dNx+Ge*dNx.transpose()*dNy);
        Gl=Gl+wi*detj*((Ke-2*Ge/3)*dNy.transpose()*dNz+Ge*dNz.transpose()*dNy);

        Ll=Ll+wi*detj*((Ke+4*Ge/3)*dNz.transpose()*dNz+Ge*dNx.transpose()*dNx+Ge*dNy.transpose()*dNy);
        Il=Il+wi*detj*((Ke-2*Ge/3)*dNz.transpose()*dNx+Ge*dNx.transpose()*dNz);
        Jl=Jl+wi*detj*((Ke-2*Ge/3)*dNz.transpose()*dNy+Ge*dNy.transpose()*dNz);

        Npore<<g1,g2,g3,g4;
        dNporeL1<<1,0,0,-1;
        dNporeL2<<0,1,0,-1;
        dNporeL3<<0,0,1,-1;
        jacobiporeLeft<<dNporeL1,dNporeL2,dNporeL3;
        jacobiporeRight<<X_coor(0,0),Y_coor(0,0),Z_coor(0,0),
                X_coor(1,0),Y_coor(1,0),Z_coor(1,0),
                X_coor(2,0),Y_coor(2,0),Z_coor(2,0),
                X_coor(3,0),Y_coor(3,0),Z_coor(3,0);
        jacobipore=jacobiporeLeft*jacobiporeRight;
        dNpore=jacobipore.inverse()*jacobiporeLeft;
        dNporex=dNpore.row(0);
        dNporey=dNpore.row(1);
        dNporez=dNpore.row(2);

        Dl=Dl-A*detj*wi*dNx.transpose()*Npore;
        Hl=Hl-A*detj*wi*dNy.transpose()*Npore;
        Ml=Ml-A*detj*wi*dNz.transpose()*Npore;
        Ql1=Ql1+Se*detj*wi*Npore.transpose()*Npore;
        Ql2=Ql2+detj*wi*((khe/gf)*dNporex.transpose()*dNporex+(kve/gf)*dNporey.transpose()*dNporey+(khe/gf)*dNporez.transpose()*dNporez);
    }

    Nl=-1*Dl.transpose();
    Ol=-1*Hl.transpose();
    Pl=-1*Ml.transpose();
    Ql=Ql1+timeInterval*Ql2;
    QQl=Nl*u0l+Ol*v0l+Pl*w0l+Ql1*p0l;

    Kl<<Al,Bl,Cl,Dl,
            El,Fl,Gl,Hl,
            Il,Jl,Ll,Ml,
            Nl,Ol,Pl,Ql;

    //Transfer QQl to F
    for (int j=0;j<4;j++)
    {
        F(index_p(j,0),0)=F(index_p(j,0),0)+QQl(j,0);
    }

    for (int j=0;j<Kl.rows();j=j+1)
    {
        int jj=DirichletAll(index_total(j,0),1);
        if(jj==1)
        {
            for (int k1=0;k1<Kl.cols();k1=k1+1)
            {
                int kk=index_total(k1,0);
                double valueBoundary=DirichletAll(index_total(j,0),2);
                F(kk,0)=F(kk,0)-Kl(k1,j)*valueBoundary;
            }
            Kl.row(j).setZero();
            Kl.col(j).setZero();
            Kl(j,j)=1;
        }
    }

    //Push back Kl to trip_total
    for (int j=0;j<Kl.rows();j++)
    {
        for (int jj=0;jj<Kl.cols();jj++)
        {
            int row_i=index_total(j,0);
            int col_j=index_total(jj,0);
            double val_ij=Kl(j,jj);
            trip_total.push_back(Trip(row_i,col_j,val_ij));
        }
    }

}

void Consolidation3D::hex20pMatrix(int &eleNum)
{
    int eleIndex=elementsIndex(eleNum,0);
    int dod=20; //Degree freedoms of displacement
    int dop=8; //Degree freedoms of pore pressure
    gauss=MatrixXd::Zero(27,4);
    gauss=gauss27;
    MatrixXi nodeIndex=MatrixXi::Zero(dod,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd Z_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd u0l=MatrixXd::Zero(dod,1); //initial x-displacment of elements
    MatrixXd v0l=MatrixXd::Zero(dod,1); //initial y-displacment of elements
    MatrixXd w0l=MatrixXd::Zero(dod,1); //initial z-displcament of elements
    MatrixXd p0l=MatrixXd::Zero(dop,1); //initial pore-pressure of elements

    MatrixXi index_u=MatrixXi::Zero(dod,1);
    MatrixXi index_v=MatrixXi::Zero(dod,1);
    MatrixXi index_w=MatrixXi::Zero(dod,1);
    MatrixXi index_p=MatrixXi::Zero(dop,1);
    MatrixXi index_total=MatrixXi::Zero(3*dod+dop,1);

    //element matrices
    MatrixXd Al = MatrixXd::Zero(dod,dod);
    MatrixXd Bl = MatrixXd::Zero(dod,dod);
    MatrixXd Cl = MatrixXd::Zero(dod,dod);
    MatrixXd Dl = MatrixXd::Zero(dod,dop);

    MatrixXd El = MatrixXd::Zero(dod,dod);
    MatrixXd Fl = MatrixXd::Zero(dod,dod);
    MatrixXd Gl = MatrixXd::Zero(dod,dod);
    MatrixXd Hl = MatrixXd::Zero(dod,dop);

    MatrixXd Il = MatrixXd::Zero(dod,dod);
    MatrixXd Jl = MatrixXd::Zero(dod,dod);
    MatrixXd Ll = MatrixXd::Zero(dod,dod);
    MatrixXd Ml = MatrixXd::Zero(dod,dop);

    MatrixXd Nl = MatrixXd::Zero(dop,dod);
    MatrixXd Ol = MatrixXd::Zero(dop,dod);
    MatrixXd Pl = MatrixXd::Zero(dop,dod);
    MatrixXd Ql1 = MatrixXd::Zero(dop,dop);
    MatrixXd Ql2 = MatrixXd::Zero(dop,dop);
    MatrixXd Ql = MatrixXd::Zero(dop,dop);
    MatrixXd QQl = MatrixXd::Zero(dop,1);
    MatrixXd Kl=MatrixXd::Zero(3*dod+dop,3*dod+dop);

    //Shape function matrices
    double g1,g2,g3,g4,wi; //gauss-point and weight point
    //For displacment field
    MatrixXd N= MatrixXd::Zero(1,dod);
    MatrixXd dNL1= MatrixXd::Zero(1,dod);
    MatrixXd dNL2= MatrixXd::Zero(1,dod);
    MatrixXd dNL3= MatrixXd::Zero(1,dod);
    MatrixXd jacobi=MatrixXd::Zero(3,3);
    MatrixXd jacobiLeft=MatrixXd::Zero(3,dod);
    MatrixXd jacobiRight=MatrixXd::Zero(dod,3);

    MatrixXd dN=MatrixXd::Zero(3,dod);
    MatrixXd dNx=MatrixXd::Zero(1,dod);
    MatrixXd dNy=MatrixXd::Zero(1,dod);
    MatrixXd dNz=MatrixXd::Zero(1,dod);

    //for pore pressure field
    MatrixXd Npore= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL1= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL2= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL3= MatrixXd::Zero(1,dop);
    MatrixXd jacobipore=MatrixXd::Zero(3,3);
    MatrixXd jacobiporeLeft=MatrixXd::Zero(3,dop);
    MatrixXd jacobiporeRight=MatrixXd::Zero(dop,3);

    MatrixXd dNpore=MatrixXd::Zero(3,dop);
    MatrixXd dNporex=MatrixXd::Zero(1,dop);
    MatrixXd dNporey=MatrixXd::Zero(1,dop);
    MatrixXd dNporez=MatrixXd::Zero(1,dop);

    //Set all matrix to zero
    u0l.setZero(); v0l.setZero();w0l.setZero();p0l.setZero();
    Al.setZero(); Bl.setZero(); Cl.setZero(); Dl.setZero();
    El.setZero(); Fl.setZero(); Gl.setZero(); Hl.setZero();
    Il.setZero(); Jl.setZero(); Ll.setZero(); Ml.setZero();
    Nl.setZero(); Ol.setZero(); Pl.setZero(); Ql.setZero();    Ql1.setZero(); Ql2.setZero();
    Kl.setZero();

    //get node coordinates
    for (int j=0;j<dod;j++)
    {
        if(j<10){nodeIndex(j,0)=elements(eleIndex,j+1)-1;}
        else{nodeIndex(j,0)=elements(eleIndex+1,j-10+1)-1;}
        index_u(j,0)=nodfmt(nodeIndex(j,0),0)+0;
        index_v(j,0)=nodfmt(nodeIndex(j,0),0)+1;
        index_w(j,0)=nodfmt(nodeIndex(j,0),0)+2;
        if(j<dop){index_p(j,0)=nodfmt(nodeIndex(j,0),0)+3;}

        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);
        Z_coor(j,0)=coordinates(nodeIndex(j,0),3);

        u0l(j,0)=X0(index_u(j,0),0);
        v0l(j,0)=X0(index_v(j,0),0);
        w0l(j,0)=X0(index_w(j,0),0);
        if(j<dop){p0l(j,0)=X0(index_p(j,0),0);}
    }
    index_total<<index_u,index_v,index_w,index_p;

    //Loop over gauss point
    for (int jj=0;jj<gauss.rows();jj++)
    {
        g1=gauss(jj,0); g2=gauss(jj,1); g3=gauss(jj,2); wi=gauss(jj,3);
        //20 nodes  shape function
        N(0,0)=(1.0/8.0)*(1-g1)*(1-g2)*(1-g3)*(-g1-g2-g3-2);
        N(0,1)=(1.0/8.0)*(1-g1)*(1-g2)*(1+g3)*(-g1-g2+g3-2);
        N(0,2)=(1.0/8.0)*(1+g1)*(1-g2)*(1+g3)*(+g1-g2+g3-2);
        N(0,3)=(1.0/8.0)*(1+g1)*(1-g2)*(1-g3)*(+g1-g2-g3-2);
        N(0,4)=(1.0/8.0)*(1-g1)*(1+g2)*(1-g3)*(-g1+g2-g3-2);
        N(0,5)=(1.0/8.0)*(1-g1)*(1+g2)*(1+g3)*(-g1+g2+g3-2);
        N(0,6)=(1.0/8.0)*(1+g1)*(1+g2)*(1+g3)*(+g1+g2+g3-2);
        N(0,7)=(1.0/8.0)*(1+g1)*(1+g2)*(1-g3)*(+g1+g2-g3-2);
        N(0,8)=(1.0/4.0)*(1-g1)*(1-g2)*(1-g3*g3);
        N(0,9)=(1.0/4.0)*(1-g1*g1)*(1-g2)*(1+g3);
        N(0,10)=(1.0/4.0)*(1+g1)*(1-g2)*(1-g3*g3);
        N(0,11)=(1.0/4.0)*(1-g1*g1)*(1-g2)*(1-g3);
        N(0,12)=(1.0/4.0)*(1-g1)*(1+g2)*(1-g3*g3);
        N(0,13)=(1.0/4.0)*(1-g1*g1)*(1+g2)*(1+g3);
        N(0,14)=(1.0/4.0)*(1+g1)*(1+g2)*(1-g3*g3);
        N(0,15)=(1.0/4.0)*(1-g1*g1)*(1+g2)*(1-g3);
        N(0,16)=(1.0/4.0)*(1-g1)*(1-g2*g2)*(1-g3);
        N(0,17)=(1.0/4.0)*(1-g1)*(1-g2*g2)*(1+g3);
        N(0,18)=(1.0/4.0)*(1+g1)*(1-g2*g2)*(1+g3);
        N(0,19)=(1.0/4.0)*(1+g1)*(1-g2*g2)*(1-g3);

        dNL1(0,0)=(1.0/8.0)*(1-g2)*(1-g3)*(-1*(-g1-g2-g3-2)-1*(1-g1));
        dNL1(0,1)=(1.0/8.0)*(1-g2)*(1+g3)*(-1*(-g1-g2+g3-2)-1*(1-g1));
        dNL1(0,2)=(1.0/8.0)*(1-g2)*(1+g3)*(+1*(+g1-g2+g3-2)+1*(1+g1));
        dNL1(0,3)=(1.0/8.0)*(1-g2)*(1-g3)*(+1*(+g1-g2-g3-2)+1*(1+g1));
        dNL1(0,4)=(1.0/8.0)*(1+g2)*(1-g3)*(-1*(-g1+g2-g3-2)-1*(1-g1));
        dNL1(0,5)=(1.0/8.0)*(1+g2)*(1+g3)*(-1*(-g1+g2+g3-2)-1*(1-g1));
        dNL1(0,6)=(1.0/8.0)*(1+g2)*(1+g3)*(+1*(+g1+g2+g3-2)+1*(1+g1));
        dNL1(0,7)=(1.0/8.0)*(1+g2)*(1-g3)*(+1*(+g1+g2-g3-2)+1*(1+g1));
        dNL1(0,8)=(-1.0/4.0)*(1-g2)*(1-g3*g3);
        dNL1(0,9)=(-g1/2)*(1-g2)*(1+g3);
        dNL1(0,10)=(+1.0/4.0)*(1-g2)*(1-g3*g3);
        dNL1(0,11)=(-g1/2)*(1-g2)*(1-g3);
        dNL1(0,12)=(-1.0/4.0)*(1+g2)*(1-g3*g3);
        dNL1(0,13)=(-g1/2.0)*(1+g2)*(1+g3);
        dNL1(0,14)=(+1.0/4.0)*(1+g2)*(1-g3*g3);
        dNL1(0,15)=(-g1/2.0)*(1+g2)*(1-g3);
        dNL1(0,16)=(-1.0/4.0)*(1-g2*g2)*(1-g3);
        dNL1(0,17)=(-1.0/4.0)*(1-g2*g2)*(1+g3);
        dNL1(0,18)=(+1.0/4.0)*(1-g2*g2)*(1+g3);
        dNL1(0,19)=(+1.0/4.0)*(1-g2*g2)*(1-g3);

        dNL2(0,0)=(1.0/8.0)*(1-g1)*(1-g3)*(-1*(-g1-g2-g3-2)-1*(1-g2));
        dNL2(0,1)=(1.0/8.0)*(1-g1)*(1+g3)*(-1*(-g1-g2+g3-2)-1*(1-g2));
        dNL2(0,2)=(1.0/8.0)*(1+g1)*(1+g3)*(-1*(+g1-g2+g3-2)-1*(1-g2));
        dNL2(0,3)=(1.0/8.0)*(1+g1)*(1-g3)*(-1*(+g1-g2-g3-2)-1*(1-g2));
        dNL2(0,4)=(1.0/8.0)*(1-g1)*(1-g3)*(+1*(-g1+g2-g3-2)+1*(1+g2));
        dNL2(0,5)=(1.0/8.0)*(1-g1)*(1+g3)*(+1*(-g1+g2+g3-2)+1*(1+g2));
        dNL2(0,6)=(1.0/8.0)*(1+g1)*(1+g3)*(+1*(+g1+g2+g3-2)+1*(1+g2));
        dNL2(0,7)=(1.0/8.0)*(1+g1)*(1-g3)*(+1*(+g1+g2-g3-2)+1*(1+g2));
        dNL2(0,8)=(-1.0/4.0)*(1-g1)*(1-g3*g3);
        dNL2(0,9)=(-1.0/4.0)*(1-g1*g1)*(1+g3);
        dNL2(0,10)=(-1.0/4.0)*(1+g1)*(1-g3*g3);
        dNL2(0,11)=(-1.0/4.0)*(1-g1*g1)*(1-g3);
        dNL2(0,12)=(1.0/4.0)*(1-g1)*(1-g3*g3);
        dNL2(0,13)=(1.0/4.0)*(1-g1*g1)*(1+g3);
        dNL2(0,14)=(1.0/4.0)*(1+g1)*(1-g3*g3);
        dNL2(0,15)=(1.0/4.0)*(1-g1*g1)*(1-g3);
        dNL2(0,16)=(-g2/2)*(1-g1)*(1-g3);
        dNL2(0,17)=(-g2/2)*(1-g1)*(1+g3);
        dNL2(0,18)=(-g2/2)*(1+g1)*(1+g3);
        dNL2(0,19)=(-g2/2)*(1+g1)*(1-g3);

        dNL3(0,0)=(1.0/8.0)*(1-g1)*(1-g2)*(-1*(-g1-g2-g3-2)-1*(1-g3));
        dNL3(0,1)=(1.0/8.0)*(1-g1)*(1-g2)*(1*(-g1-g2+g3-2)+1*(1+g3));
        dNL3(0,2)=(1.0/8.0)*(1+g1)*(1-g2)*(1*(g1-g2+g3-2)+1*(1+g3));
        dNL3(0,3)=(1.0/8.0)*(1+g1)*(1-g2)*(-1*(g1-g2-g3-2)-1*(1-g3));
        dNL3(0,4)=(1.0/8.0)*(1-g1)*(1+g2)*(-1*(-g1+g2-g3-2)-1*(1-g3));
        dNL3(0,5)=(1.0/8.0)*(1-g1)*(1+g2)*(1*(-g1+g2+g3-2)+1*(1+g3));
        dNL3(0,6)=(1.0/8.0)*(1+g1)*(1+g2)*(1*(g1+g2+g3-2)+1*(1+g3));
        dNL3(0,7)=(1.0/8.0)*(1+g1)*(1+g2)*(-1*(g1+g2-g3-2)-1*(1-g3));
        dNL3(0,8)=(-g3/2)*(1-g1)*(1-g2);
        dNL3(0,9)=(1.0/4.0)*(1-g1*g1)*(1-g2);
        dNL3(0,10)=(-g3/2)*(1+g1)*(1-g2);
        dNL3(0,11)=(-1.0/4.0)*(1-g1*g1)*(1-g2);
        dNL3(0,12)=(-g3/2)*(1-g1)*(1+g2);
        dNL3(0,13)=(1.0/4.0)*(1-g1*g1)*(1+g2);
        dNL3(0,14)=(-g3/2)*(1+g1)*(1+g2);
        dNL3(0,15)=(-1.0/4.0)*(1-g1*g1)*(1+g2);
        dNL3(0,16)=(-1.0/4.0)*(1-g1)*(1-g2*g2);
        dNL3(0,17)=(1.0/4.0)*(1-g1)*(1-g2*g2);
        dNL3(0,18)=(1.0/4.0)*(1+g1)*(1-g2*g2);
        dNL3(0,19)=(-1.0/4.0)*(1+g1)*(1-g2*g2);

        jacobiLeft<<dNL1,dNL2,dNL3;
        jacobiRight<<X_coor,Y_coor,Z_coor;
        jacobi=jacobiLeft*jacobiRight;
        double detj;
        detj=abs(jacobi.determinant());
        dN=jacobi.inverse()*jacobiLeft;
        dNx=dN.row(0);
        dNy=dN.row(1);
        dNz=dN.row(2);
        //--------------
        //8 nodes shape function
        Npore(0,0)=(1.0/8.0)*(1-g1)*(1-g2)*(1-g3);
        Npore(0,1)=(1.0/8.0)*(1-g1)*(1-g2)*(1+g3);
        Npore(0,2)=(1.0/8.0)*(1+g1)*(1-g2)*(1+g3);
        Npore(0,3)=(1.0/8.0)*(1+g1)*(1-g2)*(1-g3);
        Npore(0,4)=(1.0/8.0)*(1-g1)*(1+g2)*(1-g3);
        Npore(0,5)=(1.0/8.0)*(1-g1)*(1+g2)*(1+g3);
        Npore(0,6)=(1.0/8.0)*(1+g1)*(1+g2)*(1+g3);
        Npore(0,7)=(1.0/8.0)*(1+g1)*(1+g2)*(1-g3);

        dNporeL1(0,0)=(-1.0/8.0)*(1-g2)*(1-g3);
        dNporeL1(0,1)=(-1.0/8.0)*(1-g2)*(1+g3);
        dNporeL1(0,2)=(+1.0/8.0)*(1-g2)*(1+g3);
        dNporeL1(0,3)=(+1.0/8.0)*(1-g2)*(1-g3);
        dNporeL1(0,4)=(-1.0/8.0)*(1+g2)*(1-g3);
        dNporeL1(0,5)=(-1.0/8.0)*(1+g2)*(1+g3);
        dNporeL1(0,6)=(+1.0/8.0)*(1+g2)*(1+g3);
        dNporeL1(0,7)=(+1.0/8.0)*(1+g2)*(1-g3);

        dNporeL2(0,0)=(-1.0/8.0)*(1-g1)*(1-g3);
        dNporeL2(0,1)=(-1.0/8.0)*(1-g1)*(1+g3);
        dNporeL2(0,2)=(-1.0/8.0)*(1+g1)*(1+g3);
        dNporeL2(0,3)=(-1.0/8.0)*(1+g1)*(1-g3);
        dNporeL2(0,4)=(+1.0/8.0)*(1-g1)*(1-g3);
        dNporeL2(0,5)=(+1.0/8.0)*(1-g1)*(1+g3);
        dNporeL2(0,6)=(+1.0/8.0)*(1+g1)*(1+g3);
        dNporeL2(0,7)=(+1.0/8.0)*(1+g1)*(1-g3);

        dNporeL3(0,0)=(-1.0/8.0)*(1-g1)*(1-g2);
        dNporeL3(0,1)=(+1.0/8.0)*(1-g1)*(1-g2);
        dNporeL3(0,2)=(+1.0/8.0)*(1+g1)*(1-g2);
        dNporeL3(0,3)=(-1.0/8.0)*(1+g1)*(1-g2);
        dNporeL3(0,4)=(-1.0/8.0)*(1-g1)*(1+g2);
        dNporeL3(0,5)=(+1.0/8.0)*(1-g1)*(1+g2);
        dNporeL3(0,6)=(+1.0/8.0)*(1+g1)*(1+g2);
        dNporeL3(0,7)=(-1.0/8.0)*(1+g1)*(1+g2);
        jacobiporeLeft<<dNporeL1,dNporeL2,dNporeL3;
        jacobiporeRight<<X_coor(0,0),Y_coor(0,0),Z_coor(0,0),
                X_coor(1,0),Y_coor(1,0),Z_coor(1,0),
                X_coor(2,0),Y_coor(2,0),Z_coor(2,0),
                X_coor(3,0),Y_coor(3,0),Z_coor(3,0),
                X_coor(4,0),Y_coor(4,0),Z_coor(4,0),
                X_coor(5,0),Y_coor(5,0),Z_coor(5,0),
                X_coor(6,0),Y_coor(6,0),Z_coor(6,0),
                X_coor(7,0),Y_coor(7,0),Z_coor(7,0);

        jacobipore=jacobiporeLeft*jacobiporeRight;
        dNpore=jacobipore.inverse()*jacobiporeLeft;
        dNporex=dNpore.row(0);
        dNporey=dNpore.row(1);
        dNporez=dNpore.row(2);
        double detjpore;
        detjpore=abs(jacobipore.determinant());

        //Element matrix
        Al=Al+wi*detj*((Ke+4*Ge/3)*dNx.transpose()*dNx+Ge*dNy.transpose()*dNy+Ge*dNz.transpose()*dNz);
        Bl=Bl+wi*detj*((Ke-2*Ge/3)*dNx.transpose()*dNy+Ge*dNy.transpose()*dNx);
        Cl=Cl+wi*detj*((Ke-2*Ge/3)*dNx.transpose()*dNz+Ge*dNz.transpose()*dNx);

        Fl=Fl+wi*detj*((Ke+4*Ge/3)*dNy.transpose()*dNy+Ge*dNx.transpose()*dNx+Ge*dNz.transpose()*dNz);
        El=El+wi*detj*((Ke-2*Ge/3)*dNy.transpose()*dNx+Ge*dNx.transpose()*dNy);
        Gl=Gl+wi*detj*((Ke-2*Ge/3)*dNy.transpose()*dNz+Ge*dNz.transpose()*dNy);

        Ll=Ll+wi*detj*((Ke+4*Ge/3)*dNz.transpose()*dNz+Ge*dNx.transpose()*dNx+Ge*dNy.transpose()*dNy);
        Il=Il+wi*detj*((Ke-2*Ge/3)*dNz.transpose()*dNx+Ge*dNx.transpose()*dNz);
        Jl=Jl+wi*detj*((Ke-2*Ge/3)*dNz.transpose()*dNy+Ge*dNy.transpose()*dNz);

        Dl=Dl-A*detj*wi*dNx.transpose()*Npore;
        Hl=Hl-A*detj*wi*dNy.transpose()*Npore;
        Ml=Ml-A*detj*wi*dNz.transpose()*Npore;
        Ql1=Ql1+Se*detj*wi*Npore.transpose()*Npore;
        Ql2=Ql2+detj*wi*((khe/gf)*dNporex.transpose()*dNporex+(kve/gf)*dNporey.transpose()*dNporey+(khe/gf)*dNporez.transpose()*dNporez);
    }

    Nl=-1*Dl.transpose();
    Ol=-1*Hl.transpose();
    Pl=-1*Ml.transpose();
    Ql=Ql1+timeInterval*Ql2;
    QQl=Nl*u0l+Ol*v0l+Pl*w0l+Ql1*p0l;

    Kl<<Al,Bl,Cl,Dl,
            El,Fl,Gl,Hl,
            Il,Jl,Ll,Ml,
            Nl,Ol,Pl,Ql;
    //Transfer QQl to F
    for (int j=0;j<dop;j++)
    {
        F(index_p(j,0),0)=F(index_p(j,0),0)+QQl(j,0);
    }

    //Assign boundary condition
    for (int j=0;j<Kl.rows();j=j+1)
    {
        int jj=DirichletAll(index_total(j,0),1);
        if(jj==1)
        {
            for (int k1=0;k1<Kl.cols();k1=k1+1)
            {
                int kk=index_total(k1,0);
                double valueBoundary=DirichletAll(index_total(j,0),2);
                F(kk,0)=F(kk,0)-Kl(k1,j)*valueBoundary;
            }
            Kl.row(j).setZero();
            Kl.col(j).setZero();
            Kl(j,j)=1;
        }
    }

    //Push back Kl to trip_total
    for (int j=0;j<Kl.rows();j++)
    {
        for (int jj=0;jj<Kl.cols();jj++)
        {
            int row_i=index_total(j,0);
            int col_j=index_total(jj,0);
            double val_ij=Kl(j,jj);
            trip_total.push_back(Trip(row_i,col_j,val_ij));
        }
    }

}

void Consolidation3D::pyra13pMatrix(int &eleNum)
{
    int eleIndex=elementsIndex(eleNum,0);
    int dod=13; //Degree freedoms of displacement
    int dop=5; //Degree freedoms of pore pressure
    gauss=MatrixXd::Zero(27,4);
    gauss=gauss27;
    MatrixXi nodeIndex=MatrixXi::Zero(dod,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd Z_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd u0l=MatrixXd::Zero(dod,1); //initial x-displacment of elements
    MatrixXd v0l=MatrixXd::Zero(dod,1); //initial y-displacment of elements
    MatrixXd w0l=MatrixXd::Zero(dod,1); //initial z-displcament of elements
    MatrixXd p0l=MatrixXd::Zero(dop,1); //initial pore-pressure of elements

    MatrixXi index_u=MatrixXi::Zero(dod,1);
    MatrixXi index_v=MatrixXi::Zero(dod,1);
    MatrixXi index_w=MatrixXi::Zero(dod,1);
    MatrixXi index_p=MatrixXi::Zero(dop,1);
    MatrixXi index_total=MatrixXi::Zero(3*dod+dop,1);

    //element matrices
    MatrixXd Al = MatrixXd::Zero(dod,dod);
    MatrixXd Bl = MatrixXd::Zero(dod,dod);
    MatrixXd Cl = MatrixXd::Zero(dod,dod);
    MatrixXd Dl = MatrixXd::Zero(dod,dop);

    MatrixXd El = MatrixXd::Zero(dod,dod);
    MatrixXd Fl = MatrixXd::Zero(dod,dod);
    MatrixXd Gl = MatrixXd::Zero(dod,dod);
    MatrixXd Hl = MatrixXd::Zero(dod,dop);

    MatrixXd Il = MatrixXd::Zero(dod,dod);
    MatrixXd Jl = MatrixXd::Zero(dod,dod);
    MatrixXd Ll = MatrixXd::Zero(dod,dod);
    MatrixXd Ml = MatrixXd::Zero(dod,dop);

    MatrixXd Nl = MatrixXd::Zero(dop,dod);
    MatrixXd Ol = MatrixXd::Zero(dop,dod);
    MatrixXd Pl = MatrixXd::Zero(dop,dod);
    MatrixXd Ql1 = MatrixXd::Zero(dop,dop);
    MatrixXd Ql2 = MatrixXd::Zero(dop,dop);
    MatrixXd Ql = MatrixXd::Zero(dop,dop);
    MatrixXd QQl = MatrixXd::Zero(dop,1);
    MatrixXd Kl=MatrixXd::Zero(3*dod+dop,3*dod+dop);

    //Shape function matrices
    double g1,g2,g3,wi; //gauss-point and weight point
    //For displacment field
    MatrixXd N= MatrixXd::Zero(1,dod);
    MatrixXd dNL1= MatrixXd::Zero(1,dod);
    MatrixXd dNL2= MatrixXd::Zero(1,dod);
    MatrixXd dNL3= MatrixXd::Zero(1,dod);
    MatrixXd jacobi=MatrixXd::Zero(3,3);
    MatrixXd jacobiLeft=MatrixXd::Zero(3,dod);
    MatrixXd jacobiRight=MatrixXd::Zero(dod,3);

    MatrixXd dN=MatrixXd::Zero(3,dod);
    MatrixXd dNx=MatrixXd::Zero(1,dod);
    MatrixXd dNy=MatrixXd::Zero(1,dod);
    MatrixXd dNz=MatrixXd::Zero(1,dod);

    //for pore pressure field
    MatrixXd Npore= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL1= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL2= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL3= MatrixXd::Zero(1,dop);
    MatrixXd jacobipore=MatrixXd::Zero(3,3);
    MatrixXd jacobiporeLeft=MatrixXd::Zero(3,dop);
    MatrixXd jacobiporeRight=MatrixXd::Zero(dop,3);

    MatrixXd dNpore=MatrixXd::Zero(3,dop);
    MatrixXd dNporex=MatrixXd::Zero(1,dop);
    MatrixXd dNporey=MatrixXd::Zero(1,dop);
    MatrixXd dNporez=MatrixXd::Zero(1,dop);

    //Set all matrix to zero
    u0l.setZero(); v0l.setZero();w0l.setZero();p0l.setZero();
    Al.setZero(); Bl.setZero(); Cl.setZero(); Dl.setZero();
    El.setZero(); Fl.setZero(); Gl.setZero(); Hl.setZero();
    Il.setZero(); Jl.setZero(); Ll.setZero(); Ml.setZero();
    Nl.setZero(); Ol.setZero(); Pl.setZero(); Ql.setZero();    Ql1.setZero(); Ql2.setZero();
    Kl.setZero();

    //get node coordinates
    for (int j=0;j<dod;j++)
    {
        if(j<5){nodeIndex(j,0)=elements(eleIndex,j+1)-1;}
        else{nodeIndex(j,0)=elements(eleIndex+1,j-5+1)-1;}
        index_u(j,0)=nodfmt(nodeIndex(j,0),0)+0;
        index_v(j,0)=nodfmt(nodeIndex(j,0),0)+1;
        index_w(j,0)=nodfmt(nodeIndex(j,0),0)+2;
        if(j<dop){index_p(j,0)=nodfmt(nodeIndex(j,0
                                                ),0)+3;}

        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);
        Z_coor(j,0)=coordinates(nodeIndex(j,0),3);

        u0l(j,0)=X0(index_u(j,0),0);
        v0l(j,0)=X0(index_v(j,0),0);
        w0l(j,0)=X0(index_w(j,0),0);
        if(j<dop){p0l(j,0)=X0(index_p(j,0),0);}
    }
    index_total<<index_u,index_v,index_w,index_p;
    //Loop over gauss point
    for (int jj=0;jj<gauss.rows();jj++)
    {
        g1=gauss(jj,0); g2=gauss(jj,1); g3=gauss(jj,2); wi=gauss(jj,3);
        //20 nodes  shape function
        N(0,0)=(1.0/16.0)*(1-g1)*(1-g2)*(1-g3)*(-2-g1+g1*g3-g2+g2*g3);
        N(0,1)=(1.0/16.0)*(1+g1)*(1-g2)*(1-g3)*(-2+g1-g1*g3-g2+g2*g3);
        N(0,2)=(1.0/16.0)*(1+g1)*(1+g2)*(1-g3)*(-2+g1-g1*g3+g2-g2*g3);
        N(0,3)=(1.0/16.0)*(1-g1)*(1+g2)*(1-g3)*(-2-g1+g1*g3+g2-g2*g3);
        N(0,4)=(1.0/2.0)*(1+g3)*g3;
        N(0,5)=(1.0/8.0)*(1-g1*g1)*(1-g2)*(1-2*g3+g3*g3);
        N(0,6)=(1.0/8.0)*(1+g1)*(1-g2*g2)*(1-2*g3+g3*g3);
        N(0,7)=(1.0/8.0)*(1-g1*g1)*(1+g2)*(1-2*g3+g3*g3);
        N(0,8)=(1.0/8.0)*(1-g1)*(1-g2*g2)*(1-2*g3+g3*g3);
        N(0,9)=(1.0/4.0)*(1-g1-g2+g1*g2)*(1-g3*g3);
        N(0,10)=(1.0/4.0)*(1+g1-g2-g1*g2)*(1-g3*g3);
        N(0,11)=(1.0/4.0)*(1+g1+g2+g1*g2)*(1-g3*g3);
        N(0,12)=(1.0/4.0)*(1-g1+g2-g1*g2)*(1-g3*g3);

        dNL1(0,0)=(1.0/16.0)*(1-g2)*(1-g3)*(-1*(-2-g1+g1*g3-g2+g2*g3)+(1-g1)*(-1+g3));
        dNL1(0,1)=(1.0/16.0)*(1-g2)*(1-g3)*(+1*(-2+g1-g1*g3-g2+g2*g3)+(1+g1)*(+1-g3));
        dNL1(0,2)=(1.0/16.0)*(1+g2)*(1-g3)*(+1*(-2+g1-g1*g3+g2-g2*g3)+(1+g1)*(+1-g3));
        dNL1(0,3)=(1.0/16.0)*(1+g2)*(1-g3)*(-1*(-2-g1+g1*g3+g2-g2*g3)+(1-g1)*(-1+g3));
        dNL1(0,4)=0;
        dNL1(0,5)=(-1.0/4.0)*g1*(1-g2)*(1-2*g3+g3*g3);
        dNL1(0,6)=(+1.0/8.0)*(1-g2*g2)*(1-2*g3+g3*g3);
        dNL1(0,7)=(-1.0/4.0)*g1*(1+g2)*(1-2*g3+g3*g3);
        dNL1(0,8)=(-1.0/8.0)*(1-g2*g2)*(1-2*g3+g3*g3);
        dNL1(0,9)=(1.0/4.0)*(1-g3*g3)*(-1+g2);
        dNL1(0,10)=(1.0/4.0)*(1-g3*g3)*(+1-g2);
        dNL1(0,11)=(1.0/4.0)*(1-g3*g3)*(+1+g2);
        dNL1(0,12)=(1.0/4.0)*(1-g3*g3)*(-1-g2);

        dNL2(0,0)=(1.0/16.0)*(1-g1)*(1-g3)*(-1*(-2-g1+g1*g3-g2+g2*g3)+(1-g2)*(-1+g3));
        dNL2(0,1)=(1.0/16.0)*(1+g1)*(1-g3)*(-1*(-2+g1-g1*g3-g2+g2*g3)+(1-g2)*(-1+g3));
        dNL2(0,2)=(1.0/16.0)*(1+g1)*(1-g3)*(+1*(-2+g1-g1*g3+g2-g2*g3)+(1+g2)*(+1-g3));
        dNL2(0,3)=(1.0/16.0)*(1-g1)*(1-g3)*(+1*(-2-g1+g1*g3+g2-g2*g3)+(1+g2)*(+1-g3));
        dNL2(0,4)=0;
        dNL2(0,5)=(-1.0/8.0)*(1-g1*g1)*(1-2*g3+g3*g3);
        dNL2(0,6)=(-1.0/4.0)*g2*(1+g1)*(1-2*g3+g3*g3);
        dNL2(0,7)=(+1.0/8.0)*(1-g1*g1)*(1-2*g3+g3*g3);
        dNL2(0,8)=(-1.0/4.0)*g2*(1-g1)*(1-2*g3+g3*g3);
        dNL2(0,9)=(1.0/4.0)*(-1+g1)*(1-g3*g3);
        dNL2(0,10)=(1.0/4.0)*(-1-g1)*(1-g3*g3);
        dNL2(0,11)=(1.0/4.0)*(+1+g1)*(1-g3*g3);
        dNL2(0,12)=(1.0/4.0)*(+1-g1)*(1-g3*g3);

        dNL3(0,0)=(1.0/16.0)*(1-g1)*(1-g2)*(-1*(-2-g1+g1*g3-g2+g2*g3)+(1-g3)*(+g1+g2));
        dNL3(0,1)=(1.0/16.0)*(1+g1)*(1-g2)*(-1*(-2+g1-g1*g3-g2+g2*g3)+(1-g3)*(-g1+g2));
        dNL3(0,2)=(1.0/16.0)*(1+g1)*(1+g2)*(-1*(-2+g1-g1*g3+g2-g2*g3)+(1-g3)*(-g1-g2));
        dNL3(0,3)=(1.0/16.0)*(1-g1)*(1+g2)*(-1*(-2-g1+g1*g3+g2-g2*g3)+(1-g3)*(+g1-g2));
        dNL3(0,4)=(1.0+2*g3)/2.0;
        dNL3(0,5)=(1.0/4.0)*(1-g1*g1)*(1-g2)*(-1+g3);
        dNL3(0,6)=(1.0/4.0)*(1+g1)*(1-g2*g2)*(-1+g3);
        dNL3(0,7)=(1.0/4.0)*(1-g1*g1)*(1+g2)*(-1+g3);
        dNL3(0,8)=(1.0/4.0)*(1-g1)*(1-g2*g2)*(-1+g3);
        dNL3(0,9)=(-g3/2.0)*(1-g1-g2+g1*g2);
        dNL3(0,10)=(-g3/2.0)*(1+g1-g2-g1*g2);
        dNL3(0,11)=(-g3/2.0)*(1+g1+g2+g1*g2);
        dNL3(0,12)=(-g3/2.0)*(1-g1+g2-g1*g2);

        jacobiLeft<<dNL1,dNL2,dNL3;
        jacobiRight<<X_coor,Y_coor,Z_coor;
        jacobi=jacobiLeft*jacobiRight;
        double detj;
        detj=abs(jacobi.determinant());

        dN=jacobi.inverse()*jacobiLeft;
        dNx=dN.row(0);
        dNy=dN.row(1);
        dNz=dN.row(2);
        //--------------
        //8 nodes shape function
        Npore(0,0)=(1.0/8.0)*(1-g1)*(1-g2)*(1-g3);
        Npore(0,1)=(1.0/8.0)*(1+g1)*(1-g2)*(1-g3);
        Npore(0,2)=(1.0/8.0)*(1+g1)*(1+g2)*(1-g3);
        Npore(0,3)=(1.0/8.0)*(1-g1)*(1+g2)*(1-g3);
        Npore(0,4)=(1.0/2.0)*(1+g3);


        dNporeL1(0,0)=(-1.0/8.0)*(1-g2)*(1-g3);
        dNporeL1(0,1)=(+1.0/8.0)*(1-g2)*(1-g3);
        dNporeL1(0,2)=(+1.0/8.0)*(1+g2)*(1-g3);
        dNporeL1(0,3)=(-1.0/8.0)*(1+g2)*(1-g3);
        dNporeL1(0,4)=0;


        dNporeL2(0,0)=(-1.0/8.0)*(1-g1)*(1-g3);
        dNporeL2(0,1)=(-1.0/8.0)*(1+g1)*(1-g3);
        dNporeL2(0,2)=(+1.0/8.0)*(1+g1)*(1-g3);
        dNporeL2(0,3)=(+1.0/8.0)*(1-g1)*(1-g3);
        dNporeL2(0,4)=0;


        dNporeL3(0,0)=(-1.0/8.0)*(1-g1)*(1-g2);
        dNporeL3(0,1)=(-1.0/8.0)*(1+g1)*(1-g2);
        dNporeL3(0,2)=(-1.0/8.0)*(1+g1)*(1+g2);
        dNporeL3(0,3)=(-1.0/8.0)*(1-g1)*(1+g2);
        dNporeL3(0,4)=1.0/2.0;

        jacobiporeLeft<<dNporeL1,dNporeL2,dNporeL3;
        jacobiporeRight<<X_coor(0,0),Y_coor(0,0),Z_coor(0,0),
                X_coor(1,0),Y_coor(1,0),Z_coor(1,0),
                X_coor(2,0),Y_coor(2,0),Z_coor(2,0),
                X_coor(3,0),Y_coor(3,0),Z_coor(3,0),
                X_coor(4,0),Y_coor(4,0),Z_coor(4,0);

        jacobipore=jacobiporeLeft*jacobiporeRight;
        dNpore=jacobipore.inverse()*jacobiporeLeft;
        dNporex=dNpore.row(0);
        dNporey=dNpore.row(1);
        dNporez=dNpore.row(2);
        double detjpore;
        detjpore=abs(jacobipore.determinant());

        //Element matrix
        Al=Al+wi*detj*((Ke+4*Ge/3)*dNx.transpose()*dNx+Ge*dNy.transpose()*dNy+Ge*dNz.transpose()*dNz);
        Bl=Bl+wi*detj*((Ke-2*Ge/3)*dNx.transpose()*dNy+Ge*dNy.transpose()*dNx);
        Cl=Cl+wi*detj*((Ke-2*Ge/3)*dNx.transpose()*dNz+Ge*dNz.transpose()*dNx);

        Fl=Fl+wi*detj*((Ke+4*Ge/3)*dNy.transpose()*dNy+Ge*dNx.transpose()*dNx+Ge*dNz.transpose()*dNz);
        El=El+wi*detj*((Ke-2*Ge/3)*dNy.transpose()*dNx+Ge*dNx.transpose()*dNy);
        Gl=Gl+wi*detj*((Ke-2*Ge/3)*dNy.transpose()*dNz+Ge*dNz.transpose()*dNy);

        Ll=Ll+wi*detj*((Ke+4*Ge/3)*dNz.transpose()*dNz+Ge*dNx.transpose()*dNx+Ge*dNy.transpose()*dNy);
        Il=Il+wi*detj*((Ke-2*Ge/3)*dNz.transpose()*dNx+Ge*dNx.transpose()*dNz);
        Jl=Jl+wi*detj*((Ke-2*Ge/3)*dNz.transpose()*dNy+Ge*dNy.transpose()*dNz);

        Dl=Dl-A*detj*wi*dNx.transpose()*Npore;
        Hl=Hl-A*detj*wi*dNy.transpose()*Npore;
        Ml=Ml-A*detj*wi*dNz.transpose()*Npore;
        Ql1=Ql1+Se*detj*wi*Npore.transpose()*Npore;
        Ql2=Ql2+detj*wi*((khe/gf)*dNporex.transpose()*dNporex+(kve/gf)*dNporey.transpose()*dNporey+(khe/gf)*dNporez.transpose()*dNporez);
    }

    Nl=-1*Dl.transpose();
    Ol=-1*Hl.transpose();
    Pl=-1*Ml.transpose();
    Ql=Ql1+timeInterval*Ql2;
    QQl=Nl*u0l+Ol*v0l+Pl*w0l+Ql1*p0l;

    Kl<<Al,Bl,Cl,Dl,
            El,Fl,Gl,Hl,
            Il,Jl,Ll,Ml,
            Nl,Ol,Pl,Ql;
    //Transfer QQl to F
    for (int j=0;j<dop;j++)
    {
        F(index_p(j,0),0)=F(index_p(j,0),0)+QQl(j,0);
    }

    //Assign boundary condition
    for (int j=0;j<Kl.rows();j=j+1)
    {
        int jj=DirichletAll(index_total(j,0),1);
        if(jj==1)
        {
            for (int k1=0;k1<Kl.cols();k1=k1+1)
            {
                int kk=index_total(k1,0);
                double valueBoundary=DirichletAll(index_total(j,0),2);
                F(kk,0)=F(kk,0)-Kl(k1,j)*valueBoundary; //Wrong but nobody knows
            }
            Kl.row(j).setZero();
            Kl.col(j).setZero();
            Kl(j,j)=1;
        }
    }

    //Push back Kl to trip_total
    for (int j=0;j<Kl.rows();j++)
    {
        for (int jj=0;jj<Kl.cols();jj++)
        {
            int row_i=index_total(j,0);
            int col_j=index_total(jj,0);
            double val_ij=Kl(j,jj);
            trip_total.push_back(Trip(row_i,col_j,val_ij));
        }
    }

}

void Consolidation3D::prism15pMatrix(int &eleNum)
{
    int eleIndex=elementsIndex(eleNum,0);
    int dod=15; //Degree freedoms of displacement
    int dop=6; //Degree freedoms of pore pressure
    gauss=MatrixXd::Zero(9,4);
    gauss=gauss9;
    MatrixXi nodeIndex=MatrixXi::Zero(dod,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd Z_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd u0l=MatrixXd::Zero(dod,1); //initial x-displacment of elements
    MatrixXd v0l=MatrixXd::Zero(dod,1); //initial y-displacment of elements
    MatrixXd w0l=MatrixXd::Zero(dod,1); //initial z-displcament of elements
    MatrixXd p0l=MatrixXd::Zero(dop,1); //initial pore-pressure of elements

    MatrixXi index_u=MatrixXi::Zero(dod,1);
    MatrixXi index_v=MatrixXi::Zero(dod,1);
    MatrixXi index_w=MatrixXi::Zero(dod,1);
    MatrixXi index_p=MatrixXi::Zero(dop,1);
    MatrixXi index_total=MatrixXi::Zero(3*dod+dop,1);

    //element matrices
    MatrixXd Al = MatrixXd::Zero(dod,dod);
    MatrixXd Bl = MatrixXd::Zero(dod,dod);
    MatrixXd Cl = MatrixXd::Zero(dod,dod);
    MatrixXd Dl = MatrixXd::Zero(dod,dop);

    MatrixXd El = MatrixXd::Zero(dod,dod);
    MatrixXd Fl = MatrixXd::Zero(dod,dod);
    MatrixXd Gl = MatrixXd::Zero(dod,dod);
    MatrixXd Hl = MatrixXd::Zero(dod,dop);

    MatrixXd Il = MatrixXd::Zero(dod,dod);
    MatrixXd Jl = MatrixXd::Zero(dod,dod);
    MatrixXd Ll = MatrixXd::Zero(dod,dod);
    MatrixXd Ml = MatrixXd::Zero(dod,dop);

    MatrixXd Nl = MatrixXd::Zero(dop,dod);
    MatrixXd Ol = MatrixXd::Zero(dop,dod);
    MatrixXd Pl = MatrixXd::Zero(dop,dod);
    MatrixXd Ql1 = MatrixXd::Zero(dop,dop);
    MatrixXd Ql2 = MatrixXd::Zero(dop,dop);
    MatrixXd Ql = MatrixXd::Zero(dop,dop);
    MatrixXd QQl = MatrixXd::Zero(dop,1);
    MatrixXd Kl=MatrixXd::Zero(3*dod+dop,3*dod+dop);

    //Shape function matrices
    double g1,g2,g3,g4,wi; //gauss-point and weight point
    //For displacment field
    MatrixXd N= MatrixXd::Zero(1,dod);
    MatrixXd dNL1= MatrixXd::Zero(1,dod);
    MatrixXd dNL2= MatrixXd::Zero(1,dod);
    MatrixXd dNL3= MatrixXd::Zero(1,dod);
    MatrixXd jacobi=MatrixXd::Zero(3,3);
    MatrixXd jacobiLeft=MatrixXd::Zero(3,dod);
    MatrixXd jacobiRight=MatrixXd::Zero(dod,3);

    MatrixXd dN=MatrixXd::Zero(3,dod);
    MatrixXd dNx=MatrixXd::Zero(1,dod);
    MatrixXd dNy=MatrixXd::Zero(1,dod);
    MatrixXd dNz=MatrixXd::Zero(1,dod);

    //for pore pressure field
    MatrixXd Npore= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL1= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL2= MatrixXd::Zero(1,dop);
    MatrixXd dNporeL3= MatrixXd::Zero(1,dop);
    MatrixXd jacobipore=MatrixXd::Zero(3,3);
    MatrixXd jacobiporeLeft=MatrixXd::Zero(3,dop);
    MatrixXd jacobiporeRight=MatrixXd::Zero(dop,3);

    MatrixXd dNpore=MatrixXd::Zero(3,dop);
    MatrixXd dNporex=MatrixXd::Zero(1,dop);
    MatrixXd dNporey=MatrixXd::Zero(1,dop);
    MatrixXd dNporez=MatrixXd::Zero(1,dop);

    //Set all matrix to zero
    u0l.setZero(); v0l.setZero();w0l.setZero();p0l.setZero();
    Al.setZero(); Bl.setZero(); Cl.setZero(); Dl.setZero();
    El.setZero(); Fl.setZero(); Gl.setZero(); Hl.setZero();
    Il.setZero(); Jl.setZero(); Ll.setZero(); Ml.setZero();
    Nl.setZero(); Ol.setZero(); Pl.setZero(); Ql.setZero();    Ql1.setZero(); Ql2.setZero();
    Kl.setZero();

    //get node coordinates
    for (int j=0;j<dod;j++)
    {
        if(j<6){nodeIndex(j,0)=elements(eleIndex,j+1)-1;}
        else{nodeIndex(j,0)=elements(eleIndex+1,j-6+1)-1;}
        index_u(j,0)=nodfmt(nodeIndex(j,0),0)+0;
        index_v(j,0)=nodfmt(nodeIndex(j,0),0)+1;
        index_w(j,0)=nodfmt(nodeIndex(j,0),0)+2;
        if(j<dop){index_p(j,0)=nodfmt(nodeIndex(j,0),0)+3;}

        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);
        Z_coor(j,0)=coordinates(nodeIndex(j,0),3);

        u0l(j,0)=X0(index_u(j,0),0);
        v0l(j,0)=X0(index_v(j,0),0);
        w0l(j,0)=X0(index_w(j,0),0);
        if(j<dop){p0l(j,0)=X0(index_p(j,0),0);}
    }
    index_total<<index_u,index_v,index_w,index_p;

    //Loop over gauss point
    for (int jj=0;jj<gauss.rows();jj++)
    {
        g1=gauss(jj,0); g2=gauss(jj,1); g3=gauss(jj,2); wi=gauss(jj,3);
        g4=1.0-g1-g2;
        //20 nodes  shape function
        N(0,0)=0.5*(g1*(2*g1-1.0)*(1.0-g3)-g1*(1.0-g3*g3));
        N(0,1)=0.5*(g2*(2*g2-1.0)*(1.0-g3)-g2*(1.0-g3*g3));
        N(0,2)=0.5*(g4*(2*g4-1.0)*(1.0-g3)-g4*(1.0-g3*g3));

        N(0,3)=0.5*(g1*(2*g1-1.0)*(1.0+g3)-g1*(1.0-g3*g3));
        N(0,4)=0.5*(g2*(2*g2-1.0)*(1.0+g3)-g2*(1.0-g3*g3));
        N(0,5)=0.5*(g4*(2*g4-1.0)*(1.0+g3)-g4*(1.0-g3*g3));

        N(0,6)=2.0*g1*g2*(1-g3);
        N(0,7)=2.0*g2*g4*(1-g3);
        N(0,8)=2.0*g1*g4*(1-g3);

        N(0,9)=2.0*g1*g2*(1+g3);
        N(0,10)=2.0*g2*g4*(1+g3);
        N(0,11)=2.0*g1*g4*(1+g3);

        N(0,12)=g1*(1.0-g3*g3);
        N(0,13)=g2*(1.0-g3*g3);
        N(0,14)=g4*(1.0-g3*g3);

        dNL1(0,0)=0.5*((4.0*g1-1.0)*(1.0-g3)-(1.0-g3*g3));
        dNL1(0,1)=0.0;
        dNL1(0,2)=0.5*((4.0*g1+4.0*g2-3.0)*(1.0-g3)+(1.0-g3*g3));

        dNL1(0,3)=0.5*((4.0*g1-1.0)*(1.0+g3)-(1.0-g3*g3));
        dNL1(0,4)=0.0;
        dNL1(0,5)=0.5*((4.0*g1+4.0*g2-3.0)*(1.0+g3)+(1.0-g3*g3));

        dNL1(0,6)=2.0*g2*(1.0-g3);
        dNL1(0,7)=-2.0*g2*(1.0-g3);
        dNL1(0,8)=2.0*(1.0-2.0*g1-g2)*(1.0-g3);

        dNL1(0,9)=2.0*g2*(1.0+g3);
        dNL1(0,10)=-2.0*g2*(1.0+g3);
        dNL1(0,11)=2.0*(1.0-2.0*g1-g2)*(1.0+g3);

        dNL1(0,12)=1.0-g3*g3;
        dNL1(0,13)=0.0;
        dNL1(0,14)=-(1.0-g3*g3);

        dNL2(0,0)=0.0;
        dNL2(0,1)=0.5*((4.0*g2-1.0)*(1.0-g3)-(1.0-g3*g3));
        dNL2(0,2)=0.5*((4.0*g2+4.0*g1-3.0)*(1.0-g3)+(1.0-g3*g3));

        dNL2(0,3)=0.0;
        dNL2(0,4)=0.5*((4.0*g2-1.0)*(1.0+g3)-(1.0-g3*g3));
        dNL2(0,5)=0.5*((4.0*g2+4.0*g1-3.0)*(1.0+g3)+(1.0-g3*g3));

        dNL2(0,6)=2.0*g1*(1.0-g3);
        dNL2(0,7)=2.0*(1.0-g1-2.0*g2)*(1.0-g3);
        dNL2(0,8)=-2.0*g1*(1.0-g3);

        dNL2(0,9)=2.0*g1*(1.0+g3);
        dNL2(0,10)=2.0*(1.0-g1-2*g2)*(1.0+g3);
        dNL2(0,11)=-2.0*g1*(1.0+g3);

        dNL2(0,12)=0.0;
        dNL2(0,13)=(1.0-g3*g3);
        dNL2(0,14)=-(1.0-g3*g3);


        dNL3(0,0)=0.5*(-2.0*g1*(2*g1-1.0)+2.0*g3*g1);
        dNL3(0,1)=0.5*(-2.0*g2*(2*g2-1.0)+2.0*g3*g2);
        dNL3(0,2)=0.5*(-2.0*g4*(2*g4-1.0)+2.0*g3*g4);

        dNL3(0,3)=0.5*(+2.0*g1*(2*g1-1.0)+2.0*g3*g1);
        dNL3(0,4)=0.5*(+2.0*g2*(2*g2-1.0)+2.0*g3*g2);
        dNL3(0,5)=0.5*(+2.0*g4*(2*g4-1.0)+2.0*g3*g4);

        dNL3(0,6)=-2.0*g1*g2;
        dNL3(0,7)=-2.0*g2*g4;
        dNL3(0,8)=-2.0*g1*g4;

        dNL3(0,9)=2.0*g1*g2;
        dNL3(0,10)=2.0*g2*g4;
        dNL3(0,11)=2.0*g1*g4;

        dNL3(0,12)=-2.0*g3*g1;
        dNL3(0,13)=-2.0*g3*g2;
        dNL3(0,14)=-2.0*g3*g4;


        jacobiLeft<<dNL1,dNL2,dNL3;
        jacobiRight<<X_coor,Y_coor,Z_coor;
        jacobi=jacobiLeft*jacobiRight;
        double detj;
        detj=abs(jacobi.determinant());
        //        cout<<"Detj 13 nodes: "<<detj<<endl;
        dN=jacobi.inverse()*jacobiLeft;
        dNx=dN.row(0);
        dNy=dN.row(1);
        dNz=dN.row(2);
        //--------------
        //6 nodes shape function
        Npore(0,0)=0.5*g1*(1.0-g3);
        Npore(0,1)=0.5*g2*(1.0-g3);
        Npore(0,2)=0.5*g4*(1.0-g3);
        Npore(0,3)=0.5*g1*(1.0+g3);
        Npore(0,4)=0.5*g2*(1.0+g3);
        Npore(0,5)=0.5*g4*(1.0+g3);

        dNporeL1(0,0)=0.5*(1.0-g3);
        dNporeL1(0,1)=0.0;
        dNporeL1(0,2)=-0.5*(1.0-g3);
        dNporeL1(0,3)=0.5*(1.0+g3);
        dNporeL1(0,4)=0.0;
        dNporeL1(0,5)=-0.5*(1.0+g3);

        dNporeL2(0,0)=0;
        dNporeL2(0,1)=0.5*(1.0-g3);
        dNporeL2(0,2)=-0.5*(1.0-g3);
        dNporeL2(0,3)=0.0;
        dNporeL2(0,4)=0.5*(1.0+g3);
        dNporeL2(0,5)=-0.5*(1.0+g3);

        dNporeL3(0,0)=-0.5*g1;
        dNporeL3(0,1)=-0.5*g2;
        dNporeL3(0,2)=-0.5*g4;
        dNporeL3(0,3)=0.5*g1;
        dNporeL3(0,4)=0.5*g2;
        dNporeL3(0,5)=0.5*g4;

        jacobiporeLeft<<dNporeL1,dNporeL2,dNporeL3;
        jacobiporeRight<<X_coor(0,0),Y_coor(0,0),Z_coor(0,0),
                X_coor(1,0),Y_coor(1,0),Z_coor(1,0),
                X_coor(2,0),Y_coor(2,0),Z_coor(2,0),
                X_coor(3,0),Y_coor(3,0),Z_coor(3,0),
                X_coor(4,0),Y_coor(4,0),Z_coor(4,0),
                X_coor(5,0),Y_coor(5,0),Z_coor(5,0);

        jacobipore=jacobiporeLeft*jacobiporeRight;
        dNpore=jacobipore.inverse()*jacobiporeLeft;
        dNporex=dNpore.row(0);
        dNporey=dNpore.row(1);
        dNporez=dNpore.row(2);
        double detjpore;
        detjpore=abs(jacobipore.determinant());
        //        cout<<"det j 6 nodes: "<<detjpore<<endl;

        //Element matrix
        Al=Al+wi*detj*((Ke+4*Ge/3)*dNx.transpose()*dNx+Ge*dNy.transpose()*dNy+Ge*dNz.transpose()*dNz);
        Bl=Bl+wi*detj*((Ke-2*Ge/3)*dNx.transpose()*dNy+Ge*dNy.transpose()*dNx);
        Cl=Cl+wi*detj*((Ke-2*Ge/3)*dNx.transpose()*dNz+Ge*dNz.transpose()*dNx);

        Fl=Fl+wi*detj*((Ke+4*Ge/3)*dNy.transpose()*dNy+Ge*dNx.transpose()*dNx+Ge*dNz.transpose()*dNz);
        El=El+wi*detj*((Ke-2*Ge/3)*dNy.transpose()*dNx+Ge*dNx.transpose()*dNy);
        Gl=Gl+wi*detj*((Ke-2*Ge/3)*dNy.transpose()*dNz+Ge*dNz.transpose()*dNy);

        Ll=Ll+wi*detj*((Ke+4*Ge/3)*dNz.transpose()*dNz+Ge*dNx.transpose()*dNx+Ge*dNy.transpose()*dNy);
        Il=Il+wi*detj*((Ke-2*Ge/3)*dNz.transpose()*dNx+Ge*dNx.transpose()*dNz);
        Jl=Jl+wi*detj*((Ke-2*Ge/3)*dNz.transpose()*dNy+Ge*dNy.transpose()*dNz);

        Dl=Dl-A*detj*wi*dNx.transpose()*Npore;
        Hl=Hl-A*detj*wi*dNy.transpose()*Npore;
        Ml=Ml-A*detj*wi*dNz.transpose()*Npore;
        Ql1=Ql1+Se*detj*wi*Npore.transpose()*Npore;
        Ql2=Ql2+detj*wi*((khe/gf)*dNporex.transpose()*dNporex+(kve/gf)*dNporey.transpose()*dNporey+(khe/gf)*dNporez.transpose()*dNporez);
    }

    Nl=-1*Dl.transpose();
    Ol=-1*Hl.transpose();
    Pl=-1*Ml.transpose();
    Ql=Ql1+timeInterval*Ql2;
    QQl=Nl*u0l+Ol*v0l+Pl*w0l+Ql1*p0l;

    Kl<<Al,Bl,Cl,Dl,
            El,Fl,Gl,Hl,
            Il,Jl,Ll,Ml,
            Nl,Ol,Pl,Ql;

    //Transfer QQl to F
    for (int j=0;j<dop;j++)
    {
        F(index_p(j,0),0)=F(index_p(j,0),0)+QQl(j,0);
    }

    //Assign boundary condition
    for (int j=0;j<Kl.rows();j=j+1)
    {
        int jj=DirichletAll(index_total(j,0),1);
        if(jj==1)
        {
            for (int k1=0;k1<Kl.cols();k1=k1+1)
            {
                int kk=index_total(k1,0);
                double valueBoundary=DirichletAll(index_total(j,0),2);
                F(kk,0)=F(kk,0)-Kl(k1,j)*valueBoundary;
            }
            Kl.row(j).setZero();
            Kl.col(j).setZero();
            Kl(j,j)=1;
        }
    }

    //Push back Kl to trip_total
    for (int j=0;j<Kl.rows();j++)
    {
        for (int jj=0;jj<Kl.cols();jj++)
        {
            int row_i=index_total(j,0);
            int col_j=index_total(jj,0);
            double val_ij=Kl(j,jj);
            trip_total.push_back(Trip(row_i,col_j,val_ij));
        }
    }
}

void Consolidation3D::drain1DMatrix(int &eleNum)
{
    int eleIndex=elementsIndex(eleNum,0);
    int dod=2; //Degree freedoms of displacement
    int dop=2; //Degree freedoms of pore pressure
    MatrixXi nodeIndex=MatrixXi::Zero(dod,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd Z_coor=MatrixXd::Zero(dod,1); //X coordinates
    MatrixXd p0l=MatrixXd::Zero(dop,1); //initial pore-pressure of elements

    MatrixXi index_p=MatrixXi::Zero(dop,1);
    MatrixXi index_total=MatrixXi::Zero(dop,1);

    //element matrices
    MatrixXd Ql1 = MatrixXd::Zero(dop,dop);
    MatrixXd Ql2 = MatrixXd::Zero(dop,dop);
    MatrixXd Ql = MatrixXd::Zero(dop,dop);
    MatrixXd QQl = MatrixXd::Zero(dop,1);
    MatrixXd Kl=MatrixXd::Zero(3*dod+dop,3*dod+dop);

    //get node coordinates
    for (int j=0;j<dod;j++)
    {
        nodeIndex(j,0)=elements(eleIndex,j+1)-1;
        index_p(j,0)=nodfmt(nodeIndex(j,0),0)+3;
        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);
        Z_coor(j,0)=coordinates(nodeIndex(j,0),3);
        p0l(j,0)=X0(index_p(j,0),0);
    }
    index_total=index_p;
    double dx,dy,dz,length;
    dx=X_coor(1,0)-X_coor(0,0);
    dy=Y_coor(1,0)-Y_coor(0,0);
    dz=Z_coor(1,0)-Z_coor(0,0);
    length=sqrt(dx*dx+dy*dy+dz*dz);
    MatrixXd Ae=MatrixXd::Zero(2,2);

    Ae<<1.0f,-1.0f,
        -1.0f,1.0f;

    double TempValue;
    TempValue=(Aw*kw)/(gf*length);

    Ql2=TempValue*Ae;

    Ql=timeInterval*Ql2;
    QQl.setZero();
    Kl=Ql;

    //Transfer QQl to F
    for (int j=0;j<dop;j++)
    {
        F(index_p(j,0),0)=F(index_p(j,0),0)+QQl(j,0);
    }

    //Assign boundary condition
    for (int j=0;j<Kl.rows();j=j+1)
    {
        int jj=DirichletAll(index_total(j,0),1);
        if(jj==1)
        {
            for (int k1=0;k1<Kl.cols();k1=k1+1)
            {
                int kk=index_total(k1,0);
                double valueBoundary=DirichletAll(index_total(j,0),2);
                F(kk,0)=F(kk,0)-Kl(k1,j)*valueBoundary;
            }
            Kl.row(j).setZero();
            Kl.col(j).setZero();
            Kl(j,j)=1;
        }
    }

    //Push back Kl to trip_total
    for (int j=0;j<Kl.rows();j++)
    {
        for (int jj=0;jj<Kl.cols();jj++)
        {
            int row_i=index_total(j,0);
            int col_j=index_total(jj,0);
            double val_ij=Kl(j,jj);
            trip_total.push_back(Trip(row_i,col_j,val_ij));
        }
    }
}

void Consolidation3D::GetMeshData(meshDataNames meshName)
{
    coordFile=meshName.coordFile;
    eleFile=meshName.eleFile;
    fixXFile=meshName.fixXFile;
    fixYFile=meshName.fixYFile;
    fixZFile=meshName.fixZFile;
    fixHFile=meshName.fixHFile;
    forceYFile=meshName.forceYFile;
    folderName=meshName.folderName;
}

void Consolidation3D::GetSoilData(Soil soilData)
{
    re=soilData.re;
    ksRatio=soilData.ksRatio;
    kv=soilData.kv;
    kh=kv*re;
    ks=kh*ksRatio;
    voidRatio=soilData.voidRatio;
    K=soilData.K;
    poissionRatio=soilData.v;
    Cvr=soilData.Ch;
}

void Consolidation3D::GetPVDData(PVD PVDData)
{
    R=PVDData.unitCellRadius;
    Rs=PVDData.smearZoneRadius;
    Rw=PVDData.PVDRadius;
    kw=PVDData.PVDConductivity;
    H=0.5f*PVDData.PVDLength;
    n=R/Rw;
    s=Rs/Rw;
}

void Consolidation3D::GetAnalysisType(SolutionType typeData)
{
    analysisType=typeData.analysisType;
    Cd0=typeData.Cd0;
    timeInterval0=typeData.timeIncrement;
    p0=typeData.intialPressure;
    tol=typeData.tolerance;
    numberOfStep=typeData.numberOfStep;
}
