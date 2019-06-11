#include "writetofile.h"
#include<iostream>
#include<string>
#include<fstream>
#include<iomanip>
#include<QString>

void WriteToFile::ToFile(const Ref<const MatrixXd> matrixA)
{
    int row=matrixA.rows();
    int col=matrixA.cols();
    double val=0;
    ofstream file_;
    file_.open(fileName);

    for (int i=0;i<row;i=i+1)
    {
        for (int j=0;j<col;j=j+1)
        {
          val=matrixA(i,j);
          file_<<std::setw(13)<<std::setprecision(5)<<std::scientific<<val<<' ';
        }
        file_<<'\n';
    }
    file_.close();
    cout<<"Exported successfully"<<endl;
}

void WriteToFile::ToFile(const Ref<const MatrixXi> matrixA)
{
    int row=matrixA.rows();
    int col=matrixA.cols();
    double val=0;
    ofstream file_;
    file_.open(fileName);

    for (int i=0;i<row;i=i+1)
    {
        for (int j=0;j<col;j=j+1)
        {
          val=matrixA(i,j);
          file_<<val<<' ';
        }
        file_<<'\n';
    }
    file_.close();
    cout<<"Exported successfully"<<endl;
}

WriteToFile::WriteToFile()
{

}
