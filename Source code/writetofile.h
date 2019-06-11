#ifndef WRITETOFILE_H
#define WRITETOFILE_H
#include<string>
#include<iostream>
#include<ostream>
#include<Eigen/Dense>
#include<QString>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;

class WriteToFile
{
public:
    string fileName;
    void ToFile(const Ref<const MatrixXd>matrixA);
    void ToFile(const Ref<const MatrixXi>matrixA);
    WriteToFile();
};

#endif // WRITETOFILE_H
