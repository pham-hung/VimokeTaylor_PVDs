#ifndef GETFILE_H
#define GETFILE_H
#include <QFile>
#include <QString>
#include <Eigen/Dense>
#include <QDebug>
#include <QStringList>
#include <QTextStream>
using namespace Eigen;
using namespace std;

class GetFile
{
public:
    int col;
    int row;
    MatrixXd data_file=MatrixXd::Zero(1,1);
    QString fileName;
    QString folderName;
    QString fullFileName;
    void SetFullFileName();
    void get_dimension(QString fullFileName);
    void import_file(QString fullFileName);
    void DoGetFile();
    void DoGetFile(QString fullFileName);
    void DoGetFile(QString folderName, QString fileName);
    GetFile();
};

#endif // GETFILE_H
