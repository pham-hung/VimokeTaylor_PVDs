#include "getfile.h"

void GetFile::SetFullFileName()
{
    fullFileName=folderName+"/"+fileName;
}

void GetFile::get_dimension(QString fullFileName)
{
    row=0;
    col=0;
    QFile File(fullFileName);
    QString Line,Temp;
    QTextStream in(&File);

//    qDebug()<<"Import data from File:"<<fileName<<endl;

    if(!File.open(QFile::ReadOnly|QFile::Text))
    {
        qDebug()<<"Cannot open the file";
    }
    else
    {
        while(!in.atEnd())
        {
            Line=in.readLine();
            row=row+1;
            //qDebug()<<Line<<endl;
            if (row==1)
            {
               QStringList list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
               col=list.count();
            }
        }
    }
    File.close();
    //qDebug()<<"The numer of row is:"<<row<<endl;
    //qDebug()<<"The numer of col is:"<<col<<endl;
}

void GetFile::import_file(QString fullFileName)
{
    //qDebug()<<"Write file to Matrix"<<endl;
    //qDebug()<<"Row and col is:"<<row<<col;
    data_file.resize(row,col);
    double Value=0;
    QFile File(fullFileName);
    QString Line,Temp;
    QTextStream in(&File);

    if(!File.open(QFile::ReadOnly|QFile::Text))
    {
        qDebug()<<"cannot open the file"<<endl;
    }
    else
    {
        for(int i=0;i<row;i++)
        {
            Line=in.readLine();
            QStringList list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
            for (int j=0;j<col;j++)
            {
                Temp=list[j];
                Value=Temp.toDouble();
                data_file(i,j)=Value;
            }
        }

    }

}

void GetFile::DoGetFile()
{
    data_file.setZero();
    data_file.resize(1,1);
    SetFullFileName();
    GetFile::get_dimension(fullFileName);
    GetFile::import_file(fullFileName);
    qDebug()<<"Imported File ...."<<fileName<<endl;
}

void GetFile::DoGetFile(QString fullFileName)
{
    data_file.setZero();
    data_file.resize(1,1);
    GetFile::get_dimension(fullFileName);
    GetFile::import_file(fullFileName);
    qDebug()<<"Imported File ...."<<fullFileName<<endl;
}

void GetFile::DoGetFile(QString folderName, QString fileName)
{
    DoGetFile();
}

GetFile::GetFile()
{
}
