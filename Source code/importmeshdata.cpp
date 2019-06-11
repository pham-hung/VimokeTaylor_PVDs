#include "importmeshdata.h"
#include "ui_importmeshdata.h"

ImportMeshData::ImportMeshData(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ImportMeshData)
{
    ui->setupUi(this);
}

ImportMeshData::~ImportMeshData()
{
    delete ui;
}

void ImportMeshData::LoadFolder()
{
    meshName.folderName=QFileDialog::getExistingDirectory(nullptr,"Choose a folder");
    meshName.coordFile=meshName.folderName+"/"+"coordinates.dat";
    meshName.eleFile=meshName.folderName+"/"+"elements.dat";
    meshName.fixXFile=meshName.folderName+"/"+"fixx.dat";
    meshName.fixYFile=meshName.folderName+"/"+"fixy.dat";
    meshName.fixZFile=meshName.folderName+"/"+"fixz.dat";
    meshName.fixHFile=meshName.folderName+"/"+"fixh.dat";
    meshName.forceYFile=meshName.folderName+"/"+"forcey.dat";
}

void ImportMeshData::LoadSubFile(int fileCase)
{
    QString fileName=QFileDialog::getOpenFileName(nullptr,"Choose a file");
    if(fileCase==1)
    {
        meshName.coordFile=fileName;
    }
    else if(fileCase==2)
    {
        meshName.eleFile=fileName;
    }
    else if(fileCase==3)
    {
        meshName.fixXFile=fileName;
    }
    else if(fileCase==4)
    {
        meshName.fixYFile=fileName;
    }
    else if(fileCase==5)
    {
        meshName.fixZFile=fileName;
    }
    else if(fileCase==6)
    {
        meshName.fixHFile=fileName;
    }
    else if (fileCase==7)
    {
        meshName.forceYFile=fileName;
    }
    else
    {
        qDebug()<<"Not supported yet"<<endl;
    }
    UpdateFileName();
}

void ImportMeshData::UpdateFileName()
{
    ui->coordLine->setText(meshName.coordFile);
    ui->eleLine->setText(meshName.eleFile);
    ui->fixXLine->setText(meshName.fixXFile);
    ui->fixYLine->setText(meshName.fixYFile);
    ui->fixZLine->setText(meshName.fixZFile);
    ui->fixHLine->setText(meshName.fixHFile);
    ui->forceYLine->setText(meshName.forceYFile);
    QFileInfo mFile(meshName.coordFile);
    meshName.folderName=mFile.path();
}

void ImportMeshData::on_pushButton_10_clicked()
{
    LoadFolder();
    UpdateFileName();
    qDebug()<<"folderName: "<<meshName.folderName<<endl;
}

void ImportMeshData::on_closeButton_clicked()
{
    this->close();
}

void ImportMeshData::on_okButton_clicked()
{
    emit SendMeshData(meshName.coordFile,meshName.eleFile);
    this->close();

}

void ImportMeshData::on_coordButton_clicked()
{
    fileCase=1;
    LoadSubFile(fileCase);
}

void ImportMeshData::on_eleButton_clicked()
{
    fileCase=2;
    LoadSubFile(fileCase);
}

void ImportMeshData::on_fixXButton_clicked()
{
    fileCase=3;
    LoadSubFile(fileCase);
}

void ImportMeshData::on_fixYButton_clicked()
{
    fileCase=4;
    LoadSubFile(fileCase);
}

void ImportMeshData::on_fixZButton_clicked()
{
    fileCase=5;
    LoadSubFile(fileCase);
}

void ImportMeshData::on_fixHButton_clicked()
{
    fileCase=6;
    LoadSubFile(fileCase);
}

void ImportMeshData::on_forceYButton_clicked()
{
    fileCase=7;
    LoadSubFile(fileCase);
}
