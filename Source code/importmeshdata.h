#ifndef IMPORTMESHDATA_H
#define IMPORTMESHDATA_H

#include <QDialog>
#include <QString>
#include <QDebug>
#include <QFileDialog>
#include <QFileInfo>

using namespace std;

struct meshDataNames{
    QString coordFile;
    QString eleFile;
    QString fixXFile;
    QString fixYFile;
    QString fixZFile;
    QString fixHFile;
    QString forceYFile;
    QString folderName;
};

namespace Ui {
class ImportMeshData;
}

class ImportMeshData : public QDialog
{
    Q_OBJECT

public:
    explicit ImportMeshData(QWidget *parent = 0);
    ~ImportMeshData();
    void SendSignal(){emit SendMeshFiles(meshName);}

signals:
    void SendMeshFiles(meshDataNames meshName);
    void SendMeshData(QString, QString);

public slots:

private slots:
    void on_pushButton_10_clicked();
    void on_closeButton_clicked();
    void on_okButton_clicked();
    void on_coordButton_clicked();
    void on_eleButton_clicked();
    void on_fixXButton_clicked();
    void on_fixYButton_clicked();
    void on_fixZButton_clicked();
    void on_fixHButton_clicked();
    void on_forceYButton_clicked();

private:
    Ui::ImportMeshData *ui;
    meshDataNames meshName;
    int fileCase;

    void LoadFolder();
    void LoadSubFile(int fileCase);
    void UpdateFileName();
};

#endif // IMPORTMESHDATA_H
