#ifndef CUTTINGPLANE_H
#define CUTTINGPLANE_H

#include <QDialog>
#include <QVector4D>
#include <QVector3D>

struct BaseCuttingPlane
{
    //Setting a cutting plane
    bool clipFlag=false;
    double normX=0;
    double normY=0;
    double normZ=1;
    double pointX=0;
    double pointY=0;
    double pointZ=0;
    QVector4D normVector;
    QVector3D pointVector;
    void setToVector();
};

namespace Ui {
class CuttingPlane;
}

class CuttingPlane : public QDialog
{
    Q_OBJECT

public:
    explicit CuttingPlane(QWidget *parent = 0);
    ~CuttingPlane();

signals:
    void SendClipInformation(BaseCuttingPlane cuttingPlaneSettingObject);

private slots:
    void on_closeButton_clicked();
    void on_okButton_clicked();

private:
    Ui::CuttingPlane *ui;
    BaseCuttingPlane cuttingPlaneSettingObject;
};

#endif // CUTTINGPLANE_H
