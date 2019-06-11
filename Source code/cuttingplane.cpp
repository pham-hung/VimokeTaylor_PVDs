#include "cuttingplane.h"
#include "ui_cuttingplane.h"

CuttingPlane::CuttingPlane(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::CuttingPlane)
{
    ui->setupUi(this);
}

CuttingPlane::~CuttingPlane()
{
    delete ui;
}

void CuttingPlane::on_closeButton_clicked()
{
    this->close();
}

void CuttingPlane::on_okButton_clicked()
{
    cuttingPlaneSettingObject.normX=ui->xNormLine->text().toDouble();
    cuttingPlaneSettingObject.normY=ui->yNormLine->text().toDouble();
    cuttingPlaneSettingObject.normZ=ui->zNormLine->text().toDouble();
    cuttingPlaneSettingObject.pointX=ui->xCoordLine->text().toDouble();
    cuttingPlaneSettingObject.pointY=ui->yCoordLine->text().toDouble();
    cuttingPlaneSettingObject.pointZ=ui->zCoordLine->text().toDouble();

    if(ui->clipFlagCheckBox->isChecked())
    {
        cuttingPlaneSettingObject.clipFlag=true;
    }
    else
    {
        cuttingPlaneSettingObject.clipFlag=false;
    }

    cuttingPlaneSettingObject.setToVector();
    emit SendClipInformation(cuttingPlaneSettingObject);
    this->close();
}

void BaseCuttingPlane::setToVector()
{
    normVector.setX(normX);
    normVector.setY(normY);
    normVector.setZ(normZ);

    if(clipFlag==true)
    {
        normVector.setW(1);
    }
    else
    {
        normVector.setW(-1);
    }

    pointVector.setX(pointX);
    pointVector.setY(pointY);
    pointVector.setZ(pointZ);
}
