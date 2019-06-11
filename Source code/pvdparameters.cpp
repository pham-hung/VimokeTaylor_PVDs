#include "pvdparameters.h"
#include "ui_pvdparameters.h"

PVDParameters::PVDParameters(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PVDParameters)
{
    ui->setupUi(this);
}

PVDParameters::~PVDParameters()
{
    delete ui;
}

void PVDParameters::UpdateData()
{
    GetUserData();
    double PI=3.14159265359;
    double PVDarea;
    PVDarea=PVDData.PVDRadius*PVDData.PVDRadius*PI;
    PVDData.PVDConductivity=PVDData.dischargeCapacity/PVDarea;
    ui->kwLine->setText(QString::number(PVDData.PVDConductivity));
    ui->AwLine->setText(QString::number(PVDarea));
}

void PVDParameters::on_pushButton_2_clicked()
{
    UpdateData();
}

void PVDParameters::on_pushButton_clicked()
{
    UpdateData();
    this->close();
}

void PVDParameters::GetUserData()
{
    PVDData.unitCellRadius=ui->ReLine->text().toDouble();
    PVDData.PVDRadius=ui->RwLine->text().toDouble();
    PVDData.smearZoneRadius=ui->RsLine->text().toDouble();
    PVDData.PVDLength=ui->HLine->text().toDouble();
    PVDData.dischargeCapacity=ui->qwLine->text().toDouble();
}
