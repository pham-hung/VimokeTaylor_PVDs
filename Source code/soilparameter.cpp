#include "soilparameter.h"
#include "ui_soilparameter.h"

SoilParameter::SoilParameter(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SoilParameter)
{
    ui->setupUi(this);
}

SoilParameter::~SoilParameter()
{
    delete ui;
}

void SoilParameter::GetUserData()
{
    soilData.kv=ui->kvLine->text().toDouble();
    soilData.re=ui->ratioKhLine->text().toDouble();
    soilData.ksRatio=ui->ratioKsLine->text().toDouble();
    soilData.voidRatio=ui->voidRatioLine->text().toDouble();
    soilData.v=ui->PoissonRatioLine->text().toDouble();
    soilData.K=ui->KLine->text().toDouble();
    soilData.Ch=ui->ChLine->text().toDouble();
}

void SoilParameter::UpdateData()
{
    GetUserData();
    if(soilData.Ch>0)
    {
        double gf=9.81;
        double mv=soilData.re*soilData.kv/soilData.Ch/gf;
        soilData.K=(1/mv)/(2*(1-2*soilData.v)/(1+soilData.v)+1);
        ui->KLine->setText(QString::number(soilData.K));
    }
    else
    {
        double gf=9.81;
        double G;
        G=3.0f*soilData.K*(1-2.0f*soilData.v)/2.0f/(1+soilData.v);
        double mv;
        mv=1.0f/(soilData.K+4.0f*G/3.0f);
        soilData.Ch=soilData.re*soilData.kv/mv/gf;
        ui->ChLine->setText(QString::number(soilData.Ch));
    }
    GetUserData();
}

void SoilParameter::on_okButton_clicked()
{
    UpdateData();
    this->close();
}

void SoilParameter::on_updateButton_clicked()
{
    UpdateData();
}
