#include "contoursetting.h"
#include "ui_contoursetting.h"


ContourSetting::ContourSetting(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ContourSetting)
{
    ui->setupUi(this);
}

ContourSetting::~ContourSetting()
{
    delete ui;
}

void ContourSetting::on_okButton_clicked()
{
    contourSettingObject.noc=ui->nocLine->text().toInt();
    if(contourSettingObject.noc<1)
    {
        contourSettingObject.noc=10;
        ui->nocLine->setText("10");
    }

    contourSettingObject.contourPosition=ui->positionCombo->currentIndex();
    contourSettingObject.fontSize=ui->fontSizeLine->text().toInt();
    if(contourSettingObject.fontSize<0)
    {
        contourSettingObject.fontSize=10;
        ui->fontSizeLine->setText("10");
    }
    contourSettingObject.numericType=ui->numericCombo->currentIndex();

    if(ui->nodeCheckBox->isChecked())
    {
        contourSettingObject.nodeVisible=true;
    }
    else
    {
        contourSettingObject.nodeVisible=false;
    }

    if(ui->meshCheckBox->isChecked())
    {
        contourSettingObject.meshVisible=true;
    }
    else
    {
        contourSettingObject.meshVisible=false;
    }

    if(ui->resultCheckBox->isChecked())
    {
        contourSettingObject.resultVisible=true;
    }
    else
    {
        contourSettingObject.resultVisible=false;
    }

    if(ui->axeCheckBox->isChecked())
    {
        contourSettingObject.axeVisible=true;
    }
    else
    {
        contourSettingObject.axeVisible=false;
    }

    if(ui->dateCheckBox->isChecked())
    {
        contourSettingObject.dateVisible=true;
    }
    else
    {
        contourSettingObject.dateVisible=false;
    }


    if(ui->titleCheckBox->isChecked())
    {
        contourSettingObject.titleVisible=true;
    }
    else
    {
        contourSettingObject.titleVisible=false;
    }

    contourSettingObject.title=ui->titleLine->text();

    if(ui->maxCheckBox->isChecked())
    {
        contourSettingObject.maxVisible=true;
    }
    else
    {
        contourSettingObject.maxVisible=false;
    }
    if(ui->deformMeshCheckBox->isChecked())
    {
        contourSettingObject.deformedMeshVisible=true;
    }
    else
    {
        contourSettingObject.deformedMeshVisible=false;
    }
    emit SendColorInformation(contourSettingObject);
    this->close();
}

void ContourSetting::on_closeButton_clicked()
{
    this->close();
}
