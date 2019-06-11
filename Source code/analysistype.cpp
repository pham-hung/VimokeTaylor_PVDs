#include "analysistype.h"
#include "ui_analysistype.h"

AnalysisType::AnalysisType(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AnalysisType)
{
    ui->setupUi(this);
}

AnalysisType::~AnalysisType()
{
    delete ui;
}

void AnalysisType::on_closeButton_clicked()
{
    this->close();
}

void AnalysisType::on_okButton_clicked()
{
    GetUserData();
    this->close();
}

void AnalysisType::GetUserData()
{
    typeData.analysisType=ui->typeCombo->currentIndex();
    typeData.startValue=ui->startLine->text().toDouble();
    typeData.endValue=ui->endLine->text().toDouble();
    typeData.Cd0=ui->Cd0Line->text().toDouble();
    typeData.tolerance=ui->tolLine->text().toDouble();
    typeData.numberOfStep=ui->numberOfStepLine->text().toDouble();
    typeData.timeIncrement=ui->dtLine->text().toDouble();
    typeData.intialPressure=ui->p0Line->text().toDouble();
}
