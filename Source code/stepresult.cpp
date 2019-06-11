#include "stepresult.h"
#include "ui_stepresult.h"

StepResult::StepResult(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::StepResult)
{
    ui->setupUi(this);
    ui->resultCombo->setCurrentIndex(-1);
}

StepResult::~StepResult()
{
    delete ui;
}

void StepResult::DoneAnimation()
{
    stepResultObject.animationFlag=false;
}

void StepResult::on_animationStart_clicked()
{
    GetUserData();
    stepResultObject.animationFlag=true;
    emit SendStepResultsInformation(stepResultObject);
    this->close();
}


void StepResult::on_resultCombo_currentIndexChanged(int index)
{
    emit PickedNewField(index);
}

void StepResult::GetInformationFromMainWindow(int nos, bool validPick)
{    
    this->nos=nos;
    this->validPick=validPick;
    if(validPick==false)
    {
        QMessageBox::warning(Q_NULLPTR,"ERROR","This selection has no data, import data first");
        ui->stepCombo->clear();
        ui->startCombo->clear();
        ui->endCombo->clear();
    }
    else
    {
        ui->stepCombo->clear();
        for (int i=0;i<nos;i++)
        {
            QString newStep=QString::number(i+1,10);
            ui->stepCombo->addItem(newStep);
            ui->startCombo->addItem(newStep);
            ui->endCombo->addItem(newStep);
        }
    }
}

void StepResult::GetInformationFromMainWindow(double minVal, double maxVal)
{
    maxValAuto=maxVal;
    minValAuto=minVal;
    if(ui->maxCombo->currentIndex()==0) //Auto
    {
        ui->minLine->setText(QString::number(minValAuto,'g',6));
        ui->maxLine->setText(QString::number(maxValAuto,'g',6));
        ui->maxLine->setEnabled(false);
        ui->minLine->setEnabled(false);
    }
    on_okButton_clicked();
}

void StepResult::on_closeButton_clicked()
{
    this->close();
}

void StepResult::on_stepCombo_currentIndexChanged(int index)
{    
    int fieldIndex=ui->resultCombo->currentIndex();
    int stepIndex=ui->stepCombo->currentText().toInt();
    if(index>=0)
    {
        emit PickedNewStep(fieldIndex,stepIndex);
    }
}

void StepResult::on_maxCombo_currentIndexChanged(int index)
{
    if(ui->maxCombo->currentIndex()==1) //manualy
    {
        ui->maxLine->setEnabled(true);
        ui->maxLine->setText("100");
        ui->minLine->setEnabled(true);
        ui->minLine->setText("-100");
    }
    else
    {
        ui->maxLine->setEnabled(false);
        ui->minLine->setEnabled(false);
        ui->minLine->setText(QString::number(minValAuto,'g',6));
        ui->maxLine->setText(QString::number(maxValAuto,'g',6));
    }
}

void StepResult::on_okButton_clicked()
{    
    GetUserData();
    stepResultObject.animationFlag=false;
    emit SendStepResultsInformation(stepResultObject);
}

void StepResult::GetUserData()
{
    stepResultObject.fieldIndex=ui->resultCombo->currentIndex();
    stepResultObject.stepIndex=ui->stepCombo->currentText().toInt();
    if(ui->maxCombo->currentIndex()==0)
    {
        stepResultObject.autoMaxMin=true;
        ui->maxLine->setEnabled(false);
        ui->minLine->setEnabled(false);
        ui->minLine->setText(QString::number(minValAuto,'g',6));
        ui->maxLine->setText(QString::number(maxValAuto,'g',6));
    }
    else
    {
        stepResultObject.autoMaxMin=false;
    }
    stepResultObject.maxValue=ui->maxLine->text().toDouble();
    stepResultObject.minValue=ui->minLine->text().toDouble();
    stepResultObject.startStep=ui->startCombo->currentText().toInt();
    stepResultObject.endStep=ui->endCombo->currentText().toInt();
    stepResultObject.delay=ui->delayLine->text().toInt();

    if(ui->autoChangeValueCheckBox->isChecked())
    {
        stepResultObject.autoChangeMaxMin=true;
    }
    else
    {
        stepResultObject.autoChangeMaxMin=false;
    }
}
