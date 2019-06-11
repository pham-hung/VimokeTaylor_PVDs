#include "progess.h"
#include "ui_progess.h"

Progess::Progess(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Progess)
{
    ui->setupUi(this);
}

Progess::~Progess()
{
    delete ui;
}

void Progess::UpDateProgress()
{
    ui->errorLabel->setText("Difference between two iterations");
    ui->progressBar->setValue(int(100*step/numberOfStep));
    ui->iterationLine->setText(QString::number(iterationNumber));
    ui->errorLine->setText(QString::number(error));
    ui->runningTimeLine->setText(QString::number(runningTime));
    QCoreApplication::processEvents();
}

void Progess::UpdateFinalStep()
{
    ui->errorLine->setText(QString::number(error));
    ui->errorLabel->setText("Minimum error");
    ui->runningTimeLine->setText(QString::number(runningTime));
    QCoreApplication::processEvents();
}


void Progess::on_closeButton_clicked()
{
    this->close();
}

void Progess::GetCalculationInformation(int step, int numberOfStep)
{
    this->step=step;
    this->numberOfStep=numberOfStep;
    UpDateProgress();
}
