#include "importresults.h"
#include "ui_importresults.h"

ImportResults::ImportResults(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ImportResults)
{
    ui->setupUi(this);
}

ImportResults::~ImportResults()
{
    delete ui;
}

void ImportResults::on_closeButton_clicked()
{
    this->close();
}

void ImportResults::on_startImportButton_clicked()
{
    emit SendResultFiles(xDispFile,yDispFile,zDispFile,poreFile);
}

void ImportResults::on_browseXButton_clicked()
{
    xDispFile=QFileDialog::getOpenFileName();
    ui->xFileLine->setText(xDispFile);
}

void ImportResults::on_browseYButton_clicked()
{
    yDispFile=QFileDialog::getOpenFileName();
    ui->yFileLine->setText(yDispFile);
}

void ImportResults::on_browseZButton_clicked()
{
    zDispFile=QFileDialog::getOpenFileName();
    ui->zFileLine->setText(zDispFile);
}

void ImportResults::on_browsePButton_clicked()
{
    poreFile=QFileDialog::getOpenFileName();
    ui->poreFileLine->setText(poreFile);
}
