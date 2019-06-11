#ifndef SCALESETTING_H
#define SCALESETTING_H

#include <QDialog>
#include <iostream>

namespace Ui {
class ScaleSetting;
}

struct BaseScaleClass{
    double xScale=1;
    double yScale=1;
    double zScale=1;
    double deformScale=1;
};

class ScaleSetting : public QDialog
{
    Q_OBJECT

public:
    explicit ScaleSetting(QWidget *parent = 0);
    ~ScaleSetting();

signals:
    void SendScaleInformation(BaseScaleClass scaleObject);
    void ChooseNewField(int fieldPosition);
    void ChooseNewStep(int step);

private slots:
    void on_closeButton_clicked();
    void on_applyButton_clicked();

private:
    Ui::ScaleSetting *ui;
    BaseScaleClass scaleObject;
};

#endif // SCALESETTING_H
