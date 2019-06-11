#ifndef SOILPARAMETER_H
#define SOILPARAMETER_H

#include <QDialog>
#include <QString>
using namespace std;

struct Soil
{
    double re=1.0;
    double ksRatio=0.5;
    double kv=1e-9;
    double voidRatio=1.8;
    double K=500;
    double v=0.2;
    double Ch=5e-7;
};

namespace Ui {
class SoilParameter;
}

class SoilParameter : public QDialog
{
    Q_OBJECT

public:
    explicit SoilParameter(QWidget *parent = 0);
    ~SoilParameter();
    void SendSignal(){emit SendSoilData(soilData);}

signals:
    void SendSoilData(Soil soilData);

private slots:
    void on_okButton_clicked();

    void on_updateButton_clicked();

private:
    Ui::SoilParameter *ui;
    Soil soilData;
    void GetUserData();
    void UpdateData();
};

#endif // SOILPARAMETER_H
