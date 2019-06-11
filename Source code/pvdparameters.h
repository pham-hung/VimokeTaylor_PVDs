#ifndef PVDPARAMETERS_H
#define PVDPARAMETERS_H

#include <QDialog>
using namespace std;

struct PVD{
    double unitCellRadius=0.6;
    double PVDRadius=0.03;
    double smearZoneRadius=0.1;
    double dischargeCapacity=1e-4;
    double PVDLength=2.0;
    double PVDConductivity;
};

namespace Ui {
class PVDParameters;
}

class PVDParameters : public QDialog
{
    Q_OBJECT

public:
    explicit PVDParameters(QWidget *parent = 0);
    ~PVDParameters();
    void SendSignal(){emit SendPVDData(PVDData);}

signals:
    void SendPVDData(PVD PVDData);

private slots:
    void on_pushButton_2_clicked();
    void on_pushButton_clicked();

private:
    PVD PVDData;
    Ui::PVDParameters *ui;
    void GetUserData();
    void UpdateData();
};

#endif // PVDPARAMETERS_H
