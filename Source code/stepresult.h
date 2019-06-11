#ifndef STEPRESULT_H
#define STEPRESULT_H

#include <QDialog>
#include <QMessageBox>
#include <QDebug>
#include <Eigen/Dense>
#include <QString>

namespace Ui {
class StepResult;
}
struct BaseStepResult{
    int fieldIndex=0;
    int stepIndex=1;
    QString resultComboName="EPWP-Excess Pore Pressure";
    double minValue=0;
    double maxValue=100;
    bool autoMaxMin=true;
    int startStep=1;
    int endStep=1;
    double delay=200;
    bool lockView=false;
    bool autoChangeMaxMin=true;
    bool animationFlag=false;
    bool validInformation=false;
};

class StepResult : public QDialog
{
    Q_OBJECT

public:
    explicit StepResult(QWidget *parent = 0);
    ~StepResult();
    void DoneAnimation();

public slots:
    void GetInformationFromMainWindow(int nos,bool validPick);
    void GetInformationFromMainWindow(double minVal, double maxVal);

private slots:
    void on_animationStart_clicked();
    void on_resultCombo_currentIndexChanged(int index); //change combo pick
    void on_closeButton_clicked();
    void on_stepCombo_currentIndexChanged(int index); //change step combo
    void on_maxCombo_currentIndexChanged(int index);
    void on_okButton_clicked();

signals:
    void SendStepResultsInformation(BaseStepResult stepResultObject);
    void PickedNewField(int index);
    void PickedNewStep(int fieldIndex,int stepIndex);

private:
    void GetUserData();
    Ui::StepResult *ui;
    BaseStepResult stepResultObject;
    int nos;
    bool validPick=true;
    double minValAuto, maxValAuto;
};

#endif // STEPRESULT_H
