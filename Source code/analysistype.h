#ifndef ANALYSISTYPE_H
#define ANALYSISTYPE_H

#include <QDialog>
#include <QDebug>
using namespace std;

struct SolutionType{
    int analysisType;
    double startValue;
    double endValue;
    double Cd0;
    double tolerance;
    double numberOfStep;
    double timeIncrement;
    double intialPressure;
};

namespace Ui {
class AnalysisType;
}

class AnalysisType : public QDialog
{
    Q_OBJECT

public:
    explicit AnalysisType(QWidget *parent = 0);
    ~AnalysisType();
    void SendSignal(){emit SendTypeData(typeData);}
    SolutionType typeData;

private slots:
    void on_closeButton_clicked();
    void on_okButton_clicked();

signals:
    void SendTypeData(SolutionType typeData);
private:
    Ui::AnalysisType *ui;    
    void GetUserData();
};

#endif // ANALYSISTYPE_H
