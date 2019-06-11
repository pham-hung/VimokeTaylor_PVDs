#ifndef PROGESS_H
#define PROGESS_H

#include <QDialog>
#include <QProgressBar>
#include <QDebug>

using namespace std;

namespace Ui {
class Progess;
}

class Progess : public QDialog
{
    Q_OBJECT

public:
    explicit Progess(QWidget *parent = 0);
    ~Progess();
    int step;
    int iterationNumber=1;
    int runningTime=0;
    int numberOfStep;
    double error=999999;
    void UpDateProgress();
    void UpdateFinalStep();

public slots:
    void GetCalculationInformation(int step, int numberOfStep);
private slots:
    void on_closeButton_clicked();

private:
    Ui::Progess *ui;
};

#endif // PROGESS_H
