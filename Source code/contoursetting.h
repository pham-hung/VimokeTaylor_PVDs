#ifndef CONTOURSETTING_H
#define CONTOURSETTING_H

#include <QDialog>
#include <QString>
using namespace std;

struct BaseContourSetting
{
    int noc=10; //number of contour
    int contourPosition=0; //contourPosition, 0 is left, 1 is right, 2 is bottom
    int fontSize=10;
    int numericType=0; //1 is float (3 decimal), 1 is interger, 2 is scientific
    bool nodeVisible=false;
    bool meshVisible=false;
    bool resultVisible=true;
    bool axeVisible=true;
    bool dateVisible=true;
    bool titleVisible=true;
    bool maxVisible=true;
    bool deformedMeshVisible=true;
    QString title="Unit is m, Kpa, s";
};

namespace Ui {
class ContourSetting;
}

class ContourSetting : public QDialog
{
    Q_OBJECT

public:
    explicit ContourSetting(QWidget *parent = 0);
    ~ContourSetting();

private slots:
    void on_okButton_clicked();
    void on_closeButton_clicked();

signals:
    void SendColorInformation(BaseContourSetting contourSettingObject);

private:
    Ui::ContourSetting *ui;
    BaseContourSetting contourSettingObject;
};

#endif // CONTOURSETTING_H
