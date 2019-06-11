#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QBoxLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QLabel>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QCheckBox>
#include <QFormLayout>
#include <QTextEdit>
#include <QFile>
#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QComboBox>
#include <QRgb>
#include <QMouseEvent>
#include <QFrame>
#include <QThread>
#include <QProcess>
#include <QRadioButton>
#include <QSpacerItem>
#include <Eigen/Dense>
#include <iostream>
#include <windows.h>
#include "glwidget.h"
#include "consolidation3d.h"
#include "importmeshdata.h"
#include "soilparameter.h"
#include "pvdparameters.h"
#include "analysistype.h"
#include "progess.h"
#include "importresults.h"
#include "getfile.h"

//Widgetclass
#include "scalesetting.h"
#include "contoursetting.h"
#include "cuttingplane.h"
#include "stepresult.h"
#include "base3dviewport.h"
#include "elementinformation.h"

using namespace std;
using namespace Eigen;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void SetUpLabels();
    void SetUpLayout();
    double GoldenSectionOpimization(double a,double b);
    double CalculateError(double Cd);
    void GetResultsFromTest();  //Get result from consolidation 3D class
    void CalculateElementIndex();

private slots:
    void on_actionMesh_triggered(); //Import mesh data
    void on_actionSoil_parameters_triggered(); //Import soil parameters
    void on_actionPVD_parameters_triggered();  //Import PVD parameters
    void on_actionAnalysis_Type_triggered();   //Import analysis type
    void on_actionFind_Correction_Factor_triggered(); //Golden section method
    void on_actionAnalyze_With_Defined_Correction_Factor_triggered(); //Find Cd
    void on_actionResult_files_triggered(); //Import Result Files
    void on_actionScale_Setting_triggered(); //Open Scale Setting
    void on_actionContour_Colour_Setting_triggered(); //Open Contour Setting
    void on_actionCutting_Plane_Setting_triggered();  //Open Cutting Plane Setting
    void on_actionTime_Step_Animation_triggered();    //Open Pick Time Setting

    //Control OpenGL
    void Replot();
    void ReplotResetViewport();


public slots:
    void GetResultFileNames(QString xDispFile, QString yDispFile, QString zDispFile, QString poreFile); //Import Results
    void GetMeshData(QString coordFile, QString eleFile); //Import Mesh

    //Slots to connect external Widgets
    void GetScaleInformation(BaseScaleClass scaleObject);  //Scale information
    void GetColorInformation(BaseContourSetting contourSettingObject);  //Colour contour information
    void GetCuttingPlaneInformation(BaseCuttingPlane cuttingPlaneSettingObject);  //Cutting plane information
    void NewPickedField(int fieldIndex);  //Choose another result field
    void NewPickedStep(int fieldIndex,int stepIndex); //Chose another result field and another step
    void GetStepResultInformation(BaseStepResult stepResultObject); //Information of each step results (max, min...)

    //Slots to connect OpenGL Widget
    void GetSignalFromOpenGLWidget(); //Start OpenGL Widget
    void GetValueAtMousePosition(double mouseValue); //Getting value when click mouse
    void GetViewportInformation(Base3DViewport viewPortSettingObject); //Getting viewport information
    void ShowResultWindowPick(){stepResultSetting->show();}

signals:
    void SendBackInformationNewPickedStep(int,bool);      //send step
    void SendBackInformationNewPickedStep(double,double); //send max min values


    //signals to OpenGL Widget
    void SendShownData(Ref<MatrixXd> coordinatesScale, Ref<MatrixXd> coordinatesDeform, Ref<MatrixXd> elements, Ref<MatrixXd> shownData,ElementInfor elementInfor);
    void SendViewportInformation(Base3DViewport viewPortSettingObject);
    void SendContourSettingInformation(BaseContourSetting contourSettingObject);
    void SendCuttingPlaneInformation(BaseCuttingPlane cuttingPlaneSettingObject);
    void SendStepResultInformation(BaseStepResult stepResultObject);

private:
    Ui::MainWindow *ui;
    GLWidget *widget;

    //Layout
    QGridLayout *mainLayout=new QGridLayout;
    QVBoxLayout *glLayout = new QVBoxLayout;
    QBoxLayout *bottomLayout= new QBoxLayout(QBoxLayout::LeftToRight);

    //Bottom buttons
    QPushButton *replotButton= new QPushButton;       
    QPushButton *topButton= new QPushButton;
    QPushButton *botButton= new QPushButton;
    QPushButton *rightButton= new QPushButton;
    QPushButton *leftButton= new QPushButton;
    QPushButton *behindButton= new QPushButton;
    QPushButton *frontButton= new QPushButton;
    QComboBox *resultBox =new QComboBox;

    //Spacer
    QSpacerItem *spacer1 = new QSpacerItem(20,10,QSizePolicy::Expanding,QSizePolicy::Minimum);
    Consolidation3D *test=new Consolidation3D;

    //Import Data Classes
    ImportMeshData *importMesh=new ImportMeshData;
    SoilParameter *importSoil=new SoilParameter;
    PVDParameters *importPVD = new PVDParameters;
    AnalysisType *importType=new AnalysisType;
    Progess *calculationStatus=new Progess;
    ImportResults *importResult=new ImportResults;

    //Golden search results
    double Cd;
    double minError;

    //Timer
    QElapsedTimer timer;

    //OpenGL variables
    bool resultValid=false;

    //Results from consolidation analysis
    MatrixXd xDisp, yDisp, zDisp, Pore;
    QString xDispFile, yDispFile, zDispFile, poreFile;
    MatrixXd shownData=MatrixXd::Zero(1,1); //XDisp, or YDisp, or ZDisp, or what ever

    //Coordinates and elements
    MatrixXd coordinates=MatrixXd::Zero(1,4);
    MatrixXd coordinatesScale=MatrixXd::Zero(1,4);
    MatrixXd coordinatesDeform=MatrixXd::Zero(1,4); //after deformon
    MatrixXd elements=MatrixXd::Zero(1,1);

    //View control objects
    BaseScaleClass scaleObject;     //setting scale view
    BaseContourSetting contourSettingObject; //colorBand setting
    BaseCuttingPlane cuttingPlaneSettingObject; //Cutting Plane
    BaseStepResult stepResultObject; //Step Result
    Base3DViewport viewPortSettingObject;

    //View control widgets
    ScaleSetting *scaleView=new ScaleSetting;
    ContourSetting *colorSetting=new ContourSetting;
    CuttingPlane *cuttingPlaneSetting=new CuttingPlane;
    StepResult *stepResultSetting = new StepResult;

    //Mesh informatin
    int numberOfNodes, numberOfElements, numberOfStep;
    ElementInfor elementInfor;
    int fieldIndex;

    //private functions
    void SetUpSignalSlotConnections();
    void UpdateScaleCoordinates();
    void UpdateDeformationCoordinates();
    void UpdateShownInOpenGLData();
    void RunAnimation();
};

#endif // MAINWINDOW_H

//Variables start with normal letter
//Function start with Capital letter
