#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    SetUpLabels();
    SetUpLayout();    

    connect(importMesh,SIGNAL(SendMeshFiles(meshDataNames)),test,SLOT(GetMeshData(meshDataNames)));
    connect(importSoil,SIGNAL(SendSoilData(Soil)),test,SLOT(GetSoilData(Soil)));
    connect(importPVD,SIGNAL(SendPVDData(PVD)),test,SLOT(GetPVDData(PVD)));
    connect(importType,SIGNAL(SendTypeData(SolutionType)),test,SLOT(GetAnalysisType(SolutionType)));
    connect(test,SIGNAL(sendCalculationInfor(int,int)),calculationStatus,SLOT(GetCalculationInformation(int,int)));
    connect(importResult,SIGNAL(SendResultFiles(QString,QString,QString,QString)),this,SLOT(GetResultFileNames(QString,QString,QString,QString)));
    connect(importMesh,SIGNAL(SendMeshData(QString,QString)),this,SLOT(GetMeshData(QString,QString)));
    connect(resultBox,SIGNAL(currentIndexChanged(int)),this,SLOT(ShowResultWindowPick()));

    SetUpSignalSlotConnections();
    Replot();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::SetUpLabels()
{
    this->setWindowTitle("Finding correction factor of Vimoke-Taylor concept for PVD analysis");
    resultBox->addItem("X-Displacement");
    resultBox->addItem("Y-Displacement");
    resultBox->addItem("Z-Displacement");
    resultBox->addItem("EPWP-Excess Pore Pressure");
    resultBox->addItem("Material number");
    resultBox->addItem("Element Type");    

    replotButton->setText("Replot");
    topButton->setText("TOP");
    botButton->setText("BOT");
    rightButton->setText("RIGHT");
    leftButton->setText("LEFT");
    behindButton->setText("BEHIND");
    frontButton->setText("FRONT");

}

void MainWindow::SetUpLayout()
{
    this->setGeometry(200,200,1366,768);
    ui->centralWidget->setLayout(mainLayout);

    //Layout
    mainLayout->addLayout(glLayout,0,0,18,10);
    mainLayout->addLayout(bottomLayout,19,0,1,10);
    bottomLayout->addWidget(resultBox);
    bottomLayout->addWidget(replotButton);
    bottomLayout->addWidget(topButton);
    bottomLayout->addWidget(botButton);
    bottomLayout->addWidget(rightButton);
    bottomLayout->addWidget(leftButton);
    bottomLayout->addWidget(behindButton);
    bottomLayout->addWidget(frontButton);
    bottomLayout->addItem(spacer1);

    //Width
    resultBox->setFixedWidth(200);
    replotButton->setFixedWidth(70);
    topButton->setFixedWidth(70);
    botButton->setFixedWidth(70);
    rightButton->setFixedWidth(70);
    leftButton->setFixedWidth(70);
    behindButton->setFixedWidth(70);
    frontButton->setFixedWidth(70);

}

double MainWindow::GoldenSectionOpimization(double a, double b)
{
    //a is start value
    //b is end value
    //return Cdmin
    double tol=importType->typeData.tolerance;
    double CdMinError;

    int i=0;
    double gr=(1+sqrt(5.0))/2.0;
    double c=b-(b-a)/gr;
    double d=a+(b-a)/gr;
    double fc=CalculateError(c);
    double fd=CalculateError(d);
    i=i+2;

    while(abs(c-d)>tol)
    {
        if(fc<fd)
        {
            b=d;
            d=c;
            fd=fc;
            c=b-(b-a)/gr;
            fc=CalculateError(c);
        }
        else
        {
            a=c;
            c=d;
            fc=fd;
            d=a+(b-a)/gr;
            fd=CalculateError(d);
        }
        calculationStatus->runningTime=timer.elapsed()/1000;
        calculationStatus->error=abs(c-d);
        calculationStatus->iterationNumber=i;
        i=i+1;
    }
    CdMinError=(a+b)/2.0f;
    return CdMinError;
}

double MainWindow::CalculateError(double Cd)
{
    test->setCd(Cd);
    test->drainedSolve();
    test->calculateError();
    return test->getError();
}

void MainWindow::GetResultsFromTest()
{
    xDisp=test->GetXDisp();
    yDisp=test->GetYDisp();
    zDisp=test->GetZDisp();
    Pore=test->GetPore();
    test->clearData();
}

void MainWindow::CalculateElementIndex()
{
    numberOfElements=elements.col(0).maxCoeff();
    numberOfNodes=coordinates.rows();
    shownData.resize(numberOfNodes,1);
    shownData.setOnes();

    elementInfor.resetCountNumber();
    elementInfor.elmementIndex=MatrixXi::Zero(numberOfElements,1);
    int elementCount=0;

    for (int j=0;j<elementInfor.elmementIndex.rows();j++)
    {
        int eleIndex=elements(j,0);
        if(eleIndex!=0)
        {
            elementInfor.elmementIndex(elementCount,0)=j;
            elementCount=elementCount+1;
        }
    }

    for (int j=0; j<elementInfor.elmementIndex.rows();j++)
    {
        int eleIndex=elementInfor.elmementIndex(j,0);
        int nodeCount=elements(eleIndex,11);
        if(nodeCount==20){elementInfor.hex20Count=elementInfor.hex20Count+1;}
        else if(nodeCount==8){elementInfor.hex8Count=elementInfor.hex8Count+1;}
        else if(nodeCount==10){elementInfor.tet10Count=elementInfor.tet10Count+1;}
        else if(nodeCount==5){elementInfor.pyra5Count=elementInfor.pyra5Count+1;}
        else if(nodeCount==13){elementInfor.pyra13Count=elementInfor.pyra13Count+1;}
        else if(nodeCount==15){elementInfor.prism15Count=elementInfor.prism15Count+1;}
        else if(nodeCount==6){elementInfor.prism6Count=elementInfor.prism6Count+1;}
    }

}

//Import Mesh and boundary conditions data
void MainWindow::on_actionMesh_triggered()
{
    importMesh->show();
}

//Import Soil parameters
void MainWindow::on_actionSoil_parameters_triggered()
{
    importSoil->show();
}

//Import PVD parameters
void MainWindow::on_actionPVD_parameters_triggered()
{
    importPVD->show();
}

//Import Analysis Type
void MainWindow::on_actionAnalysis_Type_triggered()
{
    importType->show();
}

void MainWindow::on_actionFind_Correction_Factor_triggered()
{
    calculationStatus->show();
    timer.start();
    Cd=importType->typeData.Cd0;

    test->resetData();
    importMesh->SendSignal();
    importSoil->SendSignal();
    importPVD->SendSignal();
    importType->SendSignal();
    double startValue, endValue;
    startValue=importType->typeData.startValue;
    endValue=importType->typeData.endValue;

    test->importData();
    test->initialData();
    test->createBoundaryCondition();
    test->calculateAnalytical();
    test->setCd(1.0f);
    test->undrainedSolve();

    Cd=GoldenSectionOpimization(startValue,endValue);
    test->resetData();
    test->setCd(Cd);
    test->undrainedSolve();
    test->drainedSolve();
    test->calculateError();
    test->exportResults();
    calculationStatus->error=test->getError();
    calculationStatus->UpdateFinalStep();
    GetResultsFromTest();
}

void MainWindow::on_actionAnalyze_With_Defined_Correction_Factor_triggered()
{
    calculationStatus->show();
    calculationStatus->iterationNumber=1;
    timer.start();
    Cd=importType->typeData.Cd0;
    test->resetData();
    importMesh->SendSignal();
    importSoil->SendSignal();
    importPVD->SendSignal();
    importType->SendSignal();
    test->importData();
    test->initialData();
    test->createBoundaryCondition();
    test->calculateAnalytical();
    test->setCd(Cd);
    test->undrainedSolve();
    test->drainedSolve();
    test->calculateError();
    test->exportResults();
    calculationStatus->runningTime=timer.elapsed()/1000;
    calculationStatus->error=test->getError();
    calculationStatus->UpdateFinalStep();
    GetResultsFromTest();
}

void MainWindow::on_actionResult_files_triggered()
{
    importResult->show();
}

void MainWindow::GetResultFileNames(QString xDispFile, QString yDispFile, QString zDispFile, QString poreFile)
{
    this->xDispFile=xDispFile;
    this->yDispFile=yDispFile;
    this->zDispFile=zDispFile;
    this->poreFile=poreFile;

    GetFile importFile;
    importFile.DoGetFile(xDispFile);
    xDisp=MatrixXd::Zero(importFile.row,importFile.col);
    xDisp=importFile.data_file;

    importFile.DoGetFile(yDispFile);
    yDisp=MatrixXd::Zero(importFile.row,importFile.col);
    yDisp=importFile.data_file;

    importFile.DoGetFile(zDispFile);
    zDisp=MatrixXd::Zero(importFile.row,importFile.col);
    zDisp=importFile.data_file;

    importFile.DoGetFile(poreFile);
    Pore=MatrixXd::Zero(importFile.row,importFile.col);
    Pore=importFile.data_file;

    UpdateDeformationCoordinates();
    ReplotResetViewport();
}

void MainWindow::GetMeshData(QString coordFile, QString eleFile)
{
    GetFile importFile;
    importFile.DoGetFile(coordFile);
    coordinates=MatrixXd::Zero(importFile.row,importFile.col);
    coordinates=importFile.data_file;

    importFile.DoGetFile(eleFile);
    elements=MatrixXd::Zero(importFile.row,importFile.col);
    elements=importFile.data_file;

    CalculateElementIndex();
    UpdateScaleCoordinates();
    UpdateDeformationCoordinates();
    ReplotResetViewport();
}

void MainWindow::GetScaleInformation(BaseScaleClass scaleObject)
{
    this->scaleObject=scaleObject;
    UpdateScaleCoordinates();
    UpdateDeformationCoordinates();
    ReplotResetViewport();
}

void MainWindow::GetColorInformation(BaseContourSetting contourSettingObject)
{
    this->contourSettingObject=contourSettingObject;
    viewPortSettingObject.lockView=true;
    Replot();
}

void MainWindow::GetCuttingPlaneInformation(BaseCuttingPlane cuttingPlaneSettingObject)
{
    this->cuttingPlaneSettingObject=cuttingPlaneSettingObject;
    viewPortSettingObject.lockView=true;
    Replot();
}

void MainWindow::NewPickedField(int fieldIndex)
{
    bool validPick=true;
    if(fieldIndex==0) //X-disp
    {
        numberOfStep=xDisp.cols();
        if(xDisp.rows()==1)
        {
            validPick=false;
        }
    }
    else if (fieldIndex==1) //Y-Disp
    {
        numberOfStep=yDisp.cols();
        if(yDisp.rows()==1)
        {
            validPick=false;
        }
    }
    else if (fieldIndex==2) //Z-Displacement
    {
        numberOfStep=zDisp.cols();
        if(zDisp.rows()==1)
        {
            validPick=false;
        }
    }
    else if (fieldIndex==3) //EPWP
    {
        numberOfStep=Pore.cols();
        if(Pore.rows()==1)
        {
            validPick=false;
        }
    }
    else if(fieldIndex==4 || fieldIndex==5)
    {
        numberOfStep=1;
        validPick=true;
    }

    if(validPick==true)
    {
        this->resultBox->setCurrentIndex(fieldIndex);
    }
    emit SendBackInformationNewPickedStep(numberOfStep,validPick);
}

void MainWindow::NewPickedStep(int fieldIndex, int stepIndex)
{
    double minVal, maxVal;
    if(fieldIndex==0) //Xdisp
    {
        minVal=xDisp.col(stepIndex-1).minCoeff();
        maxVal=xDisp.col(stepIndex-1).maxCoeff();
    }
    else if (fieldIndex==1) //Ydisp
    {
        minVal=yDisp.col(stepIndex-1).minCoeff();
        maxVal=yDisp.col(stepIndex-1).maxCoeff();
    }
    else if (fieldIndex==2) //Z-Displacement
    {
        minVal=zDisp.col(stepIndex-1).minCoeff();
        maxVal=zDisp.col(stepIndex-1).maxCoeff();
    }
    else if (fieldIndex==3) //EPWP
    {
        minVal=Pore.col(stepIndex-1).minCoeff();
        maxVal=Pore.col(stepIndex-1).maxCoeff();
    }
    else if(fieldIndex==4)
    {
        minVal=1;
        maxVal=elements.col(12).maxCoeff();
    }
    else if(fieldIndex==5)
    {
        minVal=1;
        maxVal=elements.col(13).maxCoeff();
    }
    UpdateShownInOpenGLData();
    emit SendBackInformationNewPickedStep(minVal,maxVal);
}

void MainWindow::GetStepResultInformation(BaseStepResult stepResultObject)
{
    this->stepResultObject=stepResultObject;
    UpdateDeformationCoordinates();
    UpdateShownInOpenGLData();

    if(stepResultObject.animationFlag==false)
    {
        viewPortSettingObject.lockView=true;
        Replot();
    }
    else
    {
        viewPortSettingObject.lockView=true;
        RunAnimation();
        stepResultSetting->DoneAnimation();
    }
}

void MainWindow::GetSignalFromOpenGLWidget()
{
    //Send data to OpenGL Widget
    emit SendShownData(coordinatesScale,coordinatesDeform,elements,shownData,elementInfor);
    emit SendViewportInformation(viewPortSettingObject);
    emit SendContourSettingInformation(contourSettingObject);
    emit SendCuttingPlaneInformation(cuttingPlaneSettingObject);
    emit SendStepResultInformation(stepResultObject);

    //disconnect connection
    disconnect(widget,SIGNAL(SendSignalToMainWindow()),this,SLOT(GetSignalFromOpenGLWidget()));
    disconnect(this,SIGNAL(SendShownData(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,ElementInfor)),
               widget,SLOT(GetShownData(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,ElementInfor)));
    disconnect(this,SIGNAL(SendViewportInformation(Base3DViewport)),
               widget,SLOT(GetViewportInformation(Base3DViewport)));
    disconnect(this,SIGNAL(SendContourSettingInformation(BaseContourSetting)),
               widget,SLOT(GetContourSettingInformation(BaseContourSetting)));
    disconnect(this,SIGNAL(SendCuttingPlaneInformation(BaseCuttingPlane)),
               widget,SLOT(GetCuttingPlaneInformation(BaseCuttingPlane)));
    disconnect(this,SIGNAL(SendStepResultInformation(BaseStepResult)),
               widget,SLOT(GetStepResultInformation(BaseStepResult)));

}

void MainWindow::GetValueAtMousePosition(double mouseValue)
{

}

void MainWindow::GetViewportInformation(Base3DViewport viewPortSettingObject)
{
    this->viewPortSettingObject=viewPortSettingObject;
}

void MainWindow::SetUpSignalSlotConnections()
{
    //Connect scale widget to main window
    connect(scaleView,SIGNAL(SendScaleInformation(BaseScaleClass)),this,SLOT(GetScaleInformation(BaseScaleClass)));

    //Connect colour setting to main window
    connect(colorSetting,SIGNAL(SendColorInformation(BaseContourSetting)),this,SLOT(GetColorInformation(BaseContourSetting)));

    //connect cuttingn plane to main window
    connect(cuttingPlaneSetting,SIGNAL(SendClipInformation(BaseCuttingPlane)),this,SLOT(GetCuttingPlaneInformation(BaseCuttingPlane)));

    //connect step result to main window
    connect(stepResultSetting,SIGNAL(PickedNewField(int)),this,SLOT(NewPickedField(int)));
    connect(stepResultSetting,SIGNAL(PickedNewStep(int,int)),this,SLOT(NewPickedStep(int,int)));
    connect(stepResultSetting,SIGNAL(SendStepResultsInformation(BaseStepResult)),this,SLOT(GetStepResultInformation(BaseStepResult)));
    connect(this,SIGNAL(SendBackInformationNewPickedStep(int,bool)),stepResultSetting,SLOT(GetInformationFromMainWindow(int,bool)));
    connect(this,SIGNAL(SendBackInformationNewPickedStep(double,double)),stepResultSetting,SLOT(GetInformationFromMainWindow(double,double)));

    //Connect replot button
    connect(replotButton,SIGNAL(clicked(bool)),this,SLOT(ReplotResetViewport()));

}

void MainWindow::UpdateScaleCoordinates()
{
    coordinatesScale=coordinates;
    if(coordinates.rows()>1)
    {
        coordinatesScale.col(0)=coordinates.col(0);
        coordinatesScale.col(1)=coordinates.col(1)*scaleObject.xScale;
        coordinatesScale.col(2)=coordinates.col(2)*scaleObject.yScale;
        coordinatesScale.col(3)=coordinates.col(3)*scaleObject.zScale;
    }
}

void MainWindow::UpdateDeformationCoordinates()
{
    coordinatesDeform=coordinatesScale;
    int fieldIndex=stepResultObject.fieldIndex;
    int stepIndex=stepResultObject.stepIndex-1;
    vector<int> fieldValid={0,1,2,3};
    auto validCheck=std::find(fieldValid.begin(),fieldValid.end(),fieldIndex);

    if(validCheck != fieldValid.end()) //if results is H, deformation...
    {
        if((xDisp.rows()==yDisp.rows() && yDisp.rows()== zDisp.rows())&&(xDisp.rows()>1))
        {
            coordinatesDeform.col(1)=coordinatesDeform.col(1)+xDisp.col(stepIndex)*scaleObject.deformScale*scaleObject.xScale;
            coordinatesDeform.col(2)=coordinatesDeform.col(2)+yDisp.col(stepIndex)*scaleObject.deformScale*scaleObject.yScale;
            coordinatesDeform.col(3)=coordinatesDeform.col(3)+zDisp.col(stepIndex)*scaleObject.deformScale*scaleObject.zScale;
        }
    }
}

void MainWindow::UpdateShownInOpenGLData()
{
    int fieldIndex=stepResultObject.fieldIndex;
    int stepIndex=stepResultObject.stepIndex-1;
    vector<int> fieldValid={0,1,2,3,4,5};
    auto validCheck=std::find(fieldValid.begin(),fieldValid.end(),fieldIndex);

    bool firstCondition=bool(validCheck != fieldValid.end());
    bool secondCondition=bool((xDisp.rows()==yDisp.rows() && yDisp.rows()== zDisp.rows())&&(zDisp.rows()>1));

    if(firstCondition) //if results is H, deformation...
    {
        if(secondCondition)
        {
            if(fieldIndex==0)
            {
                shownData=xDisp.col(stepIndex);
            }
            else if(fieldIndex==1)
            {
                shownData=yDisp.col(stepIndex);
            }
            else if(fieldIndex==2)
            {
                shownData=zDisp.col(stepIndex);
            }
            else if(fieldIndex==3)
            {
                shownData=Pore.col(stepIndex);
            }         
            else if(fieldIndex==4)
            {
                shownData.resize(1,1);
                shownData(0,0)=1;
            }
            else if(fieldIndex==5)
            {
                shownData.resize(1,1);
                shownData(0,0)=2;
            }
        }
        else
        {
            if(fieldIndex==4)
            {
                shownData.resize(1,1);
                shownData(0,0)=1;
            }
            else if(fieldIndex==5)
            {
                shownData.resize(1,1);
                shownData(0,0)=2;
            }
        }
    }
    else
    {
        QMessageBox::warning(Q_NULLPTR,"ERROR","This is under developed");
        shownData.resize(1,1);
        shownData(0,0)=1;
    }
}

void MainWindow::RunAnimation()
{
    stepResultSetting->close();
    int firstStep=stepResultObject.startStep;
    int lastStep=stepResultObject.endStep;
    if(lastStep<=firstStep)
    {
        QMessageBox::warning(Q_NULLPTR,"ERROR","Chose again start and end step for animation");
    }
    else
    {
        for(int ii=firstStep;ii<lastStep+1;ii++)
        {
            stepResultObject.stepIndex=ii;
            UpdateShownInOpenGLData();
            UpdateDeformationCoordinates();
            Replot();
            QCoreApplication::processEvents();
            QThread::msleep(stepResultObject.delay);
        }
    }
}

void MainWindow::on_actionScale_Setting_triggered()
{
    scaleView->show();
}

void MainWindow::on_actionContour_Colour_Setting_triggered()
{
    colorSetting->show();
}

void MainWindow::on_actionCutting_Plane_Setting_triggered()
{
    cuttingPlaneSetting->show();
}

void MainWindow::on_actionTime_Step_Animation_triggered()
{
    stepResultSetting->show();
}

void MainWindow::Replot()
{
    //Add widget to layout and setting
    if (glLayout->count()!=0)
    {
        glLayout->removeWidget(widget);
        widget->setParent(NULL);
        delete widget;
        widget=NULL;
    }
    widget=new GLWidget;
    QSurfaceFormat format;
    format.setRenderableType(QSurfaceFormat::OpenGL);
    format.setProfile(QSurfaceFormat::OpenGLContextProfile::CompatibilityProfile);
    format.setVersion(3,3);
    widget->setFormat(format);

    //Signal and slot connections
    connect(widget,SIGNAL(SendSignalToMainWindow()),this,SLOT(GetSignalFromOpenGLWidget()));
    connect(widget,SIGNAL(SendValueAtMousePosition(double)),this,SLOT(GetValueAtMousePosition(double)));
    connect(widget,SIGNAL(SendViewportInformation(Base3DViewport)),this,SLOT(GetViewportInformation(Base3DViewport)));

    connect(this,SIGNAL(SendShownData(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,ElementInfor)),
            widget,SLOT(GetShownData(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,ElementInfor)));
    connect(this,SIGNAL(SendViewportInformation(Base3DViewport)),widget,SLOT(GetViewportInformation(Base3DViewport)));
    connect(this,SIGNAL(SendContourSettingInformation(BaseContourSetting)),widget,SLOT(GetContourSettingInformation(BaseContourSetting)));
    connect(this,SIGNAL(SendCuttingPlaneInformation(BaseCuttingPlane)),widget,SLOT(GetCuttingPlaneInformation(BaseCuttingPlane)));
    connect(this,SIGNAL(SendStepResultInformation(BaseStepResult)),widget,SLOT(GetStepResultInformation(BaseStepResult)));


    connect(topButton,SIGNAL(clicked(bool)),widget,SLOT(ZpositiveToXY()));
    connect(botButton,SIGNAL(clicked(bool)),widget,SLOT(ZnegativeToXY()));
    connect(rightButton,SIGNAL(clicked(bool)),widget,SLOT(XpositiveToYZ()));
    connect(leftButton,SIGNAL(clicked(bool)),widget,SLOT(XnegativeToYZ()));
    connect(behindButton,SIGNAL(clicked(bool)),widget,SLOT(YnegativeToXZ()));
    connect(frontButton,SIGNAL(clicked(bool)),widget,SLOT(YpositiveToXZ()));

    glLayout->addWidget(widget);
    widget->ForceToUpdateWiget();
}

void MainWindow::ReplotResetViewport()
{
    viewPortSettingObject.lockView=false;
    Replot();
}
