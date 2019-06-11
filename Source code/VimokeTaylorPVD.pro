#-------------------------------------------------
#
# Project created by QtCreator 2019-06-05T10:21:52
#
#-------------------------------------------------

QT       += core gui
CONFIG +=console

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

INCLUDEPATH += "C:/Eigen/Eigen"
INCLUDEPATH +=C:/BoostLib
INCLUDEPATH +=C:/BoostLib/stage/lib
TARGET = VTP
TEMPLATE = app

INCLUDEPATH += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/include"
INCLUDEPATH += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/include/intel64"
INCLUDEPATH += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/lib/intel64"
INCLUDEPATH += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/lib/intel64_win"

LIBS += "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.4.246\windows\mkl\lib\intel64_win\mkl_core.lib"
LIBS += "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.4.246\windows\mkl\lib\intel64_win\mkl_intel_lp64.lib"
LIBS += "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.4.246\windows\mkl\lib\intel64_win\mkl_intel_thread.lib"
LIBS += "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.4.246\windows\compiler\lib\intel64_win\libiomp5md.lib"

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
        main.cpp \
        mainwindow.cpp \
    glwidget.cpp \
    input.cpp \
    consolidation3d.cpp \
    gauss.cpp \
    getfile.cpp \
    glwidget.cpp \
    input.cpp \
    main.cpp \
    mainwindow.cpp \
    pvdlib.cpp \
    writetofile.cpp \
    importmeshdata.cpp \
    soilparameter.cpp \
    pvdparameters.cpp \
    analysistype.cpp \
    progess.cpp \
    importresults.cpp \
    analysistype.cpp \
    consolidation3d.cpp \
    contoursetting.cpp \
    cuttingplane.cpp \
    gauss.cpp \
    getfile.cpp \
    glwidget.cpp \
    importmeshdata.cpp \
    importresults.cpp \
    input.cpp \
    main.cpp \
    mainwindow.cpp \
    progess.cpp \
    pvdlib.cpp \
    pvdparameters.cpp \
    scalesetting.cpp \
    soilparameter.cpp \
    stepresult.cpp \
    writetofile.cpp \
    analysistype.cpp \
    consolidation3d.cpp \
    contoursetting.cpp \
    cuttingplane.cpp \
    gauss.cpp \
    getfile.cpp \
    glwidget.cpp \
    importmeshdata.cpp \
    importresults.cpp \
    input.cpp \
    main.cpp \
    mainwindow.cpp \
    progess.cpp \
    pvdlib.cpp \
    pvdparameters.cpp \
    scalesetting.cpp \
    soilparameter.cpp \
    stepresult.cpp \
    writetofile.cpp \
    base3dviewport.cpp

HEADERS += \
        mainwindow.h \
    glwidget.h \
    input.h \
    consolidation3d.h \
    gauss.h \
    getfile.h \
    glwidget.h \
    input.h \
    mainwindow.h \
    pvdlib.h \
    writetofile.h \
    importmeshdata.h \
    soilparameter.h \
    pvdparameters.h \
    analysistype.h \
    progess.h \
    importresults.h \
    contoursetting.h \
    analysistype.h \
    consolidation3d.h \
    contoursetting.h \
    cuttingplane.h \
    gauss.h \
    getfile.h \
    glwidget.h \
    importmeshdata.h \
    importresults.h \
    input.h \
    mainwindow.h \
    progess.h \
    pvdlib.h \
    pvdparameters.h \
    scalesetting.h \
    soilparameter.h \
    stepresult.h \
    writetofile.h \
    base3dviewport.h \
    elementinformation.h

FORMS += \
        mainwindow.ui \
    importmeshdata.ui \
    soilparameter.ui \
    pvdparameters.ui \
    analysistype.ui \
    progess.ui \
    importresults.ui \
    analysistype.ui \
    contoursetting.ui \
    cuttingplane.ui \
    importmeshdata.ui \
    importresults.ui \
    mainwindow.ui \
    progess.ui \
    pvdparameters.ui \
    scalesetting.ui \
    soilparameter.ui \
    stepresult.ui \
    analysistype.ui \
    contoursetting.ui \
    cuttingplane.ui \
    importmeshdata.ui \
    importresults.ui \
    mainwindow.ui \
    progess.ui \
    pvdparameters.ui \
    scalesetting.ui \
    soilparameter.ui \
    stepresult.ui

RESOURCES += \
    resources.qrc
