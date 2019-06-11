#ifndef IMPORTRESULTS_H
#define IMPORTRESULTS_H

#include <QDialog>
#include <QFileDialog>
#include <QString>

using namespace std;

namespace Ui {
class ImportResults;
}

class ImportResults : public QDialog
{
    Q_OBJECT

public:
    explicit ImportResults(QWidget *parent = 0);
    ~ImportResults();    
signals:
    void SendResultFiles(QString xDispFile,QString yDispFile, QString zDispFile, QString poreFile);
private slots:
    void on_closeButton_clicked();
    void on_startImportButton_clicked();
    void on_browseXButton_clicked();
    void on_browseYButton_clicked();
    void on_browseZButton_clicked();
    void on_browsePButton_clicked();

private:
    Ui::ImportResults *ui;
    QString xDispFile,yDispFile,zDispFile,poreFile;

};

#endif // IMPORTRESULTS_H
