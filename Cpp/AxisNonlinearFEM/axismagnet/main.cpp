#include "mainwindow.h"
#include "axismagnet2d.h"
#include "superlutest.h"
#include <QApplication>
#include <QDebug>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    superlumttest();
    AxisMagnet2D fem;
    fem.loadCOMSOLmeshfile("../mesh.mphtxt");
    fem.NewtonSolve();
    return a.exec();
}
