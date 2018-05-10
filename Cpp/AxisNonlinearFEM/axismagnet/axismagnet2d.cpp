#include "axismagnet2d.h"
#include <stdio.h>
#include <QDebug>

AxisMagnet2D::AxisMagnet2D()
    :pnode(nullptr)
    ,pelement(nullptr)
    ,materialList(nullptr)
    ,num_point(0)
    ,num_domain(0)
    ,num_element(0)
    ,tolerance(1e-7)
{

}

AxisMagnet2D::~AxisMagnet2D()
{
    if(pnode != nullptr) free(pnode);
    if(pelement != nullptr) free(pelement);
    if(materialList != nullptr) free(materialList);
}

bool AxisMagnet2D::loadCOMSOLmeshfile(const char fn[])
{
    qDebug()<<"Load COMSOL MESH file Starts!";

    qDebug()<<"Successfully Loaded COMSOL MESH file!";
    return true;
}

bool AxisMagnet2D::NewtonSolve()
{
    qDebug()<<"Newton Solve Starts!";

    qDebug()<<"Newton Solve Ends!";
    return true;
}
