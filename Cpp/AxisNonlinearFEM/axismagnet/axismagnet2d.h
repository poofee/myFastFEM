#ifndef AXISMAGNET2D_H
#define AXISMAGNET2D_H
#include "datatype.h"

class AxisMagnet2D
{
public:
    AxisMagnet2D();

    bool loadCOMSOLmeshfile(const char fn[]);
    bool NewtonSolve();

public:
    int num_point;
    int num_element;
    int num_domain;
    CNode * pnode;
    CElementT3 * pelement;
    char filename[256];
    CMaterial* materialList;//定义材料列表
    double tolerance;//求解精度
};

#endif // AXISMAGNET2D_H
