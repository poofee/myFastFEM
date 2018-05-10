#ifndef DATATYPE_H
#define DATATYPE_H

class CNode{
public:
    double x,y;//横坐标、纵坐标
    double A;//磁势

};

class CElementT3{
public:
    int n[3];//三个顶点编号
    double P[3];
    double Q[3];
    double AREA;
    double Bx,By,B;
    double mu;
    int domain;
    double ydot;
    bool isLinear;
};

class CMaterial{
public:
    CMaterial();
    double getMu(double B);
    double getdvDb(double B);
public:
    double mu;//磁导率
    double * Bdata;
    double * Hdata;
    int BHpoints;//BH曲线的节点数目
    double J;//电流密度
};

#endif // DATATYPE_H


