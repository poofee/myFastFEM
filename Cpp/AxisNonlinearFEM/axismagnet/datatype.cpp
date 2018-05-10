#include "datatype.h"

CMaterial::CMaterial()
    :Bdata(nullptr)
    ,Hdata(nullptr)
    ,BHpoints(0)
{
}

double CMaterial::getMu(double B)
{
    //线性材料
    if(BHpoints == 0){
        return mu;
    }else{
        //通过插值计算非线性材料mu
    }
    return mu;
}

double CMaterial::getdvDb(double B)
{
    //计算非线性材料dvdB

    return 0;
}
