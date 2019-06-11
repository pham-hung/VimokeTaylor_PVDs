#ifndef ELEMENTINFORMATION_H
#define ELEMENTINFORMATION_H

#include <Eigen/Dense>
using namespace Eigen;

struct ElementInfor{

    MatrixXi elmementIndex=MatrixXi::Zero(1,1);
    int tet10Count=0, hex8Count=0, hex20Count=0, pyra5Count=0, pyra13Count=0, prism6Count=0, prism15Count=0;

    void resetCountNumber()
    {
        elmementIndex=MatrixXi::Zero(0,0);
        tet10Count=0;
        hex8Count=0;
        hex20Count=0;
        pyra5Count=0;
        pyra13Count=0;
        prism6Count=0;
        prism15Count=0;
    }
};
#endif // ELEMENTINFORMATION_H
