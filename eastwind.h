#ifndef EASTWIND_H
#define EASTWIND_H
#include <cmath>

#endif // EASTWIND_H



class eastwind
{
private:
double H;
public:
eastwind (double H)
{
    this->H = H;
}
double WSol ()
{
    double W;
    double HW[17]={0.1 , 6, 12 , 18 , 24 , 30 , 36, 42 , 48, 54, 60, 66, 72, 78, 84, 90, 94};
    double WW[17]={1.9 , 8, 15 , 18 , 24 , 12 , 10, 8  ,  6,  5,  4,  3,  3,  2,  2,  1,  0.1};
    for (int i = 0; i < 17; i++)
    {
        if (H >= 0.1 && H < 6)  W = 1.9;
        if (H >= 6 && H < 12)  W = 8;
        if (H >= 12 && H < 18)  W = 15;
        if (H >= 18 && H < 24)  W = 18;
        if (H >= 24 && H < 30)  W = 24;
        if (H >= 30 && H < 36)  W = 12;
        if (H >= 36 && H < 42)  W = 10;
        if (H >= 42 && H < 48)  W = 8;
        if (H >= 48 && H < 54)  W = 6;
        if (H >= 54 && H < 60)  W = 5;
        if (H >= 60 && H < 66)  W = 4;
        if (H >= 66 && H < 72)  W = 3;
        if (H >= 72 && H < 84)  W = 2;
        if (H >= 84 && H < 94)  W = 1;
        if (H>=94) {W=0;}
        //if (H >= HW[i] && H < HW[i+1])  W = WW[i];
        return W;
    }

}


};


