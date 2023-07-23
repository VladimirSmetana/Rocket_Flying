#ifndef EQUATIONS_H
#define EQUATIONS_H
#include <cmath>

class equations
{
    private:
    double po;
    double S;
    double g, m, I;
    double CX, CY;
    double PENG, alpha;

    double const  R = 6371000;


    public:

    equations (double po, double S, double g, double m, double CX, double CY, double PENG, double alpha)
    {
        this->po = po;
        this->S = S;
        this->g = g;
        this->m = m;
        this->CX = CX;
        this->CY= CY;
        this->PENG = PENG;
        this->alpha = alpha;
    }


    // Баллистические уравнения
    double fdV(double vv, double ii)
    {
        double F_P = (PENG * cos((M_PI * alpha) / 180));
        double F_X = CX * S * po * pow(vv, 2) / 2;
        return  F_P/m -  F_X/m - g * sin(ii);
    }



    double fdY(double hh, double vv, double ii)
    {
        double F_P = (PENG * sin((M_PI * alpha) / 180));
        double F_Y = (CY * ((M_PI * alpha) / 180) * S * (po * pow(vv, 2)) / 2);
        return (F_P+ F_Y)/ (m*vv)  - ((g  * cos(ii))) * (1-pow(vv,2)/(g*(R+hh))) /vv;
    }

    double returndV(double vv, double ii)
    {
        double F_P = (PENG * cos((M_PI * alpha) / 180));
        double F_X = CX * S * po * pow(vv, 2) / 2;
        return  -F_P/m -  F_X/m - g * sin(ii);
    }

    double returndY(double hh, double vv, double ii)
    {
        double F_P = (PENG * sin((M_PI * alpha) / 180));
        double F_Y = (CY * ((M_PI * alpha) / 180) * S * (po * pow(vv, 2)) / 2);
        return (-F_P+ F_Y)/ (m*vv)  - ((g  * cos(ii))) * (1-pow(vv,2)/(g*(R+hh))) /vv;
    }

    double fdN(double hh, double vv, double ii)
    {
        return (vv /(R+hh)) * cos(ii);
    }

    double fdH(double vv, double ii)
    {
        return vv * sin(ii);
    }

    // Динамические уравнения



};

#endif // EQUATIONS_H



/*
double dv (double O, double w, double v, double b)
{
    Ru  = PENG;
    double Cvv = po*v_prog*S*CX/m;
    double CvO = (g*m*cos(tetta) - PENG*alpha - Ru*b_prog)/m;
    double Cvw = (PENG*alpha + Ru*b_prog)/m;
    double Cvb = (Ru*alpha)/m;
    return -(Cvv*v+CvO*O+Cvb*b+Cvw*w);
}
double dO (double O, double w, double v, double b)
{
    Ru  = PENG;
    double COO = (Ru - g*m*sin(tetta) + (po*pow(v_prog,2)/2)*S*CY -
                  Ru*b_prog*alpha)/(m*v_prog);
    double COw = (-Ru - (po*pow(v_prog,2)/2)*S*CY + Ru*b_prog*alpha)/(m*v_prog);
    double COv = -(po*v_prog*S*CY*alpha)/(m*v_prog);
    double COb = -Ru/(m*v_prog);
    return -(COO*O+COw*w+COv*v+COb*b);
}
*/
