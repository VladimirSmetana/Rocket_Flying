#ifndef ALPHA_H
#define ALPHA_H
#include <cmath>

class alpha
{
private:
    double Vel, k1, k2, time, septime, def, tv1, tv2;
public:
    alpha (double Vel, double k1, double k2, double time, double septime, int def, int tv1, int tv2)
    {
        this->Vel = Vel;
        this->k1 = k1;
        this->k2 = k2;
        this->time = time;
        this->septime = septime;
        this->def = def;
        this->tv1 = tv1;
        this->tv2 = tv2;
    }

    double A ()
    {
        double ans;
        if (Vel>50 && Vel<270 && time <= septime)
        {
            double z=M_PI*(Vel-50);
            double che=(Vel-50)+0.25*(270-Vel);
            ans = k1*pow(sin(z/che),2);
        }
        else
        if (time >= septime)
        {

            if (def == 0) ans = k2 * (time - septime);
            if (def == 180) {
                if (time - septime<60)
                {
                    double z=M_PI*(time - septime);
                    double che=(time - septime)+0.25*(septime+60-time);
                    ans = 90*pow(sin(z/che),2);
                    //ans = 3 * (time - septime);
                }
                if (time>tv1 && time<tv2)
                {
                double z=M_PI*(time - tv1);
                double che=(time - tv1)+0.25*(tv1+(tv2-tv1)-time);
                ans = 0.25*pow(sin(z/che),2); }

       //         }

                 } ;

        }
        else
        ans = 0;


        return ans;
    }
};

#endif // ALPHA_H
