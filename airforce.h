#include <QDebug>
#ifndef AIRFORCE_H
#define AIRFORCE_H
#include <cmath>

#endif // AIRFORCE_H

class airforce
{
    private:
    double CX;
    double CY;
    public:
    airforce (double Mah)
    {

        if ((Mah > 0) && (Mah <= 0.25)) CY = 2.8;
        if ((Mah > 0.25) && (Mah <= 1.1)) CY = 2.8+0.447*(Mah-0.25);
        if ((Mah > 1.1) && (Mah <= 1.6)) CY = 3.18 - 0.660 * (Mah - 1.1);
        if ((Mah > 1.6) && (Mah <= 3.6)) CY = 2.85 + 0.350 * (Mah - 1.6);
        if (Mah > 3.6) CY = 12.78/Mah;

        if (Mah <= 0.8) { CX = 0.29; };
        if (Mah > 0.8 && Mah < 1.068) { CX = Mah-0.51; };
        if (Mah >= 1.068) { CX = 0.597/Mah; };
        if (CX > 1) { CX = 1; };
    }



    double getCX()
    {
        return CX;
    }

    double getCY()
    {
        return CY;
    }
};
