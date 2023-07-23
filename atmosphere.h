#include <QDebug>
#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H
#include <cmath>

#endif // ATMOSPHERE_H

class atmosphere
    {
    private:
    double H;
    double po;
    double Bett, mstep, pp;
    double tCel, Hmas, vsred, lsred, lamb;
    double B[4];
    double A[5];

    double HT[15]{1     ,  11019,  32000 , 47350 , 71802 , 86152 , 104128, 120000 , 140000, 160000, 200000, 325000, 400000, 600000, 1200000};
    double TT[15]{288.15,  216.65, 228.65, 270.65, 214.65, 186.65, 203.81, 334.417, 559.6,  695.6,  834.4,  941.9,  984.65, 995.9,  1000};
    double TMM[8]{288.15,  216.65, 228.65, 270.65, 214.65, 186.65, 212.0 , 380.60 };
    double dHT[14];

    double Bst [11]{-0.0065, -0.0065, 0.0000 , 0.0010, 0.0028, 0, -0.0028, -0.0020, 0, 0.0030, 0.0110};
    double Hst [11] {1,  11019, 20630, 32000 , 47350 , 51412, 71802 , 86152 , 95411, 104128, 120000};

    double const
        Mc = 28.964420,
        gc = 9.80665,
        ac = 340.294,
        Hpc = 8434.5,
        nc = 25.471 * pow(10, 24),
        pc = 101325.0,
        Tc = 288.15,
        vc = 458.94,
        yc = 12.013,
        nuc = 14.607 * pow(10, -6),
        muc = 17.894 * pow(10, -6),
        lac = 25.343 * pow(10, -3),
        omegac = 6.9193 * pow(10, 9),
        poc = 1.2250,
        Na = 602.257 * pow(10, 24),
        RB = 8314.32, r = 287.05287,
        SOS = 110.4,
        BS = 1.458 * pow(10, -6),
        hi = 1.4,
        b = 0.365 * pow(10, -9),
        R = 6371000,
        Hp = (R * H) / (R + H);


        double Hs;
        double Pap;
        double Mol = Mc, g = gc, T = Tc, Tm = Tc * Mc / Mol, P = pc, n = nc, yyd = yc, a = ac, omega = omegac;
        double t;

        public:
        atmosphere (double H)
        {
            for (int i = 1; i < 15; i++)
            {
                if (H >= HT[i - 1] && H < HT[i]) { T = (TT[i-1])+(H-HT[i-1])*(TT[i]-TT[i-1])/ (HT[i] - HT[i - 1]); };
            }

            for (int i = 1; i < 8; i++)
            {
                if (H >= HT[i - 1] && H < HT[i]) { Tm = (TMM[i - 1]) + (H - HT[i - 1]) * (TMM[i] - TMM[i - 1]) / (HT[i] - HT[i - 1]); };
            }


            if (H<94000)
            {


            /// Введение атмосферы

            /// Коэффициенты полинома молярной массы
            B[0] = 46.9083;  B[1] = -29.71210 * pow(10, -5); B[2] = 12.08693 * pow(10, -10); B[3] = -1.85675 * pow(10, -15);

            /// Молярная масса
            Mol = Mc;


                Bett = (7466 * pow(H, 3) - 1262795028 * pow(H, 2) + 61597340039789 * H - 833732588564247562) * pow(10, -20);
                /// Давление
                if (abs(Bett) < 0.0000001) { pp = log(101325); Hs = H - 0.1; P = exp(pp - (0.434294 * gc / (r * T)) * (H- 0)); }  //P = exp(pp - (0.434294 * gc / (r * T)) * (H - Hs)); }
                if (abs(Bett) >= 0.0000001) { pp = log(101325); Hs = H-0.1; P = exp(pp - (gc * log((Tm + Bett *  (H - 0))  / Tm)) / (Bett * r)); }
                Pap = 101325 * exp(-gc*H*Mc/(RB*T));
                /// Плотность
                po = (P * Mol) / (RB * T);
                /// Концентрация частиц воздуха
                n = 7.243611 * pow(10, 22) * P / T;



            tCel = T - 273.15;
            yyd = po * g;
            Hmas = (RB / Mol) * (T / g);
            a = 20.046796 * sqrt(T);
            vsred = 145.50685 * T / Mol;
            lsred = 2.332376 * pow(10, -5) * T / P;
            omega = 6.238629 * pow(10, 6) * P / sqrt(T * Mol);
            lamb = (2.648151 * pow(10, -3) * pow(T, 3 / 2)) / (T + 245.4 * pow(10, -(12 / T)));




        }
            g = gc * pow(R / (R + H), 2);
        }

        double const PI =  3.1415926535;
        double D = 2.66;
        double S = PI*pow(D,2)/4;


        double get_T()
        {
            return T;
        }

        double get_n()
        {
            return n;
        }
        double get_pressure()
        {
            return P;
        }
        double get_density()
        {
            return po;
        }

        double get_AOG()
        {
            return g;
        }

        double get_SV()
        {
            return a;
        }
    };
