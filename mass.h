#ifndef MACENIN_H
#define MACENIN_H
#include <cmath>

class  mass
{
    private:
    double Ratio = 3.5;
    double po_dry = 300;

    double po_C = 440, po_O = 1150;
    double s[2];

    double mb[2];
    double m_dry[2];
    double m_fuel[2];
    double m_instrument;
    double M_Rocket;

    double Peng[2];

    double w_dry[2];
    double w_headunit;
    double w_instrument;

    double m_C[2];
    double m_O[2];


    double L_C[2], L_O[2];
    double U_C[2], U_O[2];

    double L_spC[2], L_spO[2];
    double L_b[2];

    double L_headunit = 11;
    double L_instrument = 1;

    double W_Rocket;
    double L_Rocket;

    double CIL, CON, CILCON;

    double S_dry[2], S_o[2], S_c[2], S_instrument, S_headunit, S_Rocket;
    double I_dry[2], I_o[2], I_c[2], I_instrument, I_headunit, I_Rocket;

    double m_center;



    public:
    double L_taC[2], L_taO[2];
    double W_taC[2], W_taO[2];
    double w_C[2], w_O[2];
    double L_engine[2];
    double K[14];

    double fun_S(double X_l, double X_b, double mass)
                 {
                     return 0.5 * mass * (X_l+X_b);
                 }
    double fun_I(double X_l, double X_b, double mass, double D)
                 {
                     return 0.25 * mass * ( pow(X_l+X_b, 2) + 0.333 * pow(X_b - X_l, 2)+pow((D/2),2) );
                 }

    void MCI_f
    (
            double time,double h,double payload, double D,
            double m_block1, double m_block2,
            double s1, double s2,
            double P1, double P2
    )

    {
        mb[0] = m_block1; mb[1] = m_block2;
        M_Rocket = mb[0] + mb[1] + payload;
        s[0] = s1; s[1] = s2;
        Peng[0] = P1; Peng[1] = P2;

        w_headunit = 80;

        m_instrument = 0.0008*M_Rocket;
        w_instrument = m_instrument/150;
        W_Rocket = w_headunit+w_instrument;

        L_Rocket += L_headunit + L_instrument;

        for (int i=0;i<=1;i++)
        {
            m_fuel[i] = mb[i] * (s[i] - 1) / s[i];
            m_dry[i] = mb[i]-m_fuel[i];
            w_dry[i]=m_dry[i]/po_dry;
            m_C[i]=m_fuel[i]*1/(1+Ratio);
            m_O[i]=m_fuel[i]*Ratio/(1+Ratio);
            w_C[i] = m_C[i]/po_C;
            w_O[i] = m_O[i]/po_O;

            //L_C[i] = w_C[i]*4/(M_PI*pow(D,2));
            //L_O[i] = w_O[i]*4/(M_PI*pow(D,2));
            L_taC[i] = 1.1*L_C[i];
            L_taO[i] = 1.1*L_O[i];
            //L_spC[i] = L_taC[i] - L_C[i];
            //L_spO[i] = L_taO[i] - L_O[i];


            L_b[i] = L_taC[i]+L_taO[i]+L_engine[i];

            W_taC[i] = 1.25*w_C[i];
            W_taO[i] = 1.25*w_O[i];
            W_Rocket +=  W_taC[i] + W_taO[i];
        }
            W_Rocket = W_Rocket/(1-0.15);


            //L_engine[0] = 3.2;
            //L_engine[1] = 5.4;



        //U_C[0] = time -

        // Геометрические исходные данные

        L_spO[1] = 0.6;
        L_O[1] = 4.8;
        L_spC[1] = 0.6;
        L_C[1] = 4.1+0.6;
        L_engine[1] =3.2-0.6;
        L_spO[0] = 1;
        L_O[0] = 17.3;
        L_spC[0] = 2.7;
        L_C[0] = 13.2;
        L_engine[0] = 4.9;
        // Координаты граничных точек РН
        K[0] = 0;
        K[1] = K[0] + L_headunit;
        K[2] = K[1] + L_instrument;
        K[3] = K[2] + L_spO[1];
        K[4] = K[3] + L_O[1];
        K[5] = K[4] + L_spC[1];
        K[6] = K[5] + L_C[1];
        K[7] = K[6] + L_engine[1];
        K[8] = K[7] + L_spO[0];
        K[9] = K[8] + L_O[0];
        K[10]= K[9] + L_spC[0];
        K[11]= K[10]+ L_C[0];
        K[13]= K[11]+ 1.3;
        K[12]= K[11]+ L_engine[0];



        L_Rocket =K[12];


        CIL = (L_Rocket - L_headunit) / D; CON = (L_headunit) / D; CILCON = CIL / CON;
        // Определение моментов инерции
        S_dry[0] = fun_S (K[6], K[12], m_dry[0]);
        S_o[0] = fun_S (K[8], K[9], m_O[0]);
        S_c[0] = fun_S (K[10], K[11], m_C[0]);
                S_dry[1] = fun_S (K[1], K[6], m_dry[1]);
                S_o[1] = fun_S (K[3], K[4], m_O[1]);
                S_c[1] = fun_S (K[5], K[6], m_C[1]);
        S_headunit = fun_S (K[0], K[1], 16500);
        //S_Rocket = S_instrument + S_Rocket;

        I_dry[0] = fun_I (K[6], K[12], m_dry[0], D);
        I_o[0] = fun_I (K[10], K[11], m_O[0], D);
        I_c[0] = fun_I (K[8], K[9], m_C[0], D);
                I_dry[1] = fun_I (K[1], K[6], m_dry[1], D);
                I_o[1] = fun_I (K[5], K[6], m_O[1], D);
                I_c[1] = fun_I (K[3], K[4], m_C[1], D);
        I_headunit = fun_I (K[0], K[1], 16500, D);


    m_center = S_Rocket/M_Rocket;



    return;
    }

    double Focus (double Mah, double D, double Lgo, double wgo, double CILCON)
    {
        double dXF;
        double KM = (pow(abs(1 - pow(Mah, 2)), 0.5)) / CIL;
        if (CILCON > 0.0 && CILCON <= 0.5) { if (Mah < 1) { dXF = -0.03 * pow(KM, 2) - 0.04 * KM + 0.09; } if (Mah >= 1) { dXF = -0.03 * pow(KM, 2) + 0.06 * KM + 0.09; }; }
        if (CILCON > 0.5 && CILCON <= 1.0) { if (Mah < 1) { dXF = -0.03 * pow(KM, 2) - 0.04 * KM + 0.09; } if (Mah >= 1) { dXF = -0.03 * pow(KM, 2) + 0.16 * KM + 0.10; }; }
        if (CILCON > 1.0 && CILCON <= 2.0) { if (Mah < 1) { dXF = -0.03 * pow(KM, 2) - 0.04 * KM + 0.09; } if (Mah >= 1) { dXF = -0.05 * pow(KM, 2) + 0.34 * KM + 0.05; }; }
        if (CILCON > 2.0 && CILCON <= 3.0) { if (Mah < 1) { dXF = -0.03 * pow(KM, 2) - 0.04 * KM + 0.09; } if (Mah >= 1) { dXF = -0.07 * pow(KM, 2) + 0.42 * KM + 0.01; }; }
        if (CILCON > 3.0)                  { if (Mah < 1) { dXF = -0.03 * pow(KM, 2) - 0.04 * KM + 0.09; } if (Mah >= 1) { dXF = -0.08 * pow(KM, 2) + 0.51 * KM - 0.02; }; }
        double XF = Lgo - wgo / ((M_PI * pow(D, 2)) / 4) + Lgo * dXF;
        return XF;
    }

    double get_CM()
    {
        return m_center;
    }
    double get_I()
    {
        return I_Rocket;
    }
    double get_S()
    {
        return S_Rocket;
    }
    double get_lenght()
    {
        return L_Rocket;
    }
    double get_lst1()
    {
        return L_b[0];
    }
    double get_lst2()
    {
        return L_b[1];
    }
    double get_SGO()
    {
        return S_headunit;
    }
    double get_IGO()
    {
        return I_headunit;
    }

    double get_wgo()
    {
        return w_headunit;
    }

    double get_lengo()
    {
        return L_headunit;
    }

    double get_CIL()
    {
        return CIL;
    }



};



#endif // MACENIN_H
