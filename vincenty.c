// #ifdef __MSVC__ //Patch from Nohal to make Route compile under MSVC
// #define _USE_MATH_DEFINES
#include <math.h>

const double InverseFlattening = 298.25722356366546517369570015525;
const double Flattening = 0.003352810664739999495231881355;
const double SemiMajorAxis = 6378137.0;
const double SemiMinorAxis = 6356752.314245226792991161346435546875;
const double Eps = 0.5e-15;

double DistVincenty(double, double, double, double);

double toRad(double degree)
{
    return degree * M_PI / 180;
}

double DistVincenty(double lat1, double lon1, double lat2, double lon2)
{
    double L = toRad(lon2) - toRad(lon1);
    double U1 = atan((1 - Flattening) * tan(toRad(lat1)));
    double U2 = atan((1 - Flattening) * tan(toRad(lat2)));

    double sinU1 = sin(U1);
    double cosU1 = cos(U1);
    double sinU2 = sin(U2);
    double cosU2 = cos(U2);

    double dCosU1CosU2 = cosU1 * cosU2;
    double dCosU1SinU2 = cosU1 * sinU2;

    double dSinU1SinU2 = sinU1 * sinU2;
    double dSinU1CosU2 = sinU1 * cosU2;

    double lambda = L;
    double lambdaP = (M_PI * 2);
    int iterLimit = 0;
    double cosSqAlpha;
    double sinSigma;
    double cos2SigmaM;
    double cosSigma;
    double sigma;
    double sinAlpha;
    double C;
    double sinLambda, cosLambda;

    do
    {
        sinLambda = sin(lambda);
        cosLambda = cos(lambda);
        sinSigma = sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) +
                        (dCosU1SinU2 - dSinU1CosU2 * cosLambda) * (dCosU1SinU2 - dSinU1CosU2 * cosLambda));

        if (sinSigma == 0)
        {
            return -1;
        }
        cosSigma = dSinU1SinU2 + dCosU1CosU2 * cosLambda;
        sigma = atan2(sinSigma, cosSigma);
        sinAlpha = dCosU1CosU2 * sinLambda / sinSigma;
        cosSqAlpha = 1.0 - sinAlpha * sinAlpha;

        if (cosSqAlpha < 1e-8)
        {
            cos2SigmaM = 0.0; // equatorial line: cosSqAlpha=0 (§6)
        }
        else
        {
            cos2SigmaM = cosSigma - 2.0 * dSinU1SinU2 / cosSqAlpha;
        }

        C = Flattening / 16.0 * cosSqAlpha * (4.0 + Flattening * (4.0 - 3.0 * cosSqAlpha));
        lambdaP = lambda;
        lambda = L + (1.0 - C) * Flattening * sinAlpha *
                         (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM)));
    } while (fabs(lambda - lambdaP) > Eps && ++iterLimit < 100);

    double uSq = cosSqAlpha * (SemiMajorAxis * SemiMajorAxis - SemiMinorAxis * SemiMinorAxis) /
                 (SemiMinorAxis * SemiMinorAxis);
    double A = 1.0 + uSq / 16384.0 * (4096.0 + uSq * (-768.0 + uSq * (320.0 - 175.0 * uSq)));
    double B = uSq / 1024.0 * (256.0 + uSq * (-128.0 + uSq * (74.0 - 47.0 * uSq)));
    double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4.0 * (cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM) - B / 6.0 * cos2SigmaM * (-3.0 + 4.0 * sinSigma * sinSigma) * (-3.0 + 4.0 * cos2SigmaM * cos2SigmaM)));

    if (iterLimit < 98)
    {
        return SemiMinorAxis * A * (sigma - deltaSigma);
    }
    else
    {
        return -1.0;
    } //if false function failed to converge!
}