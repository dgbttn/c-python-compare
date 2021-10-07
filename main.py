import time
from ctypes import CDLL, c_double
import pathlib
from math import sin, cos, asin, sqrt, radians, tan, atan2, pi
from geopy.distance import great_circle, geodesic


N = 1000000


class custom_great_circle(great_circle):
    """Haversine formula"""

    RADIUS = 6371.009

    def measure(self, a, b):
        lat1, lng1 = radians(a[0]), radians(a[1])
        lat2, lng2 = radians(b[0]), radians(b[1])

        delta_lat = lat2 - lat1
        delta_long = lng2 - lng1

        haversine = sin(delta_lat/2)**2 + cos(lat1) * \
            cos(lat2)*sin(delta_long/2)**2
        return asin(sqrt(haversine))*2*self.RADIUS


class custom_geodesic(geodesic):
    """Vincenty formula"""

    ELLIPSOID = (6378.137, 6356.7523142, 1 / 298.257223563)

    def measure(self, a, b):
        lat1, lng1 = radians(a[0]), radians(a[1])
        lat2, lng2 = radians(b[0]), radians(b[1])

        a, b, f = self.ELLIPSOID
        L = lng2 - lng1

        u1 = atan2((1-f) * tan(lat1), 1)
        u2 = atan2((1-f) * tan(lat2), 1)

        sinU1, cosU1 = sin(u1), cos(u1)
        sinU2, cosU2 = sin(u2), cos(u2)

        sinU1sinU2 = sinU1 * sinU2
        sinU1cosU2 = sinU1 * cosU2
        cosU1cosU2 = cosU1 * cosU2
        cosU1sinU2 = cosU1 * sinU2

        lamda, lamdaP, iter_limit = L, 2*pi, 100
        while iter_limit > 0 and abs(lamda-lamdaP) > 1e-12:
            sinLambda, cosLambda = sin(lamda), cos(lamda)
            sinSigma = sqrt((cosU2*sinLambda)**2 +
                            (cosU1sinU2-sinU1cosU2*cosLambda)**2)

            if sinSigma == 0:
                return -1

            cosSigma = sinU1sinU2 + cosU1cosU2*cosLambda
            sigma = atan2(sinSigma, cosSigma)
            sinAlpha = cosU1cosU2*sinLambda/sinSigma
            cosSqAlpha = 1.0 - sinAlpha**2

            if cosSqAlpha < 1e-15:
                cos2SigmaM = 0
            else:
                cos2SigmaM = cosSigma - 2.0 * sinU1sinU2/cosSqAlpha

            C = f/16.0 * cosSqAlpha * (4.0 + f*(4.0 - 3.0*cosSqAlpha))
            lamdaP = lamda
            lamda = L + (1.0-C)*f*sinAlpha*(
                sigma + C*sinSigma * (cos2SigmaM + C*cosSigma*(-1.0 + 2.0*cos2SigmaM**2)))
            iter_limit -= 1

        if iter_limit == 0:
            return -1

        uSq = cosSqAlpha * ((a/b)**2 - 1.0)
        A = 1.0 + uSq / 16384.0 * \
            (4096.0 + uSq * (-768.0 + uSq * (320.0 - 175.0 * uSq)))
        B = uSq / 1024.0 * (256.0 + uSq * (-128.0 + uSq * (74.0 - 47.0 * uSq)))
        deltaSigma = B * sinSigma * (
            cos2SigmaM + B / 4.0 * (
                cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM) -
                B / 6.0 * cos2SigmaM * (-3.0 + 4.0 * sinSigma * sinSigma) * (-3.0 + 4.0 * cos2SigmaM * cos2SigmaM)))
        return b * A * (sigma - deltaSigma)


def test_dist_func(func):
    p1 = (41.49008, -71.312796)
    p2 = (41.499498, -81.695391)
    d = 0
    now = time.time()
    for _ in range(N):
        if func.__name__ in ('C++ Geodesic', 'C++ Geo Great Circle'):
            d = func(*p1, *p2)
        else:
            d = func(p1, p2).meters
    print('%30s %13f %10f' % (func.__name__, d, time.time()-now))


def load_cdll(filename, funcname):
    cfile = CDLL('{}/{}.so'.format(pathlib.Path().absolute(), filename))
    func = cfile.__getattr__(funcname)
    func.argtypes = [c_double, c_double, c_double, c_double]
    func.restype = c_double
    return func


vincenty_func = load_cdll('vincenty', 'DistVincenty')
vincenty_func.__name__ = 'C++ Geodesic'
georef_func = load_cdll('georef', 'DistGreatCircle')
georef_func.__name__ = 'C++ Geo Great Circle'


if __name__ == '__main__':
    print('%30s %13s %10s' % ('Function', 'Distance (m)', 'Time (s)'))
    test_dist_func(great_circle)
    test_dist_func(geodesic)
    test_dist_func(custom_great_circle)
    test_dist_func(custom_geodesic)
    test_dist_func(georef_func)
    test_dist_func(vincenty_func)

    # with open('locations.txt', 'r') as f:
    #     locations = [list(map(float, s.split())) for s in f.readlines()]

    # funcs = [great_circle, custom_great_circle,
    #          georef_func, custom_geodesic, vincenty_func]
    # max_diff = [0.0, 0.0, 0.0, 0.0, 0.0]
    # p_max_diff = [0.0, 0.0, 0.0, 0.0, 0.0]
    # sum_diff = [0.0, 0.0, 0.0, 0.0, 0.0]
    # p_sum_diff = [0.0, 0.0, 0.0, 0.0, 0.0]

    # print(len(locations))
    # for loc in locations:
    #     p1 = (loc[0], loc[1])
    #     p2 = (loc[2], loc[3])
    #     dist = geodesic(p1, p2).meters
    #     for i, func in enumerate(funcs):
    #         if func.__name__ in ('C++ Geodesic', 'C++ Geo Great Circle'):
    #             d = func(*loc)
    #         else:
    #             d = func(p1, p2).meters
    #         diff = abs(dist-d)
    #         max_diff[i] = max(max_diff[i], diff)
    #         sum_diff[i] += diff
    #         p_max_diff[i] = max(p_max_diff[i], diff/dist*100.0)
    #         p_sum_diff[i] += diff/dist*100.0

    # print(max_diff)
    # print(*map(lambda s: s/len(locations), sum_diff))
    # print(p_max_diff)
    # print(*map(lambda s: s/len(locations), p_sum_diff))

    # # GC , custom GC, C++ GC, custom Geo, C++ Geo
