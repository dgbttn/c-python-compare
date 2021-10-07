import time
import pathlib
from ctypes import CDLL, c_double
from geopy.distance import great_circle, geodesic


N = 100000


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
    func = getattr(cfile, funcname)
    func.argtypes = [c_double, c_double, c_double, c_double]
    func.restype = c_double
    return func


vincenty_func = load_cdll('vincenty', 'DistVincenty')
vincenty_func.__name__ = 'C++ Geodesic'
georef_func = load_cdll('georef', 'DistGreatCircle')
georef_func.__name__ = 'C++ Geo Great Circle'


def test_time_consuming():
    print('%30s %13s %10s' % ('Function', 'Distance (m)', 'Time (s)'))
    test_dist_func(great_circle)
    test_dist_func(geodesic)
    test_dist_func(georef_func)
    test_dist_func(vincenty_func)


if __name__ == '__main__':
    test_time_consuming()
