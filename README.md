# C .vs Python Testing

Test with accuracy and time consuming of geolocation distance calculators

### Requirements

```
gcc, python3
```

### To start

```
$ gcc -shared -o vincenty.so -fPIC vincenty.c
$ gcc -shared -o georef.so -fPIC georef.c
$ python3 main.py
```

### Result

Time:

```
                Function    Distance (m)     Time (s)
            great_circle   864214.494339     7.813234
                geodesic   866455.432910   161.672827
     custom_great_circle   864214.494339     2.767540
         custom_geodesic   866455.432912     9.695861
    C++ Geo Great Circle   866455.432913     1.116327
            C++ Geodesic   866455.432915     1.378029
```
