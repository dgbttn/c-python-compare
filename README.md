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
