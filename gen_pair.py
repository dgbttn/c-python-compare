from random import random

with open('locations.txt', 'w+') as f:
    for _ in range(100000):
        s = '{} {} {} {}\n'.format(
            (random()*1.999 - 0.9999)*90,
            (random()*1.999 - 0.9999)*90,
            (random()*1.999 - 0.9999)*90,
            (random()*1.999 - 0.9999)*90
        )
        f.write(s)
