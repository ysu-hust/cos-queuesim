#!/usr/bin/env python

#*****************************************************************************80
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 January 2018
#
#  Author:
#
#    Su Yi (suyi@hust.edu.cn)

import random

class zipfian_generator(object):
    def __init__(self, item_min, item_max, zipfianconstant=0.99):
        self.lastVal = None
        self.items = item_max - item_min + 1
        self.base = item_min
        self.zipfianconstant = zipfianconstant
        self.theta = self.zipfianconstant
        self.alpha = 1.0 / (1.0 - self.theta)
        self.zetan = self.zetastatic(self.items, self.zipfianconstant)
        self.zeta2theta = self.zeta(2, self.theta)
        self.countforzeta = self.items
        self.eta = (1 - (2.0/self.items)**(1-self.theta)) / (1 - (self.zeta2theta/self.zetan))
        self.allowitemcountdecrease = True
        self.__next__(self.items)

    def zeta(self, n, thetaVal, st=0, initialsum=0):
        self.countforzeta = n
        return self.zetastatic(n, thetaVal, st=st, initialsum=initialsum)

    def zetastatic(self, n, theta, st=0, initialsum=0):
        s = initialsum
        for i in range(st, n):
            s += 1.0 / ((i+1)**theta)
        return s

    def __next__(self, itemcount):
        if itemcount != self.countforzeta:
            if itemcount > self.countforzeta:
                self.zetan = self.zeta(itemcount, self.theta, st=self.countforzeta, initialsum=self.zetan)
                self.eta = (1 - (2.0 / self.items)**(1 - self.theta)) / (1 - self.zeta2theta / self.zetan)
            elif itemcount < self.countforzeta and self.allowitemcountdecrease:
                self.zetan = self.zeta(itemcount, self.theta)
                self.eta = (1 - (2.0 / self.items)**(1 - self.theta)) / (1 - self.zeta2theta / self.zetan)
            self.countforzeta = itemcount
        u = random.random()
        uz = u * self.zetan
        if uz < 1.0:
            return self.base
        if uz < 1.0 + 0.5**self.theta:
            return self.base+1
        ret = self.base + long(itemcount * (self.eta * u - self.eta + 1)**self.alpha)
        self.setLastValue(ret)
        return ret

    def next(self, itemcount):
        return self.__next__(itemcount)

    def next_same(self):
        return self.__next__(self.countforzeta)

    def setLastValue(self, last):
        self.lastVal = last


if __name__ == "__main__":
    zipf = zipfian_generator(0, 1000000)
    item_counts = {}
    for i in range(1000000):
        item = zipf.next()
        item_counts[item] = item_counts.get(item, 0) + 1
    sorted_ic = sorted(item_counts.items(), key=lambda x:x[1])
    for i in range(1, 10):
        print sorted_ic[i], sorted_ic[-i]
