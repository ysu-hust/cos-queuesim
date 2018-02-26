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

import time
import random
from cache import LRUCache
from zipfian_generator import zipfian_generator
from predictor import *


def get_sertime(ser_rate, ser_type="deter"):
    # ser_type = "deter"
    # ser_type = "exp"
    # ser_type = "gamma"
    # ser_type = "cache"
    if ser_type == "deter":
        return 1.0/ser_rate
    elif ser_type == "mix":
        p = random.random()
        if p > 0.8:
            return 4.6/ser_rate
        else:
            return 0.1/ser_rate
    elif ser_type == "exp":
        return random.expovariate(ser_rate)
    elif ser_type == "gamma":
        alpha = 4.0
        alpha = 3.0
        #alpha = 1.0
        beta = 1.0/(ser_rate*alpha)
        return random.gammavariate(alpha, beta)
    elif ser_type == "cache":
        p = random.random()
        if p > 0.7:
            return 10.0/ser_rate
        else:
            return 1.0/ser_rate
    else:
        return 0


class Req(object):
    def __init__(self, obj_name=None):
        self.obj_name = obj_name
        self.st = None
        self.et = None
        self.ser_t = None
    
    def duration(self):
        if self.et is not None and self.st is not None:
            return self.et - self.st
        else:
            return None

    def start(self, st):
        self.st = st

    def finish(self, et):
        self.et = et


class simpleQueue(object):
    def __init__(self, service_rate, get_st=lambda ser_rate: get_sertime(ser_rate, ser_type="deter")):
        self.srate = service_rate
        self.cur_t = 0
        self.get_st = get_st

    def get_stime(self):
        return self.get_st(self.srate)

    def arrive(self, t):
        # print t, self.cur_t
        if self.cur_t < t:
            self.cur_t = t
        self.cur_t += self.get_stime()
        # print "\t", t, self.cur_t
        return self.cur_t

    def clear(self):
        self.cur_t = 0


class Que(object):
    def __init__(self, cache=None, generator=None, successor=None, get_st=lambda ser_rate: get_sertime(ser_rate, ser_type="deter")):
        self.q = []
        self.cur_t = {}
        self.res = []
        self.cache = cache
        self.generator = generator
        self.successor = successor
        self.ac_list = {}
        self.get_st = get_st

    def start(self, n, rate):
        cur_t = 0
        for i in range(n):
            t = random.expovariate(rate)
            cur_t += t
            if self.generator is None:
                self.q.append(cur_t)
            else:
                obj_name = self.generator.next_same()
                self.q.append((cur_t, obj_name))

    def warmup(self, n):
        for i in range(n):
            obj_name = self.generator.next_same()
            if self.cache.get(obj_name) < 0:
                self.cache.set(obj_name, 1)

    def proc(self, n, ser_rate):
        for i in range(n):
            self.cur_t[i] = 0
        for i in self.q:
            if self.generator is not None:
                i, obj_name = i
                req = Req(obj_name=obj_name)
            else:
                req = Req()
            req.start(i)
            st = self.get_st(ser_rate)
            cur_ser = -1
            cur_time = float('inf')
            for j in range(n):
                if self.cur_t[j] < cur_time:
                    cur_time = self.cur_t[j]
                    cur_ser = j
            #print i, cur_ser, cur_time
            if cur_time < i:
                cur_time = i
            self.cur_t[cur_ser] = cur_time + st
            req.finish(self.cur_t[cur_ser])
            req.ser_t = st
            #print self.cur_t
            self.res.append(req)

    def cache_proc(self, n, ser_rate):
        for i in range(n):
            self.cur_t[i] = 0
        hits = 0
        misses = 0
        for i, obj_name in self.q:
            req = Req(obj_name=obj_name)
            req.start(i)
            cur_ser = -1
            cur_time = float('inf')
            for j in range(n):
                if self.cur_t[j] < cur_time:
                    cur_time = self.cur_t[j]
                    cur_ser = j
            #print i, cur_ser, cur_time
            if cur_time < i:
                cur_time = i
            if self.cache.get(obj_name) > 0:
                st = self.get_st(ser_rate)
                hits += 1
            else:
                self.cache.set(obj_name, 1)
                st = self.successor.arrive(i) - cur_time
                misses += 1
            cur_time += st
            self.cur_t[cur_ser] = cur_time
            req.finish(self.cur_t[cur_ser])
            req.ser_t = st
            #print self.cur_t
            self.res.append(req)
        print hits, misses, hits/float(misses+hits)
        return (hits, misses)

    def acproc(self, n, ser_rate):
        for i in range(n):
            self.cur_t[i] = 0
        while self.q:
            ii = self.q.pop(0)
            if self.generator is not None:
                i, obj_name = ii
            else:
                i = ii
            cur_ser = -1
            cur_time = float('inf')
            for j in range(n):
                if self.cur_t[j] < cur_time:
                    cur_time = self.cur_t[j]
                    cur_ser = j
            if cur_time < i:
                cur_time = i
            accept_list = []
            accept_list.append(ii)
            try:
                while self.q[0] <= cur_time:
                    accept_list.append(self.q.pop(0))
            except:
                pass
            #print accept_list
            for t in accept_list:
                if self.generator is not None:
                    t, obj_name = t
                    req = Req(obj_name=obj_name)
                else:
                    req = Req()
                req.start(t)
                st = self.get_st(ser_rate)
                cur_time += st
                self.cur_t[cur_ser] = cur_time
                req.finish(self.cur_t[cur_ser])
                req.ser_t = st
                self.res.append(req)

    def cache_acproc(self, n, ser_rate):
        for i in range(n):
            self.cur_t[i] = 0
            self.ac_list[i] = []
        hits = 0
        misses = 0
        while self.q:
            ii = self.q[0]
            i, obj_name = ii
            cur_ser = -1
            cur_time = float('inf')
            for j in range(n):
                if self.cur_t[j] < cur_time:
                    cur_time = self.cur_t[j]
                    cur_ser = j
            # print ii, cur_ser, cur_time
            # print self.ac_list
            if self.ac_list[cur_ser] != []:
                # print "processing"
                t, obj_name = self.ac_list[cur_ser].pop(0)
                req = Req(obj_name=obj_name)
                req.start(t)
                if self.cache.get(obj_name) > 0:
                    st = self.get_st(ser_rate)
                    hits += 1
                else:
                    self.cache.set(obj_name, 1)
                    st = self.successor.arrive(t) - cur_time
                    misses += 1
                cur_time += st
                req.finish(cur_time)
                req.ser_t = st
                self.res.append(req)
            else:
                # print "accepting"
                if cur_time < i:
                    cur_time = i
                try:
                    while self.q[0][0] <= cur_time:
                        self.ac_list[cur_ser].append(self.q.pop(0))
                except:
                    pass
            self.cur_t[cur_ser] = cur_time
        print hits, misses, hits/float(misses+hits)
        return (hits, misses)

    def mean(self):
        s_time = 0
        n = 0
        ser_time = 0
        for r in self.res:
            s_time += r.duration()
            ser_time += r.ser_t
            n += 1
        print "---------mean--------"
        print s_time / n, ser_time / n

    def percentile(self):
        durs = []
        for r in self.res:
            durs.append(r.duration())
        durs.sort()
        print sum(durs), min(durs), max(durs)
        ress = []
        for i in [5, 10, 20, 40, 60, 80, 90, 95, 99]:
            idx = int(i/100.0*len(durs))
            #print idx
            ress.append((i, durs[idx]))
        print "---------percentile-------"
        print ress
        print "-----end----\n\n\n"
        return ress


if __name__ == '__main__':
    cache = LRUCache(2000)
    gen = zipfian_generator(0, 10000, zipfianconstant=0.9)

    srate = 100000.0
    ratio21 = 1000.0
    cache_hit_ratio = 0.65
    arate_max = srate/ratio21/(1-cache_hit_ratio)
    arate = arate_max*0.7

    successor = simpleQueue(srate/ratio21)

    no = 4

    cache.clear()
    successor.clear()
    q1 = Que(cache=cache, generator=gen, successor=successor)
    q1.warmup(50000)
    q1.start(100000, arate)
    q1.cache_proc(1, srate)
    q1.mean()
    q1.percentile()

    cache.clear()
    successor.clear()
    q2 = Que(cache=cache, generator=gen, successor=successor)
    q2.warmup(50000)
    q2.start(100000, arate)
    hits2, misses2 = q2.cache_proc(no, srate)
    q2.mean()
    ress2 = q2.percentile()

    cache.clear()
    successor.clear()
    q3 = Que(cache=cache, generator=gen, successor=successor)
    q3.warmup(50000)
    q3.start(100000, arate)
    hits3, misses3 = q3.cache_acproc(no, srate)
    q3.mean()
    ress3 = q3.percentile()


    hits = hits2
    misses = misses2
    ress = ress2

    hits = hits3
    misses = misses3
    ress = ress3


    theta = float(misses) / (hits + misses)
    lambd = arate
    mean1 = 1.0/srate
    mean2 = 1.0/(srate/ratio21)

    pdf_laplace_func1 = exp_pdf_laplace
    mean_func1 = exp_mean
    var1 = exp_var(mean1)
    pdf_laplace_func2 = exp_pdf_laplace
    var2 = exp_var(mean2)

    pdf_laplace_func1 = determin_pdf_laplace
    mean_func1 = determin_mean
    var1 = determin_var(mean1)
    pdf_laplace_func2 = determin_pdf_laplace
    var2 = determin_var(mean2)

    cache_pdf_laplace = get_cache_pdf_laplace(lambd, theta, mean1, pdf_laplace_func1, mm1k_mean(mean2, lambd*theta), mm1k_sojourn_pdf_laplace)
    cache_mean_value = cache_mean(lambd, theta, mean1, mean_func1, mean2, mm1k_mean)

    # print mm1k_mean(mean2, lambd*theta), cache_mean_value

    for per, t in ress:
        # print t, mg1_sojourn_cdf(t, cache_mean_value, lambd/10.0, cache_pdf_laplace), cache_cdf(t, lambd, theta, mean1, pdf_laplace_func1, mm1k_mean(mean2, lambd*theta), mm1k_sojourn_pdf_laplace)
        per = per / 100.0
        p_star, p_ag = cosmodel_reqbased_backend(t, lambd, theta, no, mean1, var1, pdf_laplace_func1, mean2, var2, pdf_laplace_func2)
        print "%.8f %.2f (%.4f+%.4f)=%.4f %.4f" % (t, per, p_star, p_ag, (p_star+p_ag), mg1_sojourn_cdf(t, cache_mean_value, lambd/float(no), cache_pdf_laplace))