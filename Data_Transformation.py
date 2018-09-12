# -*- coding: utf-8 -*-

import openpyxl
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def QQ(vec):
    v = np.sort(vec)
    q = stats.norm.ppf((np.arange(len(v))+1-0.375)/(len(v)+0.25))
    return v, q

def RB_fit(RB):
    RB_ord = np.sort(RB)/100
    th_q = stats.norm.ppf((np.arange(len(RB))+1-0.375)/(len(RB)+0.25))
    q = th_q
    l = np.log(RB_ord)
    h = np.log(1-RB_ord)
    d = (np.mean(l*l) - pow(np.mean(l),2))*(np.mean(h*h) - pow(np.mean(h),2)) - pow(np.mean(l*h) - np.mean(l)*np.mean(h),2)
    d1 = (np.mean(h*h) - pow(np.mean(h),2))*(np.mean(l*q) - np.mean(l)*np.mean(q))-(np.mean(l*h) - np.mean(l)*np.mean(h))*(np.mean(h*q) - np.mean(h)*np.mean(q))
    d2 = (np.mean(l*h) - np.mean(l)*np.mean(h))*(np.mean(l*q) - np.mean(l)*np.mean(q)) - (np.mean(l*l) - pow(np.mean(l),2))*(np.mean(h*q) - np.mean(h)*np.mean(q))
    if d1 < 0:
        a = 0
        b = -(np.mean(h*q) - np.mean(h)*np.mean(q))/(np.mean(h*h) - pow(np.mean(h),2))
    elif d2 < 0:
        a = (np.mean(l*q) - np.mean(l)*np.mean(q))/(np.mean(l*l) - pow(np.mean(l),2))
        b = 0
    else:
        a = d1/d
        b = d2/d
    c = np.mean(q) - a*np.mean(l) + b*np.mean(h)
    RB_t = a*l - b*h + c
    _, pv = stats.kstest(RB_t,'norm')
    return a, b, c, pv, RB_t, th_q

def TP_fit(TP):
    TP_ord = np.sort(TP)
    th_q = stats.norm.ppf((np.arange(len(TP))+1-0.375)/(len(TP)+0.25))
    q = th_q
    g = np.log(TP_ord)
    t = TP_ord
    d = (np.mean(g*g) - pow(np.mean(g),2))*(np.mean(t*t) - pow(np.mean(t),2)) - pow(np.mean(g*t) - np.mean(g)*np.mean(t),2)
    d1 = (np.mean(t*t) - pow(np.mean(t),2))*(np.mean(g*q) - np.mean(g)*np.mean(q))-(np.mean(g*t) - np.mean(g)*np.mean(t))*(np.mean(t*q) - np.mean(t)*np.mean(q))
    d2 = (np.mean(g*g) - pow(np.mean(g),2))*(np.mean(t*q) - np.mean(t)*np.mean(q))-(np.mean(t*g) - np.mean(t)*np.mean(g))*(np.mean(g*q) - np.mean(g)*np.mean(q))
    if d1 < 0:
        a = 0
        b = (np.mean(g*q) - np.mean(g)*np.mean(q))/(np.mean(g*g) - pow(np.mean(g),2))
    elif d2 < 0:
        a = (np.mean(t*q) - np.mean(t)*np.mean(q))/(np.mean(t*t) - pow(np.mean(t),2))
        b = 0
    else:
        a = d1/d
        b = d2/d
    c = np.mean(q) - a*np.mean(g) - b*np.mean(t)
    TP_t = a*g + b*t + c
    _, pv = stats.kstest(TP_t,'norm')
    return a, b, c, pv, TP_t, th_q

RD = openpyxl.load_workbook('data_prep_ver2.xlsx', data_only=True)
st = RD['Sheet1']

Data = np.zeros((1,5))
i = 2
value = st['B%d' %i].value
while value != None:
    Data = np.vstack((Data, np.array([st['B%d' %i].value, st['C%d' % i].value, st['D%d' %i].value, st['E%d' %i].value, st['FB%d' %i].value])))
    # PCI, RB, RNTI, TP, Timestamp
    i += 1
    value = st['B%d' %i].value
Data = Data[1:,:].astype(np.float64)

RB_test = Data[Data[:,2] == 4,1]
sample_size = 100
RB_rnd = RB_test[np.random.choice(len(RB_test),sample_size)]
a, b, c, pv, rt, tq = RB_fit(RB_rnd)
print('a=', a, 'b=', b, 'c=', c, 'pv=', pv)
plt.figure(figsize=(10,10))
f = plt.subplot(1, 1, 1)
f.plot(tq, rt, color = 'tab:blue', linewidth = 0, marker='o', markersize = 6)

TP_test = Data[Data[:,2] == 4,3]
sample_size = 100
TP_rnd = TP_test[np.random.choice(len(TP_test),sample_size)]
a, b, c, pv, tt, tq = TP_fit(TP_rnd)
print('a=', a, 'b=', b, 'c=', c, 'pv=', pv)
plt.figure(figsize=(10,10))
f = plt.subplot(1, 1, 1)
f.plot(tq, tt, color = 'tab:blue', linewidth = 0, marker='o', markersize = 6)