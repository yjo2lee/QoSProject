# -*- coding: utf-8 -*-

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from numpy.linalg import inv
from scipy.stats import norm
from keras.models import Sequential
from keras.layers import Dense, Activation
import keras
import tensorflow as tf
from keras.utils.generic_utils import get_custom_objects
from scipy.optimize import brentq
from matplotlib import font_manager, rc
import matplotlib

def TPtrans(raw, a, b, c):
    transTP = a*np.log(raw*1000) + b*raw*1000 + c
    return transTP

def BitTPtrans(raw, a, b, c):
    transTP = a*np.log(raw) + b*raw + c
    return transTP

def RBtrans(raw, a, b, c):
    transRB = a*np.log(raw) - b*np.log(1-raw) + c
    return transRB

def normalized_TP(transTP, cov, transRB):
    return (transTP - cov*transRB)/np.sqrt(1-cov**2)

# data: qos
def ProbQosData(Param):
    data_len = len(Param)
    rb_array = np.random.uniform(0.4, 1-1e-6, data_len)
    total_dat = np.zeros((data_len, 5))
    total_dat[:,0] = np.random.uniform(0.5, 4.0, data_len)
    for i in range(len(Param)):
        transRB = RBtrans(rb_array[i], Param[i,1], Param[i,2], Param[i,3])
        transTP = TPtrans(total_dat[i,0]*Param[i,0], Param[i,4], Param[i, 5], Param[i,6])
        x = normalized_TP(transTP, Param[i,7], transRB)
        if not np.isnan(x): 
            total_dat[i, 3] = stats.norm.cdf(x, loc = 0, scale = 1)
            total_dat[i, 1] = Param[i, 0]
            total_dat[i, 2] = rb_array[i]
            total_dat[i, 4] = x
    return total_dat

# data: spectrum efficiency

def gderivative(a, b, x):
    return(a/x + b/(1-x))

def objftn(x, params, u):
    return(params[1]*np.log(x)-params[2]*np.log(1-x) + params[3] - u)

def rbpercent(params):
    u = np.random.normal(0,1,1)[0]
    if params[1] == 0:
        while u <= params[3]:
            u = np.random.normal(0,1,1)[0]
    elif (params[1] > 0) and (params[2] == 0):
        while u >= params[3]:
            u = np.random.normal(0,1,1)[0]
    else:
        None
    return u

def SEData(Param):
    data_len = len(Param)
    total_dat = np.zeros((data_len, 6))
    for i in range(len(Param)):
        transrb = rbpercent(Param[i])   #u
        alpha = np.random.uniform(0, 10.0)   # alpha
        beta = alpha*1e7/1e3
        #beta = 2*alpha*1e7/1e3   #20MHz
        r = brentq(lambda x: Param[i,1]*np.log(x)-Param[i,2]*np.log(1-x) + Param[i,3] - transrb, 1e-9, 1-1e-9)
        arg = BitTPtrans(beta*r, Param[i,4], Param[i, 5], Param[i,6])
        x = normalized_TP(arg, Param[i,7], transrb)
        total_dat[i, 0] = Param[i, 0]           # rnti
        total_dat[i, 1] = alpha                 # alpha
        total_dat[i, 2] = transrb               # u
        total_dat[i, 3] = x                     # output
        total_dat[i, 4] = 1-stats.norm.cdf(x)   # prob
        total_dat[i, 5] = r
    return total_dat

def ConstructModel(num_layer, weight, batch_size):
    keras.backend.clear_session()
    model = Sequential()
    inputlayer = keras.layers.Input(shape = (3,), name = 'input')
    dense = keras.layers.Dense(weight, activation='relu')(inputlayer)
    for k in range(num_layer-1):
        dense = keras.layers.Dense(weight, activation='relu')(dense)
    outputlayer = keras.layers.Dense(1, name = 'output')(dense)
    model = keras.models.Model(inputs = inputlayer, outputs = outputlayer)
    model.compile(loss='mean_squared_error', optimizer='adam')
    model.summary()
    return model

def hbar(alpha, n, rep, model):
    testinput = np.hstack((n, alpha, np.random.normal(0,1,rep).reshape(-1,1)))
    return np.mean(1-stats.norm.cdf(model.predict(testinput)))


def train():
    # load data
    raw_data_10= np.load('Raw_data_10MHz.npy')
    Param = raw_data_10
    total_dat = ProbQosData(Param)
    total_SE = SEData(Param)

    samplesize = len(total_SE)
    train_size = int(samplesize*0.6)
    train_set, valid_set = total_SE[0:train_size,:], total_SE[train_size:samplesize,:]
    patience = 10
    weight = 20
    batch_size = 30
    num_layer = 5

    train_model = ConstructModel(num_layer, weight, batch_size)
    callbacks = [keras.callbacks.EarlyStopping(monitor='val_loss', patience=patience),
                    keras.callbacks.ModelCheckpoint(filepath='se_model_%d_%d_%d_ver2_bw10.h5' % (num_layer, weight, batch_size), monitor='val_loss', save_best_only=True)]
    # callbacks = [keras.callbacks.EarlyStopping(monitor='val_loss', patience=patience),
    #                 keras.callbacks.ModelCheckpoint(filepath='se_model_%d_%d_%d_ver2_bw20.h5' % (num_layer, weight, batch_size), monitor='val_loss', save_best_only=True)]
    train_model.fit(train_set[:,0:3],train_set[:,3],epochs = 1000, callbacks = callbacks, batch_size = batch_size, validation_data=(valid_set[:,0:3], valid_set[:,3]), verbose = 1)

def test():
    resultlen = 22
    alpha_plot = np.linspace(0.0, 7.0, 29)
    test_model = keras.models.load_model('se_model_%d_%d_%d_ver2_bw10.h5' % (num_layer, weight, batch_size))
    #model = keras.models.load_model('se_model_%d_%d_%d_ver2_bw20.h5' % (num_layer, weight, batch_size))
    rnti = np.linspace(3,resultlen+2,resultlen)
    plot_se_data = np.zeros((resultlen, 29))
    rep = 500; i = 0

    for i in range(resultlen):
        n = np.array([rnti[i]]*rep).reshape(-1,1)
        j = 0
        for j in range(len(alpha_plot)):
            tempalpha = np.array([alpha_plot[j]]*rep).reshape(-1,1)
            plot_se_data[i,j] = hbar(tempalpha, n, rep, test_model)

def main():
    train()
    test()

if __name__ == "__main__":
    main()
    
