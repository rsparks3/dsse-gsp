# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 03:18:17 2020

@author: sssaha
"""


from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from keras.layers.core import Dense, Activation, Dropout
from keras.layers import Embedding
from keras.layers.recurrent import LSTM
from keras.models import Sequential
from keras.models import load_model
from collections import Counter 
import copy
from keras.callbacks import EarlyStopping,ModelCheckpoint
from sklearn.preprocessing import MinMaxScaler
from scipy.io import loadmat
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from statistics import mean
 



from random import randint


parser = argparse.ArgumentParser()
    # parse the configuration file, otherwise throw an error  
parser.add_argument('--i', type=int, default = 20)
args = parser.parse_args()


INFORMATION_AVAILABLE = args.i

locations = loadmat('locations_{}.mat'.format(INFORMATION_AVAILABLE+6))['loc']

X_test = loadmat('input_data_2020.mat')['input_data_2020']
y_test = loadmat('output_data_2020.mat')['output_data_2020']

X_test = X_test[:,list(locations[0])]

# in_dim = X.shape[1]
# out_dim = Y.shape[1]

# X_train, X_test, y_train, y_test=train_test_split(X, Y, test_size=0.2)


# TRAIN_ROW = 49500

# X_train = demand[0:TRAIN_ROW,locations]
# X_test = demand[TRAIN_ROW:,locations]
# y_train = V[0:TRAIN_ROW,BUS_CHOSEN]
# y_test = V[TRAIN_ROW:,BUS_CHOSEN]


# model = Sequential()
# model.add(Dense(1000, input_dim=in_dim, activation='tanh'))
# # model.add(Dense(750, activation='tanh'))
# model.add(Dense(500, activation='tanh'))
# model.add(Dense(250, activation='tanh'))
# model.add(Dense(125, activation='tanh'))
# model.add(Dense(50, activation='tanh'))
# model.add(Dense(out_dim, activation='tanh'))

# model.compile(loss="mse", optimizer="adam")
# es = EarlyStopping(monitor='val_loss',
#                               min_delta=0,
#                               patience=10,
#                               verbose=1, mode='min')

# train the model
# model.fit(X_train, y_train, batch_size=500, epochs=300, validation_split=0.25, verbose=1,shuffle=False,callbacks=[es,mc])

model = load_model('Tong_paper_{}.h5'.format(INFORMATION_AVAILABLE+6))
# model = load_model('Tong_paper.h5')
ypred = model.predict(X_test)
mse_list=[]

for i in range(0,y_test.shape[0]):
    mse_list.append(mean_squared_error(y_test[i,:], ypred[i,:]))
    # print("y1 MSE:%.4f" % mean_squared_error(y_test[i,:], ypred[i,:])) 
    
    # print("y2 MSE:%.4f" % mean_squared_error(y_test[:,1], ypred[:,1]))

print(max(mse_list)) 


# ASE Calculation
ase = []

for i in range(96):
    error_data = []
    for j in range(0,7393,96):
        error_data.append( mean_squared_error(y_test[i+j,:],ypred[i+j,:]))
    ase.append(mean(error_data))


min_mse_location = np.argmin(mse_list)
plt.figure(1,dpi=300)
x_ax = range(ypred.shape[1])
plt.subplot(2, 1, 1)
bus_numbers = range(1,86)
plt.plot(bus_numbers, y_test[min_mse_location,0:85], label="Original Voltage Magnitude")
plt.plot(bus_numbers, ypred[min_mse_location,0:85], label="Predicted Voltage Magnitude")
# plt.xlabel('Bus Numbers')
plt.ylabel('Voltage Magnitude (pu)')
plt.legend()

locs, labels = plt.xticks()
plt.xticks(np.arange(1, 86, step=5))
plt.subplot(2, 1, 2)
plt.plot(bus_numbers, y_test[min_mse_location,85:], label="Original Voltage Angle")
plt.plot(bus_numbers, ypred[min_mse_location,85:], label="Predicted Voltage Angle")
locs, labels = plt.xticks()
plt.xticks(np.arange(1, 86, step=5))
plt.xlabel('Bus Numbers')
plt.ylabel('Voltage Angle (Degree)')
plt.legend()
# plt.title('Min Error Case')
plt.show()
# plt.savefig('MinErrorCase.png')



# max_mse_location = np.argmax(mse_list)
# plt.figure(2)
# x_ax = range(ypred.shape[1])
# plt.plot(x_ax, y_test[max_mse_location,:], label="Original_"+str(args.i))
# plt.plot(x_ax, ypred[max_mse_location,:], label="Predicted_"+str(args.i))
# plt.legend()
# plt.title('Max Error Case')
# plt.show()


plt.figure(3,dpi=300)
plt.plot(ase)
plt.title('Hourly Average Squarred Error')
plt.show()


# test_mse = model.evaluate(X_test, y_test, verbose=1)
# print(f'Test MSE:{test_mse}')

# predicted_values = model.predict(X_test)

# num_test_samples = len(predicted_values)
# predicted_values = np.reshape(predicted_values, (num_test_samples,1))

# fig = plt.figure()
# plt.plot(y_test)
# plt.plot(predicted_values)
# plt.plot(y_test + shifted_value)
# plt.plot(predicted_values + shifted_value)
# plt.legend(['Main Value','predicted value'])
# plt.xlabel('Scenario')
# plt.ylabel('Bus Voltage of {}'.format(BUS_CHOSEN))
# plt.title('With Information from {} Buses'.format(INFORMATION_AVAILABLE))
# plt.show()
