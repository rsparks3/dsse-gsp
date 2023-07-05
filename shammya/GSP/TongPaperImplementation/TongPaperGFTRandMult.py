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
from keras.layers import Dense, Activation, Dropout
from keras.layers import LSTM
from keras.models import Sequential
from keras.models import load_model
from collections import Counter 
import copy
from keras.callbacks import EarlyStopping,ModelCheckpoint
from sklearn.preprocessing import MinMaxScaler
from scipy.io import loadmat,savemat
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

from random import randint


locations = loadmat('sensor_locations_85bus_excluded.mat')['sensor_locations']
locations = locations.flatten('F') - 1
q_locations = locations + 84
locations = list(locations) + list(q_locations)


mc = ModelCheckpoint('Tong_paper_85bus_included_randmult.h5', monitor='val_loss', mode='min', save_best_only=True)



X = loadmat('input_data_randmult.mat')['input_data_2020']
Y = loadmat('output_data_randmult.mat')['output_data_2020']

X = X[:,locations]
# Y = Y[:,locations]
in_dim = X.shape[1]
out_dim = Y.shape[1]

X_train, X_test, y_train, y_test=train_test_split(X, Y, test_size=0.2)


TRAIN_ROW = 38000

# X_train = demand[0:TRAIN_ROW,locations]
# X_test = demand[TRAIN_ROW:,locations]
# y_train = V[0:TRAIN_ROW,BUS_CHOSEN]
# y_test = V[TRAIN_ROW:,BUS_CHOSEN]


model = Sequential()
model.add(Dense(1000, input_dim=in_dim, activation='tanh'))
model.add(Dense(750, activation='tanh'))
model.add(Dense(500, activation='tanh'))
model.add(Dense(250, activation='tanh'))
model.add(Dense(125, activation='tanh'))
model.add(Dense(50, activation='tanh'))
model.add(Dense(out_dim, activation='linear'))

model.compile(loss="mse", optimizer="adam")
es = EarlyStopping(monitor='val_loss', min_delta=0, patience=10, verbose=1, mode='min')

# train the model
model.fit(X_train, y_train, batch_size=500, epochs=1000, validation_split=0.25, verbose=1,shuffle=False,callbacks=[es,mc])


ypred = model.predict(X_test)
print("y1 MSE:%.4f" % mean_squared_error(y_test[-1,:], ypred[-1,:])) 

print("y2 MSE:%.4f" % mean_squared_error(y_test[0,:], ypred[0,:]))

 
test_mse = model.evaluate(X_test,y_test)
print("MSE:%.4f" % test_mse) 


mse_list = []
for i in range(0,y_test.shape[0]):
    mse_list.append(mean_squared_error(y_test[i,:], ypred[i,:]))

mse_list = np.asarray(mse_list)


# plt.figure(1)
# x_ax = range(y_test.shape[1])
# plt.plot(x_ax, y_test[-1,:], label="y2-test")
# plt.plot(x_ax, ypred[-1,:], label="y2-pred")
# plt.legend()
# plt.show()




plt.figure(1,dpi=300)
x_ax = range(y_test.shape[1])
plt.subplot(2, 1, 1)
bus_numbers = range(1,86)
plt.plot(bus_numbers, y_test[-1,0:85], label="Original Voltage Magnitude")
plt.plot(bus_numbers, ypred[-1,0:85], label="Predicted Voltage Magnitude")
# plt.xlabel('Bus Numbers')
plt.ylabel('Voltage Magnitude (pu)')
plt.legend()

locs, labels = plt.xticks()
plt.xticks(np.arange(1, 86, step=5))
plt.subplot(2, 1, 2)
plt.plot(bus_numbers, y_test[-1,85:], label="Original Voltage Angle")
plt.plot(bus_numbers, ypred[-1,85:], label="Predicted Voltage Angle")
locs, labels = plt.xticks()
plt.xticks(np.arange(1, 86, step=5))
plt.xlabel('Bus Numbers')
plt.ylabel('Voltage Angle (Degree)')
plt.legend()
# plt.title('Min Error Case')
plt.show()



min_mse_location = np.argmin(mse_list)
plt.figure(2)
x_ax = range(ypred.shape[1])
plt.plot(x_ax, y_test[min_mse_location,:], label="ytest")
plt.plot(x_ax, ypred[min_mse_location,:], label="ypred")
plt.legend()
plt.title('Min Error Case')
plt.show()

max_mse_location = np.argmax(mse_list)
plt.figure(3)
x_ax = range(ypred.shape[1])
plt.plot(x_ax, y_test[max_mse_location,:], label="ytest")
plt.plot(x_ax, ypred[max_mse_location,:], label="ypred")
plt.legend()
plt.title('Max Error Case')
plt.show()

savemat('ypred.mat',{'loc':ypred[-1,:],'label':'ypred'})
savemat('ytest.mat',{'loc':y_test[-1,:],'label':'ytest'})


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
