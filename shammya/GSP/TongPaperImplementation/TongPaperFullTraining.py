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
# import tensorflow.keras as keras
from tensorflow.keras.layers import Dense, Activation, Dropout
from tensorflow.keras.layers import Embedding
from tensorflow.keras.layers import LSTM
from tensorflow.keras.models import Sequential
from tensorflow.keras.models import load_model
from collections import Counter 
import copy
from tensorflow.keras.callbacks import EarlyStopping,ModelCheckpoint
from sklearn.preprocessing import MinMaxScaler
from scipy.io import loadmat,savemat
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error



from random import randint




parser = argparse.ArgumentParser()
    # parse the configuration file, otherwise throw an error  
parser.add_argument('--i', type=int, default = 20)
args = parser.parse_args()


INFORMATION_AVAILABLE = args.i
x = np.array([randint(0, 85) for p in range(0, INFORMATION_AVAILABLE)])
current_measurement = np.array([randint(171, 184) for p in range(0, 6)])
locations = list(x) + list(x+84) + list(current_measurement)

# savemat('locations_{}.mat'.format(INFORMATION_AVAILABLE+6),{'loc':(locations),'label':str(INFORMATION_AVAILABLE+6)})

mc = ModelCheckpoint('Tong_paper.h5', monitor='val_loss', mode='min', save_best_only=True)



X = loadmat('input_data.mat')['input_data']
Y = loadmat('output_data.mat')['output_data']

# X = X[:,locations]
# Y = Y[:,locations]
in_dim = X.shape[1]
out_dim = Y.shape[1]

X_train, X_test, y_train, y_test=train_test_split(X, Y, test_size=0.2)


TRAIN_ROW = 49500

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
model.add(Dense(out_dim, activation='tanh'))

model.compile(loss="mse", optimizer="adam")
es = EarlyStopping(monitor='val_loss',
                              min_delta=0,
                              patience=10,
                              verbose=1, mode='min')

# train the model
model.fit(X_train, y_train, batch_size=500, epochs=300, validation_split=0.25, verbose=1,shuffle=False,callbacks=[es,mc])


ypred = model.predict(X_test)
print("y1 MSE:%.4f" % mean_squared_error(y_test[:,0], ypred[:,0])) 

print("y2 MSE:%.4f" % mean_squared_error(y_test[:,1], ypred[:,1]))

 
test_mse = model.evaluate(X_test,y_test)
print("MSE:%.4f" % test_mse) 



x_ax = range(len(X_test))
plt.plot(x_ax, y_test[:,1], label="y2-test")
plt.plot(x_ax, ypred[:,1], label="y2-pred")
plt.legend()
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
