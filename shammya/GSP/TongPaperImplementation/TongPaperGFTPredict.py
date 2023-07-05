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
from tensorflow.keras.models import load_model

from tensorflow.keras.initializers import glorot_uniform

# loaded_model = tf.keras.models.load_model("pruned.h5",custom_objects={'GlorotUniform': glorot_uniform()})


from scipy.io import loadmat
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from statistics import mean
from scipy.io import loadmat,savemat




locations = loadmat('sensor_locations_85bus_included.mat')['sensor_locations']
locations = locations.flatten('F') - 1 
q_locations = locations + 85
locations = list(locations) + list(q_locations)

X_test = loadmat('input_data_2020.mat')['input_data_2020']
Y_test = loadmat('output_data_2020.mat')['output_data_2020']



x_test = loadmat('base_demand.mat')['demand'][0]
y_test = loadmat('base_state.mat')['state'][0]

for i in range(0,Y_test.shape[1]):
    Y_test[0,i] = y_test[i]

for i in range(0,x_test.shape[0]):
    X_test[0,i] = x_test[i]

X_test = X_test[:,locations]


parser = argparse.ArgumentParser()
    # parse the configuration file, otherwise throw an error  
parser.add_argument('--i', type=int, default = 20)
args = parser.parse_args()


INFORMATION_AVAILABLE = args.i

locations = loadmat('locations_{}.mat'.format(INFORMATION_AVAILABLE+6))['loc']




# X_test = np.reshape(X_test, (1,X_test.shape[0]))

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

model = load_model('Tong_paper_85bus_included.h5',custom_objects={'GlorotUniform': glorot_uniform()})
# model = load_model('Tong_paper.h5')
ypred = model.predict(X_test)
mse_list=[]

print(Y_test[0,:])

plt.figure(1)
x_ax = range(y_test.shape[0])
plt.plot(x_ax, y_test, label="Original_states")
plt.plot(x_ax, ypred[0,:], label="Predicted_states")
plt.legend()
plt.title('State Error Case')
plt.show()


savemat('ytest.mat',{'ytest':y_test,'label':'ytest'})
savemat('ypred.mat',{'ypred':ypred[0,:],'label':'ypred'})

# for i in range(0,y_test.shape[0]):
#     mse_list.append(mean_squared_error(y_test[i,:], ypred[i,:]))
#     # print("y1 MSE:%.4f" % mean_squared_error(y_test[i,:], ypred[i,:])) 
    
#     # print("y2 MSE:%.4f" % mean_squared_error(y_test[:,1], ypred[:,1]))

# print(max(mse_list)) 


# # ASE Calculation
# ase = []

# for i in range(96):
#     error_data = []
#     for j in range(0,7393,96):
#         error_data.append( mean_squared_error(y_test[i+j,:],ypred[i+j,:]))
#     ase.append(mean(error_data))


# min_mse_location = np.argmin(mse_list)
# plt.figure(1)
# x_ax = range(ypred.shape[1])
# plt.plot(x_ax, y_test[min_mse_location,:], label="Original_"+str(args.i))
# plt.plot(x_ax, ypred[min_mse_location,:], label="Predicted_"+str(args.i))
# plt.legend()
# plt.title('Min Error Case')
# plt.show()

# max_mse_location = np.argmax(mse_list)
# plt.figure(2)
# x_ax = range(ypred.shape[1])
# plt.plot(x_ax, y_test[max_mse_location,:], label="Original_"+str(args.i))
# plt.plot(x_ax, ypred[max_mse_location,:], label="Predicted_"+str(args.i))
# plt.legend()
# plt.title('Max Error Case')
# plt.show()


# plt.figure(3)
# plt.plot(ase,label= 'ASE_'+str(args.i))
# plt.legend()
# plt.title('ASE')
# plt.show()


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
