# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 20:42:00 2021

@author: asus
"""

from PIL import Image
import numpy as np

ipath = "./plots/features/feature_transform_%d.png"

import imageio
images = []
for i in range(11):
    images.append(imageio.imread(ipath%(i*1000)))
imageio.mimsave('./images/feature_transform_57.gif', images, fps=2)


ipath = "./plots/features/feature_transform_centroid_%d.png"

images = []
for i in range(11):
    images.append(imageio.imread(ipath%(i*1000)))
imageio.mimsave('./images/feature_transform_centroid.gif', images, fps=2)

ipath = "./plots/iterations/iteration_%d.png"

images = []
for i in range(11):
    images.append(imageio.imread(ipath%(i*1000)))
imageio.mimsave('./images/iterations_57.gif', images, fps=2)
