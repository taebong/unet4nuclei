
# coding: utf-8

# # Step 02
# # Training a U-Net model

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'inline')


# In[5]:


import sys
import os

import numpy as np
import skimage.io

import tensorflow as tf

import keras.backend
import keras.callbacks
import keras.layers
import keras.models
import keras.optimizers

import utils.model_builder
import utils.data_provider
import utils.metrics
import utils.objectives
import utils.dirtools

# Uncomment the following line if you don't have a GPU
#os.environ['CUDA_VISIBLE_DEVICES'] = ''


# # Configuration

# In[6]:


from config import config_vars


# In[7]:


#experiment_name = '03'   #03: more training data (75% of all), do augmentation
experiment_name = '01_MEF'     # training data 75%, val 25%, augmentation

config_vars = utils.dirtools.setup_experiment(config_vars, experiment_name)

data_partitions = utils.dirtools.read_data_partitions(config_vars)

config_vars


# # Initiate data generators

# In[8]:


# build session running on GPU 1
configuration = tf.ConfigProto()
configuration.gpu_options.allow_growth = True
configuration.gpu_options.visible_device_list = "0"
session = tf.Session(config = configuration)

# apply session
keras.backend.set_session(session)

train_gen = utils.data_provider.random_sample_generator(
    config_vars["normalized_images_dir"],
    config_vars["boundary_labels_dir"],
    data_partitions["training"],
    config_vars["batch_size"],
    config_vars["pixel_depth"],
    config_vars["crop_size"],
    config_vars["crop_size"],
    config_vars["rescale_labels"]
)

val_gen = utils.data_provider.single_data_from_images(
     config_vars["normalized_images_dir"],
     config_vars["boundary_labels_dir"],
     data_partitions["validation"],
     config_vars["val_batch_size"],
     config_vars["pixel_depth"],
     config_vars["crop_size"],
     config_vars["crop_size"],
     config_vars["rescale_labels"]
)


# # Build model

# In[9]:


# build model
model = utils.model_builder.get_model_3_class(config_vars["crop_size"], config_vars["crop_size"], activation=None)
model.summary()

#loss = "categorical_crossentropy"
loss = utils.objectives.weighted_crossentropy

metrics = [keras.metrics.categorical_accuracy, 
           utils.metrics.channel_recall(channel=0, name="background_recall"), 
           utils.metrics.channel_precision(channel=0, name="background_precision"),
           utils.metrics.channel_recall(channel=1, name="interior_recall"), 
           utils.metrics.channel_precision(channel=1, name="interior_precision"),
           utils.metrics.channel_recall(channel=2, name="boundary_recall"), 
           utils.metrics.channel_precision(channel=2, name="boundary_precision"),
          ]

optimizer = keras.optimizers.RMSprop(lr=config_vars["learning_rate"])

model.compile(loss=loss, metrics=metrics, optimizer=optimizer)

# Performance logging
callback_csv = keras.callbacks.CSVLogger(filename=config_vars["csv_log_file"])

callbacks=[callback_csv]


# # Training 

# In[10]:


# TRAIN
statistics = model.fit_generator(
    generator=train_gen,
    steps_per_epoch=config_vars["steps_per_epoch"],
    epochs=config_vars["epochs"],
    validation_data=val_gen,
    validation_steps=int(len(data_partitions["validation"])/config_vars["val_batch_size"]),
    callbacks=callbacks,
    verbose = 1
)

model.save_weights(config_vars["model_file"])

print('Done! :)')

