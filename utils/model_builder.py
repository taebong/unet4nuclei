import keras.layers
import keras.models
import tensorflow as tf

CONST_DO_RATE = 0.5

option_dict_conv = {"activation": "relu", "padding": "same"}
option_dict_bn = {"mode": 0, "momentum" : 0.9}


# returns a core model from gray input to 64 channels of the same size
def get_core(dim1, dim2):
    
    x = keras.layers.Input(shape=(dim1, dim2, 1))

    a = keras.layers.Conv2D(64, (3,3), **option_dict_conv)(x)  
    a = keras.layers.BatchNormalization(**option_dict_bn)(a)

    a = keras.layers.Conv2D(64, (3,3), **option_dict_conv)(a)
    a = keras.layers.BatchNormalization(**option_dict_bn)(a)

    
    y = keras.layers.MaxPooling2D()(a)

    b = keras.layers.Conv2D(128, (3,3), **option_dict_conv)(y)
    b = keras.layers.BatchNormalization(**option_dict_bn)(b)

    b = keras.layers.Conv2D(128, (3,3), **option_dict_conv)(b)
    b = keras.layers.BatchNormalization(**option_dict_bn)(b)

    
    y = keras.layers.MaxPooling2D()(b)

    c = keras.layers.Conv2D(256, (3,3), **option_dict_conv)(y)
    c = keras.layers.BatchNormalization(**option_dict_bn)(c)

    c = keras.layers.Conv2D(256, (3,3), **option_dict_conv)(c)
    c = keras.layers.BatchNormalization(**option_dict_bn)(c)

    
    y = keras.layers.MaxPooling2D()(c)

    d = keras.layers.Conv2D(512, (3,3), **option_dict_conv)(y)
    d = keras.layers.BatchNormalization(**option_dict_bn)(d)

    d = keras.layers.Conv2D(512, (3,3), **option_dict_conv)(d)
    d = keras.layers.BatchNormalization(**option_dict_bn)(d)

    
    # UP

    d = keras.layers.UpSampling2D()(d)

    #crop to match size
    ch, cw = get_crop_shape(c, d)
    crop_c = keras.layers.Cropping2D(cropping=(ch,cw))(c)    

    y = keras.layers.merge.concatenate([d, crop_c], axis=3)

    e = keras.layers.Conv2D(256, (3,3), **option_dict_conv)(y)
    e = keras.layers.BatchNormalization(**option_dict_bn)(e)

    e = keras.layers.Conv2D(256, (3,3), **option_dict_conv)(e)
    e = keras.layers.BatchNormalization(**option_dict_bn)(e)

    e = keras.layers.UpSampling2D()(e)

    ch, cw = get_crop_shape(b, e)           
    crop_b = keras.layers.Cropping2D(cropping=(ch,cw))(b)
    
    y = keras.layers.merge.concatenate([e, crop_b], axis=3)

    f = keras.layers.Conv2D(128, (3,3), **option_dict_conv)(y)
    f = keras.layers.BatchNormalization(**option_dict_bn)(f)

    f = keras.layers.Conv2D(128, (3,3), **option_dict_conv)(f)
    f = keras.layers.BatchNormalization(**option_dict_bn)(f)

    f = keras.layers.UpSampling2D()(f)

    ch, cw = get_crop_shape(a, f)
    crop_a = keras.layers.Cropping2D(cropping=(ch,cw))(a)

    y = keras.layers.merge.concatenate([f, crop_a], axis=3)

    y = keras.layers.Conv2D(64, (3,3), **option_dict_conv)(y)
    y = keras.layers.BatchNormalization(**option_dict_bn)(y)

    y = keras.layers.Conv2D(64, (3,3), **option_dict_conv)(y)
    y = keras.layers.BatchNormalization(**option_dict_bn)(y)

    ch, cw = get_crop_shape(x,y)
    y = keras.layers.ZeroPadding2D(padding=((ch[0],ch[1]), (cw[0],cw[1])))(y)
    y = keras.layers.Conv2D(3,(1,1))(y)

    return [x, y]


def get_model_3_class(dim1, dim2, activation="softmax"):
    
    [x, y] = get_core(dim1, dim2)

    y = keras.layers.Conv2D(3, (1,1), **option_dict_conv)(y)

    if activation is not None:
        y = keras.layers.Activation(activation)(y)

    model = keras.models.Model(x, y)
    
    return model


def get_crop_shape(target, refer):
    # width, the 3nd dimension
    cw = (target.get_shape()[2] - refer.get_shape()[2]).value
    assert (cw >= 0)
    if cw % 2 != 0:
        cw1, cw2 = int(cw/2), int(cw/2) + 1
    else:
        cw1, cw2 = int(cw/2), int(cw/2)
    # height, the 2st dimension
    ch = (target.get_shape()[1] - refer.get_shape()[1]).value
    assert (ch >= 0)
    if ch % 2 != 0:
        ch1, ch2 = int(ch/2), int(ch/2) + 1
    else:
        ch1, ch2 = int(ch/2), int(ch/2)

    return (ch1, ch2), (cw1, cw2)
