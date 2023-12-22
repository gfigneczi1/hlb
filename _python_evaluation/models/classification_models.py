import tensorflow as tf
from tensorflow.keras.optimizers import schedules, SGD, Adam
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import *
from tensorflow import keras
from tensorflow.keras.regularizers import l2
'''
very basic CNN model
'''
def simpleCNN():
    channel = 32
    model = Sequential()
    #model.add(Conv2D(filters=channel, kernel_size=3, activation='relu', input_shape=(120, 14, 1)))
    model.add(Conv2D(filters=channel, kernel_size=5, activation='relu', padding='valid', kernel_initializer='he_normal', input_shape=(42, 40, 1)))
    model.add(BatchNormalization())
    #model.add(MaxPooling2D(2))
    #model.add(Dropout(0.2))
    model.add(Conv2D(filters=channel*2, kernel_size=3, activation='relu', padding='valid', kernel_initializer='he_normal'))
    model.add(BatchNormalization())
    #model.add(MaxPooling2D(2))
    #model.add(Dropout(0.2))
    model.add(Conv2D(filters=channel*2, kernel_size=3, activation='relu', padding='valid', kernel_initializer='he_normal'))
    model.add(BatchNormalization())
    #model.add(Conv2D(filters=channel*3, kernel_size=3, activation='relu'))
    #model.add(BatchNormalization())
    #model.add(MaxPooling2D(2))
    #model.add(Dropout(0.25))
    model.add(Flatten())
    model.add(Dense(channel*2))
    model.add(Dropout(0.4))
    model.add(Activation('relu'))
    model.add(Dense(1, activation='sigmoid'))
    return model
'''
very basic CNN model that uses sequential Zero Padding
'''
def sam_cnn():
    model = Sequential()
    model.add(ZeroPadding2D(((1,1),(2,0))))
    model.add(Conv2D(16, 3, activation='sigmoid'))
    model.add(BatchNormalization())
    model.add(MaxPooling2D(2))
    model.add(Dropout(0.2))
    model.add(ZeroPadding2D(((1,1),(2,0))))
    model.add(Conv2D(32, 3, activation='sigmoid'))
    model.add(BatchNormalization())
    model.add(MaxPooling2D(2))
    model.add(Dropout(0.2))
    model.add(Flatten())
    model.add(Dense(1, activation='sigmoid'))
    return model


def LeNet5(input_shape = (36, 35, 1), classes = 2):
    """
    Implementation of a modified LeNet-5.
    Modified Architecture -- ConvNet --> Pool --> ConvNet --> Pool --> (Flatten) --> FullyConnected --> FullyConnected --> Softmax 

    Arguments:
    input_shape -- shape of the images of the dataset
    classes -- integer, number of classes

    Returns:
    model -- a Model() instance in Keras
    """
    
    model = Sequential([
        
    # Layer 1
    Conv2D(filters = 16, kernel_size = 5, strides = 1, activation = 'relu', input_shape = input_shape, kernel_regularizer=l2(0.0005), name = 'convolution_1'),
    MaxPooling2D(pool_size = 2, strides = 2, name = 'max_pool_1'),
    Dropout(0.25, name = 'dropout_1'),
        
    # Layer 2
    Conv2D(filters = 32, kernel_size = 5, strides = 1, activation = 'relu', kernel_regularizer=l2(0.0005), name = 'convolution_2'),
    MaxPooling2D(pool_size = 2, strides = 2, name = 'max_pool_2'),
    Dropout(0.25, name = 'dropout_2'),
        
    # Layer 3
    Flatten(name = 'flatten'),
    Dense(units = 120, activation = 'relu', name = 'fully_connected_1'),
    Dropout(0.4, name = 'dropout_3'),
        
    # Layer 4
    Dense(units = 84, activation = 'relu', name = 'fully_connected_2'),
    
    # Output
    Dense(units = 1, activation = 'sigmoid', name = 'output')
        
    ])
    
    model._name = 'LeNet5'
    

    return model

'''
This model is a ResNet model, a form of Convolutional Neural Network
Model employed currently in the classification step
'''
def regularized_padded_conv(*args, **kwargs):
    return tf.keras.layers.Conv2D(*args, **kwargs, padding='same', kernel_regularizer=_regularizer,
                                  kernel_initializer='he_normal', use_bias=False)


def bn_relu(x):
    x = tf.keras.layers.BatchNormalization()(x)
    return tf.keras.layers.ReLU()(x)


def shortcut(x, filters, stride, mode):
    if x.shape[-1] == filters:
        return x
    elif mode == 'B':
        return regularized_padded_conv(filters, 1, strides=stride)(x)
    elif mode == 'B_original':
        x = regularized_padded_conv(filters, 1, strides=stride)(x)
        return tf.keras.layers.BatchNormalization()(x)
    elif mode == 'A':
        return tf.pad(tf.keras.layers.MaxPool2D(1, stride)(x) if stride>1 else x,
                      paddings=[(0, 0), (0, 0), (0, 0), (0, filters - x.shape[-1])])
    else:
        raise KeyError("Parameter shortcut_type not recognized!")
    

def original_block(x, filters, stride=1, **kwargs):
    c1 = regularized_padded_conv(filters, 3, strides=stride)(x)
    c2 = regularized_padded_conv(filters, 3)(bn_relu(c1))
    c2 = tf.keras.layers.BatchNormalization()(c2)
    
    mode = 'B_original' if _shortcut_type == 'B' else _shortcut_type
    x = shortcut(x, filters, stride, mode=mode)
    return tf.keras.layers.ReLU()(x + c2)
    
    
def preactivation_block(x, filters, stride=1, preact_block=False):
    flow = bn_relu(x)
    if preact_block:
        x = flow
        
    c1 = regularized_padded_conv(filters, 3, strides=stride)(flow)
    if _dropout:
        c1 = tf.keras.layers.Dropout(_dropout)(c1)
        
    c2 = regularized_padded_conv(filters, 3)(bn_relu(c1))
    x = shortcut(x, filters, stride, mode=_shortcut_type)
    return x + c2


def bootleneck_block(x, filters, stride=1, preact_block=False):
    flow = bn_relu(x)
    if preact_block:
        x = flow
         
    c1 = regularized_padded_conv(filters//_bootleneck_width, 1)(flow)
    c2 = regularized_padded_conv(filters//_bootleneck_width, 3, strides=stride)(bn_relu(c1))
    c3 = regularized_padded_conv(filters, 1)(bn_relu(c2))
    x = shortcut(x, filters, stride, mode=_shortcut_type)
    return x + c3


def group_of_blocks(x, block_type, num_blocks, filters, stride, block_idx=0):
    global _preact_shortcuts
    preact_block = True if _preact_shortcuts or block_idx == 0 else False
    
    x = block_type(x, filters, stride, preact_block=preact_block)
    for i in range(num_blocks-1):
        x = block_type(x, filters)
    return x


def Resnet(input_shape, n_classes, l2_reg=1e-4, group_sizes=(2, 2, 2), features=(16, 32, 64), strides=(1, 2, 2),
           shortcut_type='B', block_type='preactivated', first_conv={"filters": 16, "kernel_size": 3, "strides": 1},
           dropout=0, cardinality=1, bootleneck_width=4, preact_shortcuts=True):
    
    global _regularizer, _shortcut_type, _preact_projection, _dropout, _cardinality, _bootleneck_width, _preact_shortcuts
    _bootleneck_width = bootleneck_width # used in ResNeXts and bootleneck blocks
    _regularizer = tf.keras.regularizers.l2(l2_reg)
    _shortcut_type = shortcut_type # used in blocks
    _cardinality = cardinality # used in ResNeXts
    _dropout = dropout # used in Wide ResNets
    _preact_shortcuts = preact_shortcuts
    
    block_types = {'preactivated': preactivation_block,
                   'bootleneck': bootleneck_block,
                   'original': original_block}
    
    selected_block = block_types[block_type]
    inputs = tf.keras.layers.Input(shape=input_shape)
    flow = regularized_padded_conv(**first_conv)(inputs)
    
    if block_type == 'original':
        flow = bn_relu(flow)
    
    for block_idx, (group_size, feature, stride) in enumerate(zip(group_sizes, features, strides)):
        flow = group_of_blocks(flow,
                               block_type=selected_block,
                               num_blocks=group_size,
                               block_idx=block_idx,
                               filters=feature,
                               stride=stride)
    
    if block_type != 'original':
        flow = bn_relu(flow)
    
    flow = tf.keras.layers.GlobalAveragePooling2D()(flow)
    outputs = tf.keras.layers.Dense(n_classes, activation='sigmoid', kernel_regularizer=_regularizer)(flow)
    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    return model


def load_weights_func(model, model_name):
    try: model.load_weights(os.path.join('saved_models', model_name + '.tf'))
    except tf.errors.NotFoundError: print("No weights found for this model!")
    return model


def cifar_resnet20(input_shape=(36, 30, 1), block_type='original', shortcut_type='A', l2_reg=1e-4, load_weights=False):
    model = Resnet(input_shape=input_shape, n_classes=2, l2_reg=l2_reg, group_sizes=(3, 3, 3), features=(16, 32, 64),
                   strides=(1, 2, 2), first_conv={"filters": 16, "kernel_size": 3, "strides": 1}, shortcut_type=shortcut_type, 
                   block_type=block_type, preact_shortcuts=False)
    if load_weights: model = load_weights_func(model, 'cifar_resnet20')
    return model