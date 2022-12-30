from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys
import glob
import concurrent.futures
import corrections
import math
import random

#Returns a new file name based on the input file name by adding "_newk" at the end, with k being the smallest possible number
def name_gen(file_name):
    i=0
    check_name = file_name
    while exists(check_name):
        i+=1
        check_name = file_name[:-5]+"_new"+str(i)+".fits"
    return check_name

#Returns a lol with the coordinates of the pixels that go over a certain value
def get_hotpixels(data, k=0.057):
    limit = np.mean(data) + k*np.std(data)
    return np.argwhere(data >= limit)

#Returns an array with the values of the neighbouring pixels
def get_neighbours(data,x_res,y_res,line,collumn,width):
    return (np.array([data[line+i][collumn+j] for i in range(-width, width+1) for j in range(-width, width+1) if ((line+i>=0 and collumn+j>=0 and line+i<x_res and collumn+j<y_res) and (not(i==0 and j==0)))]))

#Returns a modified lol where the hot pixels' value have been changed by the median of the neighbouring pixels'
def remove_hotpixels(data, hot_pixels, width=1):
    x_res,y_res = data.shape
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        neighbours = get_neighbours(data,x_res,y_res,line,collumn,width)
        mediana = np.median(neighbours)
        data[line][collumn] = mediana


def remove_hotpixels_0(data, darkframe_data, hot_pixels, width=1):
    x_res,y_res = data.shape
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        data[line][collumn] = max(0,data[line][collumn]-darkframe_data[line][collumn])


#Returns the absolute difference average between a pixel and its neighbours
def get_ADA(value, neighbours):
    total = 0
    for neighbour in neighbours:
        if (value >= neighbour):
            total += value - neighbour
        else:
            total += neighbour - value
    return (total/len(neighbours))


#Goes through every pixel in a list of potential hot pixels and, if its absolute difference average is above alpha, changes its value with the median of the neighbouring pixels'
def remove_resistant_hotpixels(data_copy, min=500,width=1):
    x_res,y_res = data_copy.shape

    hot_pixels = get_hotpixels(data_copy)
    #hot_pixels = get_background_value(data_copy)

    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        neighbours = get_neighbours(data_copy,x_res,y_res,line,collumn,width)
        average_abs_diff = get_ADA(data_copy[line][collumn], neighbours)
        if (average_abs_diff >= min):
            mediana = np.median(neighbours)
            data_copy[line][collumn] = mediana

    return data_copy

#Goes through every pixel and, if its absolute difference average is above alpha, changes its value with the median of the neighbouring pixels'
def remove_all_resistant_hotpixels(data_copy,width=1, alpha=100):
    x_res,y_res = data_copy.shape
    all_pixels = np.argwhere(data_copy >= 0)

    for k in all_pixels:
        line = k[0]
        collumn = k[1]
        neighbours = get_neighbours(data_copy,x_res,y_res,line,collumn,width)
        average_abs_diff = get_ADA(data_copy[line][collumn], neighbours)
        if (average_abs_diff >= alpha):
            mediana = np.median(neighbours)
            data_copy[line][collumn] = mediana

    return data_copy

#Calculates and returns the peak signal to noise ratio between an image with no hot pixels and a copy of it with noise inserted into it that went through the hot pixels removing algorithm
def get_PSNR(data1, data2, max=65535):
    x_res, y_res = data1.shape
    temp_mse = 0
    for x in range(0, x_res):
        for y in range(0, y_res):
            temp_mse += (float(data1[x][y]) - float(data2[x][y]))**2.0
    mse = temp_mse / (x_res*y_res)

    if mse == 0:
        return 0

    return (20*math.log(65535,10) - 10*math.log(mse,10))

def get_PSNR_0(data1, data2, hot_pixels, max=65535):
    temp_mse = 0
    for hot in hot_pixels:
        x = hot[0]
        y = hot[1]
        temp_mse += (float(data1[x][y]) - float(data2[x][y]))**2.0
    mse = temp_mse / len(hot_pixels)
    if mse == 0:
        return 0

    return (20*math.log(max,10) - 10*math.log(mse,10))

def get_background_value(data):
    counts = np.bincount(data.flatten())
    moda = np.argmax(counts)
    return moda

def gen_noise(background_value, x_res, y_res, noise_quantity=random.choice([0.005,0.006,0.007,0.008,0.009,0.01]), max_value=65535):
    num_noise = int(x_res * y_res * noise_quantity)
    possible_noise = range(background_value, max_value)
    noise_list = []
    x = []
    y = []
    for i in range(num_noise):
        noise_list.append(random.choice(possible_noise))

    for i in range(num_noise):
        x.append(random.choice(range(x_res)))
        y.append(random.choice(range(y_res)))

    noise_x_indexes = np.array(x, dtype = "uint16")
    noise_y_indexes = np.array(y, dtype = "uint16")
    noise_list = np.array(noise_list, dtype = "uint16")

    return noise_x_indexes, noise_y_indexes, noise_list

def gen_noisy_image_data(x_index, y_index, noise, perfect_image_data):
    for i in range(0,len(x_index)):
        if perfect_image_data[x_index[i]][y_index[i]] < noise[i]:
            perfect_image_data[x_index[i]][y_index[i]] = noise[i]
    return perfect_image_data

def gen_noisy_image_data_0(hot_pixels, perfect_image_data, darkframe_data, k=None):
    for hot in hot_pixels:
        x = hot[0]
        y = hot[1]
        #if x > 500 and x < 3500 and y > 500 and y < 5500:
        if k==None:
            perfect_image_data[x][y] = max(perfect_image_data[x][y],darkframe_data[x][y])
        else:
            noise = random.gauss(1,k)
            perfect_image_data[x][y] = max(perfect_image_data[x][y],noise*darkframe_data[x][y])
    #perfect_image_data = np.array(perfect_image_data, dtype = "uint16")
    return perfect_image_data

def sub_image(data, darkframe_data):
    x_res, y_res = data.shape
    for x in range(0,x_res):
        for y in range(0,y_res):
            if (data[x][y]>darkframe_data[x][y]):
                data[x][y] = data[x][y]-darkframe_data[x][y]
            else:
                data[x][y] = 0

def flat_field_correction(data, flatfield_data):
    x_res, y_res = data.shape
    scale_factor = np.max(flatfield_data.flatten())

    for x in range(0,x_res):
        for y in range(0,y_res):
            data[x][y] = round((float(data[x][y])*scale_factor) / (float(flatfield_data[x][y])))
    #data = np.array(data, dtype = "uint16")


"""def scatter_plot(hotpixels_values,average_abs_diff, k):
    plt.scatter(hotpixels_values,average_abs_diff)
    plt.xlabel('Hot values')
    plt.ylabel('Absolute difference average')
    plt.title("K = {}".format(k))
    plt.show()"""

"""def get_stacked_image(darkframe, imgs_direct, std_mult=2, width=1):
    print(f"Multiplicador de desvio = {std_mult}")
    print(f"Width = {width}")
    #calcular hot_pixels
    darkframe_data = np.copy(darkframe["PRIMARY"].data)
    hot_pixels = get_hotpixels(darkframe_data)

    imgs = glob.glob(imgs_direct+"*.fits") #Obter lista com o path das imagens
    final_image = np.zeros(shape=darkframe_data.shape, dtype=darkframe_data.dtype) #Inicializa o np.array que vai acumular os valores das varias imagens
    for img in imgs:
        with fits.open(img) as image:
            print(f"A processar {img}...")
            final_image += get_filtered_image_data(image["PRIMARY"].data, hot_pixels, width) #Adicionar a final_image os valores de cada imagem filtrada
    final_image = final_image//len(imgs) #Calcular a media de todas as imagens filtradas
    #Guardar os dados num ficherio fits
    hdu = fits.PrimaryHDU(final_image)
    hdu.writeto(imgs_direct+"stacked.fits", overwrite=True)
    print(f"Ficheiro stacked gerado: stacked.fits")"""
