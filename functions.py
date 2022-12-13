from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys
import glob
import concurrent.futures
import corrections

def name_gen(file_name):
    i=0
    check_name = file_name
    while exists(check_name):
        i+=1
        check_name = file_name[:-5]+"_new"+str(i)+".fits"
    return check_name

def get_hotpixels(hdul, limit):
    data = hdul["PRIMARY"].data
    return np.argwhere(data >= limit)

def get_neighbours(data,x_res,y_res,line,collumn,width):
    return (np.array([data[line+i][collumn+j] for i in range(-width, width+1) for j in range(-width, width+1) if ((line+i>=0 and collumn+j>=0 and line+i<x_res and collumn+j<y_res) and (not(i==0 and j==0)))]))

def scatter_plot(hotpixels_values,average_abs_diff, k):
    plt.scatter(hotpixels_values,average_abs_diff)
    plt.xlabel('Hot values')
    plt.ylabel('Absolute difference average')
    plt.title("K = {}".format(k))
    plt.show()

def remove_hotpixels(hdu, hot_pixels, width=1):
    a = np.copy(hdu["PRIMARY"].data)
    x_res,y_res = a.shape
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        neighbours = get_neighbours(a,x_res,y_res,line,collumn,width)
        mediana = np.median(neighbours)
        a[line][collumn] = mediana

    file_name = hdu.filename()
    new_name = name_gen(file_name)

    try:
        hdul = fits.PrimaryHDU(a)
        hdul.writeto(new_name,overwrite=True)
        print(f"Novo ficheiro gerado: {new_name}")

    except Exception as e:
        print(e)

def get_filtered_image_data(hdu_data, hot_pixels, width):
    a = np.copy(hdu_data)
    x_res,y_res = a.shape
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        neighbours = get_neighbours(a,x_res,y_res,line,collumn,width)
        mediana = np.median(neighbours)
        a[line][collumn] = mediana

    return a

def get_stacked_image(darkframe, imgs_direct, std_mult=2, width=1):
    print(f"Multiplicador de desvio = {std_mult}")
    print(f"Width = {width}")
    #calcular hot_pixels
    darkframe_data = np.copy(darkframe["PRIMARY"].data)
    limit = np.mean(darkframe_data) + std_mult*np.std(darkframe_data)
    hot_pixels = get_hotpixels(darkframe, limit)

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
    print(f"Ficheiro stacked gerado: stacked.fits")


def get_ADA(value, neighbours):
    total = 0
    for neighbour in neighbours:
        if (value >= neighbour):
            total += value - neighbour
        else:
            total += neighbour - value
    return (total/len(neighbours))


def remove_resistant_hotpixels(hdu,width=1):
    data_copy = np.copy(hdu["PRIMARY"].data)
    x_res,y_res = data_copy.shape

    #Obter lista de potenciais pixeis quentes não presentes no darkframe
    limit = np.mean(hdu["PRIMARY"].data) + 0.057*np.std(hdu["PRIMARY"].data)
    hot_pixels = get_hotpixels(hdu, limit)

    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        print(f"{line}:{collumn}")
        neighbours = get_neighbours(data_copy,x_res,y_res,line,collumn,width)
        average_abs_diff = get_ADA(data_copy[line][collumn], neighbours)
        if (average_abs_diff >= 100):
            mediana = np.median(neighbours)
            data_copy[line][collumn] = mediana

    file_name = hdu.filename()
    new_name = name_gen(file_name)

    try:
        hdul = fits.PrimaryHDU(data_copy)
        hdul.writeto(new_name,overwrite=True)
        print(f"Novo ficheiro gerado: {new_name}")

    except Exception as e:
        print(e)

def remove_all_resistant_hotpixels(hdu,width=1):
    data_copy = np.copy(hdu["PRIMARY"].data)
    x_res,y_res = data_copy.shape

    all_pixels = np.argwhere(data_copy >= 0)

    for k in all_pixels:
        line = k[0]
        collumn = k[1]
        print(f"{line}:{collumn}")
        neighbours = get_neighbours(data_copy,x_res,y_res,line,collumn,width)
        average_abs_diff = get_ADA(data_copy[line][collumn], neighbours)
        if (average_abs_diff >= 100):
            mediana = np.median(neighbours)
            data_copy[line][collumn] = mediana

    file_name = hdu.filename()
    new_name = name_gen(file_name)

    try:
        hdul = fits.PrimaryHDU(data_copy)
        hdul.writeto(new_name,overwrite=True)
        print(f"Novo ficheiro gerado: {new_name}")

    except Exception as e:
        print(e)
