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
    print(f"Novo ficheiro gerado: {check_name}")
    return check_name


def get_hotpixels(hdul, limit):
    data = hdul["PRIMARY"].data
    return np.argwhere(data >= limit)


def remove_hotpixels(hdu, hot_pixels, width=1):
    a = np.copy(hdu["PRIMARY"].data)
    x_res,y_res = a.shape
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        neighbours = np.array([a[line+i][collumn+j] for i in range(-width, width+1) for j in range(-width, width+1) if ((line+i>=0 and collumn+j>=0 and line+i<x_res and collumn+j<y_res) and (not(i==0 and j==0)))])
        mediana = np.median(neighbours)
        a[line][collumn] = mediana

    file_name = hdu.filename()
    new_name = name_gen(file_name)

    try:
        hdul = fits.PrimaryHDU(a)
        hdul.writeto(new_name,overwrite=True)

    except Exception as e:
        print(e)


def get_filtered_image_data(hdu_data, hot_pixels, width):
    a = np.copy(hdu_data)
    x_res,y_res = a.shape
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        neighbours = np.array([a[line+i][collumn+j] for i in range(-width, width+1) for j in range(-width, width+1) if ((line+i>=0 and collumn+j>=0 and line+i<x_res and collumn+j<y_res) and (not(i==0 and j==0)))])
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

def get_true_hotpixels(hdul, hot_pixels, media, desvio, width=1):
    a = np.copy(hdul["PRIMARY"].data)
    true_hotpixels = []
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]

        neighbours = np.array([abs(a[line][collumn] - a[line+i][collumn+j]) for i in range(-width, width+1) for j in range(-width, width+1) if ((line+i>=0 and collumn+j>=0 and line+i<len(a) and collumn+j<len(a[0])) and (not(i==0 and j==0)))])
        total = np.sum(neighbours)
        n_neighbours = len(neighbours)

        ADA = total/n_neighbours
        if ADA >= (media - desvio):
            true_hotpixels.append([line,collumn])

    true_hotpixels = np.array(true_hotpixels)

    return true_hotpixels
