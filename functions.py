from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys
import glob

#1
def name_gen(file_name):
    i=0
    check_name = file_name
    while exists(check_name):
        i+=1
        check_name = file_name[:-5]+"_new"+str(i)+".fits"
    print(f"Novo ficheiro gerado: {check_name}")
    return check_name

#2
def get_hotpixels(hdul, limit):
    data = hdul["PRIMARY"].data
    return np.argwhere(data >= limit)

#3
def remove_hotpixels(hdu, hot_pixels, raio=1): #Utiliza 1
    a = np.copy(hdu["PRIMARY"].data)
    x_res,y_res = a.shape
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        vizinhos = []
        for i in range(-raio, raio+1):
            for j in range(-raio, raio+1):
                    if ((line+i>=0 and collumn+j>=0 and line+i<x_res and collumn+j<y_res) and (not(i==0 and j==0))):
                        vizinhos.append(a[line+i][collumn+j])

        vizinhos = np.array(vizinhos)
        mediana = np.median(vizinhos)

        a[line][collumn] = mediana

    file_name = hdu.filename()
    new_name = name_gen(file_name)

    try:
        hdul = fits.PrimaryHDU(a)
        hdul.writeto(new_name,overwrite=True)

    except Exception as e:
        print(e)

#4
def get_filtered_image_data(hdu_data, hot_pixels, raio):
    a = np.copy(hdu_data)
    x_res,y_res = a.shape
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        vizinhos = []
        for i in range(-raio, raio+1):
            for j in range(-raio, raio+1):
                    if ((line+i>=0 and collumn+j>=0 and line+i<x_res and collumn+j<y_res) and (not(i==0 and j==0))):
                        vizinhos.append(a[line+i][collumn+j])

        vizinhos = np.array(vizinhos)
        mediana = np.median(vizinhos)

        a[line][collumn] = mediana

    return a

#5
def get_ADA_HV(hdu, hot_pixels, raio):
    a = np.copy(hdu["PRIMARY"].data)
    abs_diff_average = []
    hotpixels_values = []
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        hotpixels_values.append(a[line][collumn])
        total = 0
        n_vizinhos = 0
        for i in range(-raio, raio+1):
            for j in range(-raio, raio+1):
                    if ((line+i>=0 and collumn+j>=0 and line+i<len(a) and collumn+j<len(a[0])) and (not(i==0 and j==0))):
                        total += abs(a[line][collumn] - a[line+i][collumn+j])
                        n_vizinhos += 1
        abs_diff_average.append(total/n_vizinhos)

    return abs_diff_average, hotpixels_values

#6
def scatter_plot(hotpixels_values,average_abs_diff, k):
    plt.scatter(hotpixels_values,average_abs_diff)
    plt.xlabel('Hot values')
    plt.ylabel('Absolute difference average')
    plt.title("K = {}".format(k))
    plt.show()

#7
def get_stacked_image(darkframe, imgs_direct, std_mult=2, raio=1): #2,4
    print(f"Multiplicador de desvio = {std_mult}")
    print(f"Raio = {raio}")
    #calcular hot_pixels
    darkframe_data = np.copy(darkframe["PRIMARY"].data)
    limit = np.mean(darkframe_data) + std_mult*np.std(darkframe_data)
    hot_pixels = get_hotpixels(darkframe, limit)

    imgs = glob.glob(imgs_direct+"*.fits") #Obter lista com o path das imagens
    final_image = np.zeros(shape=darkframe_data.shape, dtype=darkframe_data.dtype) #Inicializa o np.array que vai acumular os valores das varias imagens
    for img in imgs:
        with fits.open(img) as image:
            print(f"A processar {img}...")
            final_image += get_filtered_image_data(image["PRIMARY"].data, hot_pixels, raio) #Adicionar a final_image os valores de cada imagem filtrada
    final_image = final_image//len(imgs) #Calcular a media de todas as imagens filtradas
    #Guardar os dados num ficherio fits
    hdu = fits.PrimaryHDU(final_image)
    hdu.writeto(imgs_direct+"stacked.fits", overwrite=True)

#8
def plot_ADA(hdu, darkframe, multiplicadores= [0, 0.5, 1, 1.5, 2, 2.5, 3],raio=1): #2,5
    """average_abs_diff = []
    hotpixels_values = []"""
    for i in multiplicadores:
        #calcular hot_pixels
        darkframe_copy = np.copy(darkframe["PRIMARY"].data)
        limit = np.mean(darkframe_copy) + i*np.std(darkframe_copy)
        hot_pixels = get_hotpixels(darkframe, limit)

        avg, hot_values = get_ADA_HV(hdu, hot_pixels, raio) #Obtem os AVG e os valores dos hot_pixels para cada multiplicador
        """
        average_abs_diff.append(avg)
        hotpixels_values.append(hot_values)"""
        scatter_plot(hot_values,avg, i)

def get_true_hotpixels(hdul, hot_pixels, media, desvio, raio=1):
    a = np.copy(hdul["PRIMARY"].data)
    true_hotpixels = []
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        total = 0
        n_vizinhos = 0
        for i in range(-raio, raio+1):
            for j in range(-raio, raio+1):
                    if ((line+i>=0 and collumn+j>=0 and line+i<len(a) and collumn+j<len(a[0])) and (not(i==0 and j==0))):
                        total += abs(a[line][collumn] - a[line+i][collumn+j])
                        n_vizinhos += 1
        ADA = total/n_vizinhos
        if ADA >= (media - desvio):
            true_hotpixels.append([line,collumn])

    true_hotpixels = np.array(true_hotpixels)

    return true_hotpixels
