from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys

def usage():
    print("""
Usage:
     -f   Indica o ficheiro dark frame a utilizar para criar o mapa.
     -t   Indica o ficheiro dark frame a utilizar para criar o mapa e é aplicada a máscara a este ficheiro.
     -d   Indica o diretorio com os ficheiros dark frame a utilizar na criação do mapa.
     -dt  Mesma função que "-t" mas para um diretorio de ficheiros.
     -i   Cria a imagem invertida.
     -m   Aplica a "mascara", removendo os pixeis quentes.

Examples: python main.py -f inputFile
          python main.py -t inputFile
          python main.py -d inputDirectory
          python main.py -dt inputDirectory
          python main.py -i inputFile
          python main.py -m darkframe originalfile
    """)
    sys.exit()

def name_gen(file_name):
    i=0
    check_name = file_name
    while exists(check_name):
        i+=1
        check_name = file_name[:-5]+"_new"+str(i)+".fits"
    print(f"Novo ficheiro gerado: {check_name}")
    return check_name

def inverted(hdu):
    file_name = hdu.filename()
    new_name = name_gen(file_name)
    n_array = np.full_like(hdu[0].data, 65535)
    n_array -= hdu[0].data
    hdu["PRIMARY"].data = n_array
    try:
        hdu.writeto(new_name)
    except Exception as e:
        print(e)

def get_hotpixels(hdul, limit):
    data = hdul["PRIMARY"].data
    return np.argwhere(data >= limit)

def remove_hotpixels(hdu, hot_pixels):

    a = hdu["PRIMARY"].data
    for k in hot_pixels:
        line = k[0]
        collumn = k[1]
        vizinhos = []
        for i in range(-1,2):
            for j in range(-1,2):
                    if ((line+i>=0 and collumn+j>=0 and line+i<len(a) and collumn+j<len(a[0])) and (not(i==0 and j==0))):
                        vizinhos.append(a[line+i][collumn+j])

        vizinhos = np.array(vizinhos)
        mediana = np.median(vizinhos)

        hdu["PRIMARY"].data[line][collumn] = mediana

    file_name = hdu.filename()
    new_name = name_gen(file_name)

    try:
        hdu.writeto(new_name)
    except Exception as e:
        print(e)

def main():

    options = ["-f", "-d", "-i", "-m", "-dt", "-t"]

    #conversao = lambda x : (1-x)*3

    if len(sys.argv) > 1 and sys.argv[1] in options:
#############################################################################################################
        if sys.argv[1] == "-f" or sys.argv[1] == "-d":
            if sys.argv[1] == "-f":
                with fits.open(sys.argv[1]) as hdul:
                    x_res,y_res = hdul["PRIMARY"].data.shape
                    n_pixels = x_res * y_res
                    limit = np.mean(hdul["PRIMARY"].data + 2*np.std(hdul["PRIMARY"].data))
                    hot_pixels = get_hotpixels(hdul, limit)
                    remove_hotpixels(hdul, hot_pixels)
                    print(f"{((len(hot_pixels)/n_pixels)*100):.2f}%")
                    hdul.close()
            else:
                pass
#############################################################################################################
        elif sys.argv[1] == "-i":
            with fits.open(sys.argv[1]) as hdul:
                inverted(hdul)
                hdul.close()
#############################################################################################################
        elif sys.argv[1] == "-m":
            with fits.open(sys.argv[3]) as hdul:
                with fits.open(sys.argv[2]) as darkframe:
                    limit = np.mean(darkframe["PRIMARY"].data + 2*np.std(darkframe["PRIMARY"].data))
                    hot_pixels = get_hotpixels(darkframe, limit)
                    remove_hotpixels(hdul, hot_pixels)
                hdul.close()
#############################################################################################################
        elif sys.argv[1] == "-dt":
            pass
#############################################################################################################
        else:
            with fits.open(sys.argv[2]) as hdul: # -t
                x_res,y_res = hdul["PRIMARY"].data.shape
                n_pixels = x_res * y_res
                limit = np.mean(hdul["PRIMARY"].data + 2*np.std(hdul["PRIMARY"].data))
                hot_pixels = get_hotpixels(hdul, limit)
                remove_hotpixels(hdul,hot_pixels)
                hdul.close()
#############################################################################################################
    else:
        usage()

if __name__ == "__main__":
    main()
