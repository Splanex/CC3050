from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys
import glob
import functions

def filter_darkframe():
    with fits.open(sys.argv[1]) as darkframe:
        temp[temp < 494] = 0
        img = fits.PrimaryHDU(temp)
        img.writeto("test18.fits", overwrite=True)

def test_make_filtered_image():
    with fits.open(sys.argv[1]) as darkframe:
        hot_pixels = np.argwhere(darkframe["PRIMARY"].data != 0)
        with fits.open(sys.argv[2]) as hdu:
            functions.remove_hotpixels(hdu, hot_pixels, width=2)

def make_filtered_image():
    with fits.open(sys.argv[1]) as darkframe:
        #limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
        limit = np.mean(darkframe["PRIMARY"].data) + 0.057*np.std(darkframe["PRIMARY"].data)
        hot_pixels = functions.get_hotpixels(darkframe,limit)
        with fits.open(sys.argv[2]) as hdu:
            functions.remove_hotpixels(hdu, hot_pixels, width=2)

def make_stacked_image():
    with fits.open(sys.argv[1]) as darkframe:
        directory = sys.argv[2]
        if directory[-1] != '/':
            directory += '/'
            functions.get_stacked_image(darkframe, directory, width=2)

def make_gamma():
    with fits.open(sys.argv[1]) as hdu:
            final_data = functions.corrections.gammaCorrection(hdu)
            img = fits.PrimaryHDU(final_data)
            img.writeto(functions.name_gen("Gamma.fits"), overwrite=True)

def make_histogramSpread():
    with fits.open(sys.argv[1]) as hdu:
            final_data = functions.corrections.histogramSpread(hdu)
            img = fits.PrimaryHDU(final_data)
            img.writeto(functions.name_gen("HistogramSpread.fits"), overwrite=True)

def make_gamma_from_filtered_image(width=1):
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            #limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
            limit = np.mean(darkframe["PRIMARY"].data) + 0.057*np.std(darkframe["PRIMARY"].data)
            hot_pixels = functions.get_hotpixels(darkframe, limit)
            new_data = functions.get_filtered_image_data(hdu["PRIMARY"].data, hot_pixels, width=width)

            hdul = fits.PrimaryHDU(new_data)
            new_name = functions.name_gen("Gamma_from_filtered_image.fits")
            hdul.writeto(new_name, overwrite=True)

            with fits.open(new_name) as temp:
                final_data = functions.corrections.gammaCorrection(temp)
                img = fits.PrimaryHDU(final_data)
                img.writeto(new_name, overwrite=True)

def make_histogramSpread_from_filtered_image(width=1):
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            #limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
            limit = np.mean(darkframe["PRIMARY"].data) + 0.057*np.std(darkframe["PRIMARY"].data)
            hot_pixels = functions.get_hotpixels(darkframe, limit)
            new_data = functions.get_filtered_image_data(hdu["PRIMARY"].data, hot_pixels, width=width)

            hdul = fits.PrimaryHDU(new_data)
            new_name = functions.name_gen("HistogramSpread_from_filtered_image.fits")
            hdul.writeto(new_name, overwrite=True)

            with fits.open(new_name) as temp:
                final_data = functions.corrections.histogramSpread(temp)
                img = fits.PrimaryHDU(final_data)
                img.writeto(new_name, overwrite=True)

def remove_resistant():
    with fits.open(sys.argv[1]) as hdu:
        functions.remove_resistant_hotpixels(hdu)

def remove_all_resistant():
    with fits.open(sys.argv[1]) as hdu:
        functions.remove_all_resistant_hotpixels(hdu)

def PSNR():
    with fits.open(sys.argv[1]) as hdu1:
        with fits.open(sys.argv[2]) as hdu2:
            hdu1_copy = np.copy(hdu1["PRIMARY"].data)
            hdu2_copy = np.copy(hdu2["PRIMARY"].data)
            print(functions.get_PSNR(hdu1_copy,hdu2_copy))

def make_noisy_file():
    with fits.open(sys.argv[1]) as hdu:
        data_copy = np.copy(hdu["PRIMARY"].data)
        background_value = functions.get_background_value(data_copy)
        x_res, y_res = data_copy.shape

        x_index, y_index, noise = functions.gen_noise(background_value, x_res, y_res)

        noisy_image_data = functions.gen_noisy_image_data(x_index, y_index, noise, data_copy)

        new_name = "Noisy_image.fits"

        try:
            hdul = fits.PrimaryHDU(noisy_image_data)
            hdul.writeto(new_name,overwrite=True)
            print(f"Novo ficheiro gerado: {new_name}")

        except Exception as e:
            print(e)

if __name__ == "__main__":
    #make_filtered_image() X
    #test_make_filtered_image() X
    #make_stacked_image() X
    #make_gamma() X
    #make_histogramSpread() X
    #make_gamma_from_filtered_image()  X
    #make_histogramSpread_from_filtered_image() X
    #filter_darkframe() X
    #remove_resistant() X
    #remove_all_resistant()X
    #PSNR()
    #make_noisy_file()
