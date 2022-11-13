from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys
import glob
import functions

def true_hotpixels():
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            darkframe_data = np.copy(darkframe["PRIMARY"].data)
            media = np.mean(darkframe_data)
            desvio = 2*np.std(darkframe_data)
            limit = media + desvio
            hot_pixels = functions.get_hotpixels(darkframe, limit)
            return functions.get_true_hotpixels(hdu, hot_pixels, media, desvio, width=2)

def make_filtered_image():
    with fits.open(sys.argv[1]) as darkframe:
        limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
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
            limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
            hot_pixels = functions.get_hotpixels(darkframe, limit)
            new_data = functions.get_filtered_image_data(hdu["PRIMARY"].data, hot_pixels, width=width)

            hdul = fits.PrimaryHDU(new_data)
            new_name = functions.name_gen("Gamma.fits")
            hdul.writeto(new_name, overwrite=True)

            with fits.open(new_name) as temp:
                final_data = functions.corrections.gammaCorrection(temp)
                img = fits.PrimaryHDU(final_data)
                img.writeto(new_name, overwrite=True)

def make_histogramSpread_from_filtered_image(width=1):
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
            hot_pixels = functions.get_hotpixels(darkframe, limit)
            new_data = functions.get_filtered_image_data(hdu["PRIMARY"].data, hot_pixels, width=width)

            hdul = fits.PrimaryHDU(new_data)
            new_name = functions.name_gen("HistogramSpread.fits")
            hdul.writeto(new_name, overwrite=True)

            with fits.open(new_name) as temp:
                final_data = functions.corrections.histogramSpread(temp)
                img = fits.PrimaryHDU(final_data)
                img.writeto(new_name, overwrite=True)

if __name__ == "__main__":
    #make_filtered_image()
    #make_stacked_image()
    #print(true_hotpixels())
    make_gamma()
    #make_histogramSpread()
    #make_gamma_from_filtered_image()
    #make_histogramSpread_from_filtered_image()
