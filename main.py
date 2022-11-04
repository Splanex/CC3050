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
            return functions.get_true_hotpixels(hdu, hot_pixels, media, desvio, raio=2)

def make_filtered_image():
    with fits.open(sys.argv[1]) as darkframe:
        limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
        hot_pixels = functions.get_hotpixels(darkframe,limit)
        with fits.open(sys.argv[2]) as hdu:
            functions.remove_hotpixels(hdu, hot_pixels, raio=2)

def make_stacked_image():
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as direct:
            functions.get_stacked_image(darkframe, direct, raio=5)

def make_gamma():
    with fits.open(sys.argv[1]) as hdu:
            final_data = functions.gammaCorrection(hdu)
            img = fits.PrimaryHDU(final_data)
            img.writeto(functions.new_name("gamma.fits"), overwrite=True)

def make_histogramSpread():
    with fits.open(sys.argv[1]) as hdu:
            final_data = functions.histogramSpread(hdu)
            img = fits.PrimaryHDU(final_data)
            img.writeto(functions.new_name("gamma.fits"), overwrite=True)

def make_gamma_from_filtered_image(raio=1):
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
            hot_pixels = functions.get_hotpixels(darkframe, limit)
            new_data = functions.get_filtered_image_data(hdu["PRIMARY"].data, hot_pixels, raio=raio)

            hdul = fits.PrimaryHDU(new_data)
            hdul.writeto("test.fits", overwrite=True)

            with fits.open("test.fits") as temp:
                final_data = functions.gammaCorrection(temp)
                img = fits.PrimaryHDU(final_data)
                img.writeto("test.fits", overwrite=True)

def make_histogramSpread_from_filtered_image(raio=1):
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            limit = np.mean(darkframe["PRIMARY"].data) + 2*np.std(darkframe["PRIMARY"].data)
            hot_pixels = functions.get_hotpixels(darkframe, limit)
            new_data = functions.get_filtered_image_data(hdu["PRIMARY"].data, hot_pixels, raio=raio)

            hdul = fits.PrimaryHDU(new_data)
            hdul.writeto("test.fits", overwrite=True)

            with fits.open("test.fits") as temp:
                final_data = functions.histogramSpread(temp)
                img = fits.PrimaryHDU(final_data)
                img.writeto("test.fits", overwrite=True)

def make_plots(multiplicadores=[0, 0.5, 1,1.5, 2, 2.5, 3], raio=1):
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            functions.plot_ADA(darkframe, hdu)

if __name__ == "__main__":
    make_plots()
