from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys
import glob
import functions
"""
name_gen(file_name)
get_hotpixels(hdul, limit)
remove_hotpixels(hdu, hot_pixels, raio=1)
get_filtered_image_data(hdu_data, hot_pixels, raio=1)
get_ADA_HV(hdu, hot_pixels, raio=1)
get_stacked_image(darkframe, imgs_direct, std_mult=2, raio=1)
get_true_hotpixels(hdul, hot_pixels, media, desvio, raio=1)
scatter_plot(hotpixels_values,average_abs_diff)
"""

def true_hotpixels():
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            darkframe_data = np.copy(darkframe["PRIMARY"].data)
            media = np.mean(darkframe_data)
            desvio = 2*np.std(darkframe_data)
            limit = media + desvio
            hot_pixels = functions.get_hotpixels(darkframe, limit)
            return functions.get_true_hotpixels(hdu, hot_pixels, media, desvio, raio=2)

if __name__ == "__main__":
    hot_pixels = true_hotpixels()
    """with fits.open(sys.argv[2]) as darkframe:
        functions.get_stacked_image(darkframe, sys.argv[3], raio=5)"""
    """with fits.open(sys.argv[2]) as hdu:
        functions.remove_hotpixels(hdu, hot_pixels, raio=2)"""
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            #functions.plot_ADA(hdu, darkframe, multiplicadores=[1.5, 2, 2.5],raio=2)
            darkframe_data = np.copy(darkframe["PRIMARY"].data)
            limit = np.mean(darkframe_data) + 2*np.std(darkframe_data)
            hot_pixels = functions.get_hotpixels(darkframe, limit)

            avg, hot_values = functions.get_ADA_HV(hdu, hot_pixels, raio=2)

            for i in range(10):
                print(f"{hot_values[i]} {avg[i]}")
