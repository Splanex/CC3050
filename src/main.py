from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys, glob, functions

def filter_darkframe():
    with fits.open(sys.argv[1]) as darkframe:
        temp = np.copy(darkframe["PRIMARY"].data)
        temp[temp < 494] = 0
        img = fits.PrimaryHDU(temp)
        img.writeto("testX.fits", overwrite=True)

def test_make_filtered_image():
    with fits.open(sys.argv[1]) as darkframe:
        hot_pixels = np.argwhere(darkframe["PRIMARY"].data != 0)
        with fits.open(sys.argv[2]) as hdu:
            new_data = np.copy(hdu["PRIMARY"].data)
            functions.remove_hotpixels(new_data, hot_pixels, width=2)

            file_name = hdu.filename()
            new_name = functions.name_gen(file_name)

            try:
                hdul = fits.PrimaryHDU(new_data)
                hdul.writeto(new_name,overwrite=True)
                print(f"Novo ficheiro gerado: {new_name}")

            except Exception as e:
                print(e)

def make_filtered_image():
    with fits.open(sys.argv[1]) as darkframe:
        darkframe_data = darkframe["PRIMARY"].data
        hot_pixels = functions.get_hotpixels(darkframe_data)
        with fits.open(sys.argv[2]) as hdu:
            new_data = np.copy(hdu["PRIMARY"].data)
            functions.remove_hotpixels(new_data, hot_pixels, width=2)

            file_name = hdu.filename()
            new_name = functions.name_gen(file_name)

            try:
                hdul = fits.PrimaryHDU(new_data)
                hdul.writeto(new_name,overwrite=True)
                print(f"Novo ficheiro gerado: {new_name}")

            except Exception as e:
                print(e)

def make_gamma():
    with fits.open(sys.argv[1]) as hdu:
            data = np.copy(hdu["PRIMARY"].data)
            final_data = functions.corrections.gammaCorrection(data)
            img = fits.PrimaryHDU(final_data)
            img.writeto(functions.name_gen("Gamma.fits"), overwrite=True)

def make_histogramSpread():
    with fits.open(sys.argv[1]) as hdu:
            data = np.copy(hdu["PRIMARY"].data)
            final_data = functions.corrections.histogramSpread(data)
            img = fits.PrimaryHDU(final_data)
            img.writeto(functions.name_gen("HistogramSpread.fits"), overwrite=True)

def make_gamma_from_filtered_image(width=1):
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            darkframe_data = darkframe["PRIMARY"].data
            hot_pixels = functions.get_hotpixels(darkframe_data)
            new_data = hdu["PRIMARY"].data
            functions.remove_hotpixels(new_data, hot_pixels, width=width)

            hdul = fits.PrimaryHDU(new_data)
            new_name = functions.name_gen("Gamma_from_filtered_image.fits")
            hdul.writeto(new_name, overwrite=True)

            with fits.open(new_name) as temp:
                data = np.copy(temp["PRIMARY"].data)
                final_data = functions.corrections.gammaCorrection(data)
                img = fits.PrimaryHDU(final_data)
                img.writeto(new_name, overwrite=True)

def make_histogramSpread_from_filtered_image(width=1):
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            darkframe_data = darkframe["PRIMARY"].data
            hot_pixels = functions.get_hotpixels(darkframe_data)
            new_data = hdu["PRIMARY"].data
            functions.remove_hotpixels(new_data, hot_pixels, width=width)

            hdul = fits.PrimaryHDU(new_data)
            new_name = functions.name_gen("HistogramSpread_from_filtered_image.fits")
            hdul.writeto(new_name, overwrite=True)

            with fits.open(new_name) as temp:
                data = np.copy(temp["PRIMARY"].data)
                final_data = functions.corrections.histogramSpread(data)
                img = fits.PrimaryHDU(final_data)
                img.writeto(new_name, overwrite=True)

def remove_resistant():
    with fits.open(sys.argv[1]) as hdu:
        data = np.copy(hdu["PRIMARY"].data)
        functions.remove_resistant_hotpixels(data)

        file_name = hdu.filename()
        new_name = functions.name_gen(file_name)

        try:
            hdul = fits.PrimaryHDU(data)
            hdul.writeto(new_name,overwrite=True)
            print(f"Novo ficheiro gerado: {new_name}")

        except Exception as e:
            print(e)

def remove_all_resistant():
    with fits.open(sys.argv[1]) as hdu:
        data = np.copy(hdu["PRIMARY"].data)
        functions.remove_all_resistant_hotpixels(data)

        file_name = hdu.filename()
        new_name = functions.name_gen(file_name)

        try:
            hdul = fits.PrimaryHDU(data)
            hdul.writeto(new_name,overwrite=True)
            print(f"Novo ficheiro gerado: {new_name}")

        except Exception as e:
            print(e)

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

        new_name = "Noisy_image3.fits"

        try:
            hdul = fits.PrimaryHDU(noisy_image_data)
            hdul.writeto(new_name,overwrite=True)
            print(f"Novo ficheiro gerado: {new_name}")

        except Exception as e:
            print(e)

def make_noisy_file_0():
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            data_copy = np.copy(hdu["PRIMARY"].data)

            darkframe_data = darkframe["PRIMARY"].data
            hot_pixels = functions.get_hotpixels(darkframe_data)

            noisy_image_data = functions.gen_noisy_image_data_0(hot_pixels, data_copy, darkframe_data, 0.5)

            new_name = "Noisy_image3.fits"

            try:
                hdul = fits.PrimaryHDU(noisy_image_data)
                hdul.writeto(new_name,overwrite=True)
                print(f"Novo ficheiro gerado: {new_name}")

            except Exception as e:
                print(e)

def evaluate():
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            darkframe_data = darkframe["PRIMARY"].data
            hot_pixels = functions.get_hotpixels(darkframe_data)
            data = hdu["PRIMARY"].data

            data_copy = np.copy(data)

            noisy_image_data = functions.gen_noisy_image_data_0(hot_pixels, data_copy, darkframe_data, 0.5)

            data_copy2 = np.copy(noisy_image_data)
            functions.sub_image(data_copy2, darkframe_data)

            try:
                hdul = fits.PrimaryHDU(noisy_image_data)
                hdul.writeto("Noisy_image3.fits",overwrite=True)

            except Exception as e:
                print(e)

            functions.remove_hotpixels(noisy_image_data, hot_pixels, width=2)
            new_data = noisy_image_data
            functions.remove_resistant_hotpixels(new_data, min=1000)
            corrected_data = new_data

            try:
                hdul2 = fits.PrimaryHDU(corrected_data)
                hdul2.writeto("Corrected_image3.fits",overwrite=True)

            except Exception as e:
                print(e)

            try:
                hdul3 = fits.PrimaryHDU(data_copy2)
                hdul3.writeto("Sub_image3.fits",overwrite=True)

            except Exception as e:
                print(e)

            print(functions.get_PSNR_0(data,corrected_data,hot_pixels))
            print(functions.get_PSNR(data,corrected_data))
            print(functions.get_PSNR(data,data_copy2))

def make_darkframe_subtration_image():
    with fits.open(sys.argv[1]) as darkframe:
        with fits.open(sys.argv[2]) as hdu:
            darkframe_data = darkframe["PRIMARY"].data
            data = hdu["PRIMARY"].data
            data_copy = np.copy(data)
            functions.sub_image(data_copy, darkframe_data)

            try:
                hdul = fits.PrimaryHDU(data_copy)
                hdul.writeto("Sub_image3.fits",overwrite=True)

            except Exception as e:
                print(e)

def make_flatfield_corrected_image(darkframe, hdu, flatfield):
    darkframe_data = darkframe["PRIMARY"].data
    data_copy = np.copy(hdu["PRIMARY"].data)
    flatfield_data = np.copy(flatfield["PRIMARY"].data)

    hot_pixels = functions.get_hotpixels(darkframe_data, k=2)

    functions.remove_hotpixels(data_copy, hot_pixels, width=2)
    functions.flat_field_correction(data_copy, flatfield_data)

    try:
        file_name = hdu.filename()
        new_name = functions.name_gen(file_name)
        hdul = fits.PrimaryHDU(data_copy)
        hdul.writeto(new_name,overwrite=True)

    except Exception as e:
        print(e)

def make_masterlight(width=1):
    with fits.open(sys.argv[1]) as darkframe:
        darkframe_data = darkframe["PRIMARY"].data
        x_res, y_res = darkframe_data.shape
        total = np.zeros((x_res,y_res))
        lights_dir = sys.argv[2]
        with fits.open(sys.argv[3]) as flatfield:
            flatfield_data = np.copy(flatfield["PRIMARY"].data)
            lights = glob.glob(lights_dir+"/NGC*")
            for light in lights:
                with fits.open(light) as hdu:
                    data_copy = np.copy(hdu["PRIMARY"].data)
                    hot_pixels = functions.get_hotpixels(darkframe_data)

                    functions.remove_hotpixels(data_copy, hot_pixels,width=width)
                    corrected_data = data_copy
                    functions.remove_resistant_hotpixels(corrected_data,width=width)
                    functions.flat_field_correction(corrected_data, flatfield_data)
                    total += corrected_data

        print(total[0][0])
        total = total / len(lights)
        total = np.array(total, dtype = "uint16")

        try:
            hdul = fits.PrimaryHDU(total)
            hdul.writeto(lights_dir+"/masterLight.fits",overwrite=True)

        except Exception as e:
            print(e)


if __name__ == "__main__":

    with fits.open(sys.argv[1]) as perfect:
        with fits.open(sys.argv[2]) as my_perfect:
            perfect_data = perfect["PRIMARY"].data
            my_perfect_data = my_perfect["PRIMARY"].data
            print(functions.get_PSNR(perfect_data,my_perfect_data))

    #make_masterlight(width=2)
    #make_flatfield_corrected_image()

    #evaluate()

    #make_filtered_image() #
    #remove_resistant() #
    #remove_all_resistant() #
    #make_noisy_file() #
    #make_noisy_file_0() #
    #PSNR() #

    #test_make_filtered_image() #
    #filter_darkframe() #

    #make_gamma() #
    #make_histogramSpread() #
    #make_gamma_from_filtered_image() #
    #make_histogramSpread_from_filtered_image() #

    #make_stacked_image()
