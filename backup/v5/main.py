from astropy.io import fits
import numpy as np
import astroalign as aa
from os.path import exists
import sys, glob, functions, scipy

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

def make_masterlight(width=1, sub=None):
    with fits.open(sys.argv[1]) as darkframe:
        darkframe_data = darkframe["PRIMARY"].data
        darkframe_data = np.array(darkframe_data, dtype = "float32")
        x_res, y_res = darkframe_data.shape
        total = np.zeros((x_res,y_res))
        total = np.array(total, dtype = "float32")
        lights_dir = sys.argv[2]
        with fits.open(sys.argv[3]) as flatfield:
            flatfield_data = np.copy(flatfield["PRIMARY"].data)
            flatfield_data = np.array(flatfield_data, dtype = "float32")
            lights = glob.glob(lights_dir+"/*.fit")
            for i in range(0,len(lights)):
                print(i)
                print(f"{(i/len(lights))*100}%")
                if i == 0: #primeiro ficheiro light
                    with fits.open(lights[i]) as hdu:
                        data_copy = np.copy(hdu["PRIMARY"].data)
                        data_copy = np.array(data_copy, dtype = "float32")

                        hot_pixels = functions.get_hotpixels(darkframe_data)
                        print("")

                        if sub == None:
                            functions.remove_hotpixels(data_copy, hot_pixels,width=width)
                            #functions.remove_resistant_hotpixels(data_copy,width=width)
                        else:
                            functions.sub_image(data_copy, darkframe_data)

                        functions.flat_field_correction(data_copy, flatfield_data)

                        total += data_copy

        total = total / len(lights)
        #total = np.round_(total)
        #total = np.array(total, dtype = "uint16")

        try:
            hdul = fits.PrimaryHDU(total)
            hdul.writeto(lights_dir+"/Masterlightnew4.fits",overwrite=True)
            print(lights_dir+"/Masterlightnew.fits")

        except Exception as e:
            print(e)

        print(np.median(total)/np.std(total))

def make_aligned_masterlight(width=1, sub=None):
    with fits.open(sys.argv[1]) as darkframe:
        darkframe_data = darkframe["PRIMARY"].data
        darkframe_data = np.array(darkframe_data, dtype = "float32")
        x_res, y_res = darkframe_data.shape
        total = np.zeros((x_res,y_res))
        total = np.array(total, dtype = "float32")
        lights_dir = sys.argv[2]
        with fits.open(sys.argv[3]) as flatfield:
            flatfield_data = np.copy(flatfield["PRIMARY"].data)
            flatfield_data = np.array(flatfield_data, dtype = "float32")
            lights = glob.glob(lights_dir+"/*.fit")
            for i in range(0,len(lights)):
                print(i)
                print(f"{(i/len(lights))*100}%")
                if i == 0: #primeiro ficheiro light
                    with fits.open(lights[i]) as hdu:
                        data_copy = np.copy(hdu["PRIMARY"].data)
                        data_copy = np.array(data_copy, dtype = "float32")

                        hot_pixels = functions.get_hotpixels(darkframe_data, k=0.9)

                        if sub == None:
                            functions.remove_hotpixels(data_copy, hot_pixels,width=width)
                            #functions.remove_resistant_hotpixels(data_copy,width=width)
                        else:
                            functions.sub_image(data_copy, darkframe_data)

                        functions.flat_field_correction(data_copy, flatfield_data)

                        total += data_copy
                else:
                    with fits.open(lights[i]) as hdu:
                        with fits.open(lights[0]) as first_hdu:
                            data_copy = np.copy(hdu["PRIMARY"].data)
                            data_copy = np.array(data_copy, dtype = "float32")

                            first_data = first_hdu["PRIMARY"].data
                            first_data = np.array(first_data, dtype = "float32")

                            hot_pixels = functions.get_hotpixels(darkframe_data, k=0.9)

                            if sub == None:
                                functions.remove_hotpixels(data_copy, hot_pixels,width=width)
                                #functions.remove_resistant_hotpixels(data_copy,width=width)
                            else:
                                functions.sub_image(data_copy, darkframe_data)

                            functions.flat_field_correction(data_copy, flatfield_data)

                            aligned_image, footprint = aa.register(data_copy, first_data)
                            aligned_image[footprint] = 0.0

                            total += aligned_image

        total = total / len(lights)
        #total = np.round_(total)
        #total = np.array(total, dtype = "uint16")

        try:
            hdul = fits.PrimaryHDU(total)
            hdul.writeto(lights_dir+"/Masterlightnew2.fits",overwrite=True)

        except Exception as e:
            print(e)

        print(np.average(total)/np.std(total))

if __name__ == "__main__":

    make_aligned_masterlight(width=2, sub=1)
    #make_masterlight(width=2, sub=1)
    #make_flatfield_corrected_image()
    #make_darkframe_subtration_image
    #evaluate()
