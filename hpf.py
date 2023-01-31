from astropy.io import fits
import astroalign
import math
import numpy as np
import scipy.ndimage as ndi
import sys

DEFAULT_K = 0.057

# I/O
def load(arg):
    if type(arg) == str:
        # Load from file
        return np.array(fits.getdata(arg), dtype="float64")
    if type(arg) == list:
        # Load list of images
        return [ load(a) for a in arg]
    # Already an image (it has been loaded already)
    return arg

def save(filename, image):
    if image.dtype != 'uint16':
        image = np.array(np.clip(image, 0.0, 65535.0), dtype='uint16')
    fits.PrimaryHDU(image).writeto(filename, overwrite=True)

def fits_info(file):
    print('---', file, '---')
    with fits.open(file) as hdul:
        hdul.info()
        for (h,v) in hdul[0].header.items():
            print(h,":", v)

# Elementary operations

def snr(image):
    image = load(image)
    return np.mean(image) / np.std(image)

def sub(a, b):
    a, b = load([a, b])
    return np.maximum(a - b, np.zeros(1, dtype=a.dtype))

def diff(a, b):
    a, b = load([a, b])
    return np.abs(a - b)

def average(image_list):
    return np.sum(load(image_list), axis=0) / len(image_list)

def psnr(a, b):
    a, b = load([a, b])
    mse = np.mean((a - b) ** 2)
    return (0 if mse == 0 else
        (20*math.log(65535,10) - 10*math.log(mse,10)))

# Flat field correction
def flat_field_correction(flatfield, image):
    image, flatfield = load([image, flatfield])
    return (image *  np.max(flatfield)) / flatfield

# Dark-frame correction algorithm (single image)
def darkframe_correction(darkframe, flatfield, image):
    image = sub(image, darkframe)
    image = flat_field_correction(image=image, flatfield=flatfield)
    return image

# Dark-frame correction algorithm (master-light variant)
def darkframe_correction_masterlight(darkframe, flatfield, image_list):
    darkframe, flatfield, image_list = load([darkframe, flatfield, image_list])
    return average(register([
        darkframe_correction(image=image, darkframe=darkframe, flatfield=flatfield)
        for image in image_list
    ]))

def get_hotpixels(image, k=DEFAULT_K):
    image = load(image)
    limit = np.median(image) + k * np.std(image)
    return np.ndarray.nonzero(image >= limit)

def highlight_hotpixels(image, k=DEFAULT_K):
    image = load(image)
    hhp = np.zeros(image.shape, dtype="uint16")
    """x = get_hotpixels(image, k)
                hhp[x] = 65535
                print(len(x[0]))
                print((len(x[0])/(len(image[0])*len(image[1])))*100)
                sys.exit()"""
    hhp[get_hotpixels(image, k)] = 65535
    return hhp

def hotpixel_correction(darkframe, flatfield, image, k=DEFAULT_K, width=1):
    hot_pixels = get_hotpixels(darkframe, k)
    print(hot_pixels[0])
    print(hot_pixels[1])
    sys.exit()    
    image = load(image).copy()
    median_image = ndi.median_filter(image, mode='reflect', size=width+1)
    image[hot_pixels] = median_image[hot_pixels]
    image = flat_field_correction(image=image, flatfield=flatfield)
    return image

def hotpixel_correction_masterlight(darkframe, flatfield, image_list,k=DEFAULT_K,width=1):    
    darkframe, flatfield, image_list = load([darkframe, flatfield, image_list])
    return average(register([
        hotpixel_correction(image=img, darkframe=darkframe, flatfield=flatfield,
                            k=k, width=width)
        for img in image_list
    ]))

def register(images):
    for i in range(1, len(images)):
        images[i],footprint = astroalign.register(images[i], images[0])
        (images[i])[footprint]= 0.0
    return images
