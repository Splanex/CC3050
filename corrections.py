from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys
import glob
import concurrent.futures

def axu_gammaCorrection(data, gamma, max):
    x_res, y_res = data.shape
    for line in range(x_res):
        for collumn in range(y_res):
            data[line][collumn] = max * ((data[line][collumn] / max) ** (1/gamma)) #data[line][collumn] ** gamma
    return data


def gammaCorrection(hdu, gamma=0.4):
    data = np.copy(hdu["PRIMARY"].data)
    max = data.max()
    k = int(len(data[0])/12)
    inicio = 0
    fim = k
    data_final = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        process1 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process2 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process3 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process4 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process5 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process6 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process7 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process8 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process9 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process10 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process11 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        inicio += k
        fim += k
        process12 = executor.submit(axu_gammaCorrection, data[inicio:fim], gamma, max)
        data_final = np.concatenate((process1.result(), process2.result(),process3.result(),process4.result(),process5.result(),process6.result(),process7.result(),process8.result(),process9.result(),process10.result(),process11.result(),process12.result()), axis=0)
        return data_final


def aux_histogramSpread(data, blackPoint, whitePoint, max):
    x_res, y_res = data.shape
    for line in range(x_res):
        for collumn in range(y_res):
            if (data[line][collumn] < blackPoint):
                z = 0
            elif (data[line][collumn] > whitePoint):
                z = 1
            else:
                z = (data[line][collumn] - blackPoint)  / (whitePoint - data[line][collumn])
            data[line][collumn] = z * max
    return data


def histogramSpread(hdu, blackPoint=0, whitePoint=65535):
    data = np.copy(hdu["PRIMARY"].data)
    max = data.max()
    k = int(len(data[0])/12)
    inicio = 0
    fim = k
    data_final = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        process1 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process2 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process3 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process4 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process5 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process6 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process7 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process8 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process9 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process10 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process11 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        inicio += k
        fim += k
        process12 = executor.submit(aux_histogramSpread, data[inicio:fim], blackPoint, whitePoint, max)
        data_final = np.concatenate((process1.result(), process2.result(),process3.result(),process4.result(),process5.result(),process6.result(),process7.result(),process8.result(),process9.result(),process10.result(),process11.result(),process12.result()), axis=0)
        return data_final
