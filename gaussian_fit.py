import scipy.optimize as opt
import numpy as np
import math
from matplotlib import pyplot as plt
from plot_slider import PageSlider

def get_gaussian_params(sigma_x, sigma_y, theta):
    theta = math.radians(theta)
    a = np.cos(theta)**2 / (2 * sigma_x**2) + np.sin(theta)**2 / (2 * sigma_y**2)
    b = -(np.sin(theta)*np.cos(theta)) / (2 * sigma_x**2) + (np.sin(theta)*np.cos(theta)) / (2 * sigma_y**2)
    c = np.sin(theta)**2 / (2 * sigma_x**2) + np.cos(theta)**2 / (2 * sigma_y**2)
    return (a, b, c)


def gaussian_2d(x_y, amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    x, y = x_y
    a, b, c = get_gaussian_params(sigma_x, sigma_y, theta)
    gaussian = offset + amplitude*np.exp(-(a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2))
    return gaussian.ravel()


def best_fit(im, guess):
    h, w = im.shape
    data = im.ravel()

    x = np.linspace(0.5, w-0.5, w)
    y = np.linspace(0.5, h-0.5, h)
    x, y = np.meshgrid(x, y)

    initial_guess = (im.max()-im.min(), guess[0], guess[1], 2, 2, 0, im.min())
    bounds = ([0, 0, 0, 0, 0, 0, 0], [255, w, h, w/2, h/2, 360, 255])

    popt, _ = opt.curve_fit(gaussian_2d, (x, y), data, p0=initial_guess, bounds=bounds) 
    return popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]


def plot_contours(params):
    num_pages = len(params)

    fig, ax = plt.subplots(num="Gaussian fits")
    fig.subplots_adjust(bottom=0.18)
    
    im, title, (amp, x0, y0, sigma_x, sigma_y, theta) = params[0]
    a, b, c = get_gaussian_params(sigma_x, sigma_y, theta)
    x = np.linspace(0.5, im.shape[0]-0.5, im.shape[0])
    y = np.linspace(0.5, im.shape[1]-0.5, im.shape[1])
    x, y = np.meshgrid(x,y)
    x0 -= 0.5
    y0 -= 0.5
    z = amp * np.e ** (-0.5 * (a*(x-[x0])**2 + 2*b*(x-[x0])*(y-[y0]) + c*(y-[y0])**2))
    im = ax.imshow(im)
    ax.set_title(title)
    ax.contour(x, y, z, levels=[0, 0.15*amp, 0.3*amp, 0.45*amp, 0.6*amp, 0.75*amp, 0.9*amp], colors="red")

    ax_slider = fig.add_axes([0.1, 0.05, 0.8, 0.04])
    slider = PageSlider(ax_slider, 'Page', num_pages, activecolor="orange")

    def update(val):
        i = int(slider.val)
        im, title, (amp, x0, y0, sigma_x, sigma_y, theta) = params[i]
        x0 -= 0.5
        y0 -= 0.5
        a, b, c = get_gaussian_params(sigma_x, sigma_y, theta)
        x = np.linspace(0.5, im.shape[0]-0.5, im.shape[0])
        y = np.linspace(0.5, im.shape[1]-0.5, im.shape[1])
        x, y = np.meshgrid(x,y)
        z = amp * np.e ** (-0.5 * (a*(x-[x0])**2 + b*(x-[x0])*(y-[y0]) + c*(y-[y0])**2))
        ax.clear()
        im = ax.imshow(im)
        ax.set_title(title)
        ax.contour(x, y, z, levels=[0, 0.15*amp, 0.3*amp, 0.45*amp, 0.6*amp, 0.75*amp, 0.9*amp], colors="red")

    slider.on_changed(update)
    plt.show()

