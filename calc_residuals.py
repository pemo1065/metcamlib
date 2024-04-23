import math
from astropy.coordinates import angular_separation
from math import degrees, radians
import find_stars
import numpy as np
from RMS.Formats import Platepar, StarCatalog
from RMS.Astrometry.ApplyAstrometry import xyToRaDecPP, raDecToXYPP
from RMS.Astrometry.Conversions import date2JD
from matplotlib import pyplot

def residual_pixels(x1, y1, x2, y2):
    return math.sqrt((x1-x2)**2 + (y1-y2)**2)

def residual_angle(d1, d2, ra1, ra2):
    return math.acos(math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(ra1-ra2))

def deg_to_rad(deg):
    return deg/360*2*math.pi

max_dist_px = 10

catalog_stars = []

def pair_stars(catalog_stars, image_stars):
    pairings = []
    print("Image stars: %s" % len(image_stars))
    print("catalog stars: %s" % len(catalog_stars))
    for ims in image_stars:
        min_dist = 1000
        best_match = None
        for cs in catalog_stars:
            dist = residual_pixels(ims[0], ims[1], cs[0], cs[1])
            if dist < min_dist:
                min_dist = dist
                best_match = cs
        if min_dist <= max_dist_px:
            pairings.append({
                "image_star": ims,
                "catalog_star": best_match,
                "residual_px": min_dist,
                "residual_x": best_match[0] - ims[0],
                "residual_y": best_match[1] - ims[1],
                "residual_angle": degrees(angular_separation(radians(best_match[2]), radians(best_match[3]), radians(ims[2]), radians(ims[3]))) #degrees(residual_angle(radians(best_match[2]), radians(ims[2]), radians(best_match[3]), radians(ims[3])))
            })
        else:
            print("Found no match for image star at coordinates (%s, %s)" % (ims[0], ims[1]))

    return pairings

def plot_residuals(pairs, im_h, im_w):
    fig, ((ax1, ax2)) = pyplot.subplots(1, 2, num="Residuals")
    ax1.set_xlabel("Radius")
    ax1.set_ylabel("Residual (deg)")
    ax2.set_xlabel("X residual (px)")
    ax2.set_ylabel("Y residual (px)")
    center = (im_h/2+0.5, im_w/2+0.5)
    residual_angles = []
    radii = []
    res_x = []
    res_y = []

    sum_x = 0
    sum_y = 0
    for pair in pairs:
        ims = pair["image_star"]
        radii.append(math.sqrt((ims[0] - center[0])**2 + (ims[1] - center[1])**2))
        residual_angles.append(pair["residual_angle"])
        res_x.append(pair["residual_x"])
        res_y.append(pair["residual_y"])
        sum_x += pair["residual_x"]
        sum_y += pair["residual_y"]
    px_avg = (sum_x/len(pairs), sum_y/len(pairs))
    ax1.scatter(radii, residual_angles)
    ax2.scatter(res_x, res_y)
    ax2.scatter([px_avg[0]], [px_avg[1]], c="r", marker="x")
    fig.show()

def get_normalized_vector(x1, y1, x2, y2):
    x_step = x1 - x2
    y_step = y1 - y2
    normalization_factor = math.sqrt(x_step**2 + y_step**2)
    x_step /= normalization_factor
    y_step /= normalization_factor
    return (x_step, y_step)

def calibrate(platepar, cat_stars, img_stars):
    platepar.fitAstrometry(platepar.JD, img_stars, cat_stars, first_platepar_fit=True, fit_only_pointing=False, fixed_scale=False)
    platepar.fitAstrometry(platepar.JD, img_stars, cat_stars, first_platepar_fit=False, fit_only_pointing=False, fixed_scale=False)
    return

def show_residuals(image_file, pairs):
    pyplot.figure("Sky image")
    img = pyplot.imread(image_file)
    pyplot.imshow(img)
    for pair in pairs:
        img_star = pair["image_star"]
        cat_star = pair["catalog_star"]
        x_step, y_step = get_normalized_vector(img_star[0], img_star[1], cat_star[0], cat_star[1])

        # Subtract 0.5 to account for pyplot having the grid points in the center of the pixels,
        # while the rest of the script assumes the grid goes along the pixel edges, i.e. the top
        # left corner of the image having coordinates (0, 0)
        x_start = cat_star[0] - 0.5
        y_start = cat_star[1] - 0.5

        # Amplify the residual error for clearer visualization

        res_line_length = pair["residual_angle"]*500
        print("[%s, %s] Residual angle: %s, length: %s" % (x_start, y_start, pair["residual_angle"], res_line_length))
        pyplot.plot([x_start, x_start + x_step*res_line_length], [y_start, y_start + y_step*res_line_length], "white")

    img_stars_x = [p["image_star"][0]-0.5 for p in pairs]
    img_stars_y = [p["image_star"][1]-0.5 for p in pairs]
    cat_stars_x = [p["catalog_star"][0]-0.5 for p in pairs]
    cat_stars_y = [p["catalog_star"][1]-0.5 for p in pairs]
    pyplot.scatter(img_stars_x, img_stars_y, c='r', marker='.', linewidths=1)
    pyplot.scatter(cat_stars_x, cat_stars_y, c='b', marker='x', linewidths=1)

    pyplot.show()

def find_closest_cat_match(star_list, cat_stars):
    for i in range(len(star_list)):
        star = star_list[i]
        min_dist = 1000
        best_match = None
        ra, dec = star[2], star[3]
        for cstar in cat_stars:
            cra, cdec = cstar[0], cstar[1]
            dist = residual_pixels(ra, dec, cra, cdec)
            if dist < min_dist:
                min_dist = dist
                best_match = cstar
        star_list[i] = (star[0], star[1], best_match[0], best_match[1])
    return star_list
        

def calc_residuals(platepar, image_file, catalog, star_list):
    image_stars = find_stars.detect(image_file, "mask.png", star_list)
    print("Image stars: %s" % image_stars)
    star_list = find_closest_cat_match(star_list, catalog)
    calibrate(platepar, np.array([(np.float64(ims[2]), np.float64(ims[3]), 0) for ims in image_stars]), np.array([(ims[0], ims[1], 0) for ims in image_stars]))
    image_stars = find_stars.detect(image_file, "mask.png")
    cra = np.array([cs[0] for cs in catalog])
    cdec = np.array([cs[1] for cs in catalog])
    #jd = date2JD(2024, 3, 7, 21, 5, 45)
    jd = platepar.JD
    x, y = raDecToXYPP(cra, cdec, jd, platepar)

    #_, cra_data, cdec_data, _ = xyToRaDecPP(len(x)*[jd], x, y, len(x)*[1], platepar, extinction_correction=False, jd_time=True)
    #for i in range(len(cra)):
    #    if x[i] < 0 or x[i] > 1920 or y[i] < 0 or y[i] > 1080:
    #        continue
    #    print("Orig: [%s, %s], New: [%s, %s], Diff: [%s, %s]" % (cra[i], cdec[i], cra_data[i], cdec_data[i], abs(cra[i]-cra_data[i]), abs(cdec[i]-cdec_data[i])))
    #exit()
    catalog_stars = list(zip(x, y, cra, cdec))

    catalog_stars = [(s[0], s[1], s[2], s[3]) for s in catalog_stars if s[0] >= 0 and s[0] <= 1920 and s[1] >= 0 and s[1] <= 1080]
    print("Stars: %s" % len(catalog_stars))

    if image_stars[0][3] == 0:
        ix = [ims[0] for ims in image_stars]
        iy = [ims[1] for ims in image_stars]
        _, ra_data, dec_data, _ = xyToRaDecPP(len(ix)*[jd], ix, iy, len(ix)*[1], platepar, extinction_correction=False, jd_time=True)

        image_stars = list(zip(ix, iy, ra_data, dec_data))

    pairs = pair_stars(catalog_stars, image_stars)

    rms_sum = 0
    for p in pairs:
        ims = p["image_star"]
        cs = p["catalog_star"]
        rms_sum += (p["residual_angle"]*60)**2
        print("pair: [%s, %s, %s, %s], [%s, %s, %s, %s]. Angular res: %s, Px res: %s" % (ims[0], ims[1], ims[2], ims[3], cs[0], cs[1], cs[2], cs[3], p["residual_angle"], p["residual_px"]))
    print("Found %s pairs" % len(pairs))
    rms = math.sqrt(rms_sum/len(pairs))
    print("RMS: %s arcmin" % rms)
    plot_residuals(pairs, 1920, 1080)
    show_residuals(image_file, pairs)

if __name__ == "__main__":
    platepar = Platepar.Platepar()
    platepar.read("parameters.json", "json")
    star_list = [(s[1], s[2], s[4], s[5]) for s in platepar.star_list]
    print("Star list: %s" % star_list)
    catalog, _, _ = StarCatalog.readStarCatalog(".", "BSC5.cat", 5.0)
    calc_residuals(platepar, "image.png", catalog, star_list)
