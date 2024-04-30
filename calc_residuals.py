import cv2
import ephem
import importlib
import math
from astropy.coordinates import angular_separation
from math import degrees, radians
import find_stars
import numpy as np
from RMS.Formats import StarCatalog
from matplotlib import pyplot

def residual_pixels(x1, y1, x2, y2):
    return math.sqrt((x1-x2)**2 + (y1-y2)**2)

def residual_angle(d1, d2, ra1, ra2):
    return math.acos(math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(ra1-ra2))

def deg_to_rad(deg):
    return deg/360*2*math.pi

catalog_stars = []

def pair_stars(catalog_stars, image_stars, max_dist_px=10):
    pairings = []
    print("Image stars: %s" % len(image_stars))
    print("catalog stars: %s" % len(catalog_stars))
    for ims in image_stars:
        dists = []
        min_dist = 1000
        best_match = None
        for cs in catalog_stars:
            dist = residual_pixels(ims[0], ims[1], cs[0], cs[1])
            dists.append(dist)
            if dist < min_dist:
                min_dist = dist
                best_match = cs
        dists = sorted(dists)
        if min_dist <= max_dist_px:
            #if dists[1]/dists[0] < 1.5:
            #    print("Ambiguous pairing for star at (%s, %s). Found two close matches with dists %s and %s. Skipping" % (ims[0], ims[1], dists[0], dists[1]))
            #    continue
            pairings.append({
                "image_star": ims,
                "catalog_star": best_match,
                "residual_px": min_dist,
                "residual_x": best_match[0] - ims[0],
                "residual_y": best_match[1] - ims[1],
                "residual_angle": degrees(angular_separation(radians(best_match[2]), radians(best_match[3]), radians(ims[2]), radians(ims[3]))) #degrees(residual_angle(radians(best_match[2]), radians(ims[2]), radians(best_match[3]), radians(ims[3])))
            })
        else:
            print("Found no match for image star at coordinates (%s, %s). Min dist: %s" % (ims[0], ims[1], min_dist))

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

def show_residuals(image, pairs, catalog_stars, image_stars, mask_file="mask.png"):
    mask = None
    if mask_file is not None:
        mask = cv2.imread(mask_file, 0)
    pyplot.figure("Sky image")
    image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    pyplot.imshow(image)
    pyplot.subplots_adjust(left=0.05, right=0.995, top=0.995, bottom=0.05)
    print("\n== Matched stars and residuals ==")
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
        res_line_length = pair["residual_px"]*40
        print("  [%s, %s] Residual angle: %s, residual pixels: %s" % (x_start, y_start, pair["residual_angle"], pair["residual_px"]))
        pyplot.plot([x_start, x_start + x_step*res_line_length], [y_start, y_start + y_step*res_line_length], "lightgray")
    img_stars_x = []
    img_stars_y = []
    for ims in image_stars:
        paired = False
        for p in pairs:
            if p["image_star"][0] == ims[0] and p["image_star"][1] == ims[1]:
                paired = True
        if not paired:
            img_stars_x.append(ims[0]-0.5)
            img_stars_y.append(ims[1]-0.5)

    paired_img_stars_x = [p["image_star"][0]-0.5 for p in pairs]
    paired_img_stars_y = [p["image_star"][1]-0.5 for p in pairs]
    paired_cat_stars_x = [cs[0]-0.5 for cs in catalog_stars if mask is None or mask[min(1079, max(0, int(cs[1])))][min(1919, max(0, int(cs[0])))] != 0]
    paired_cat_stars_y = [cs[1]-0.5 for cs in catalog_stars if mask is None or mask[min(1079, max(0, int(cs[1])))][min(1919, max(0, int(cs[0])))] != 0]
    pyplot.scatter(paired_img_stars_x, paired_img_stars_y, c='magenta', marker='x', linewidths=1)
    pyplot.scatter(paired_cat_stars_x, paired_cat_stars_y, c='b', marker='+', linewidths=1)
    pyplot.scatter(img_stars_x, img_stars_y, c='r', marker='.', linewidths=1)

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

def open_image(image_file):
    im = cv2.imread(image_file)
    im = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    return im

def calib_iteration(image, catalog, calib, mask_file="mask.png", max_iter=5, initial_max_dist=15, max_dist_reduction=2):
    star_list = calib.get_star_list()
    image_stars = find_stars.detect(image, mask_file, star_list)
    calib.calibrate(np.array([(np.float64(ims[2]), np.float64(ims[3]), 0) for ims in image_stars]), np.array([(np.float64(ims[0]), np.float64(ims[1]), 0) for ims in image_stars]))

    cra = np.array([cs[0] for cs in catalog])
    cdec = np.array([cs[1] for cs in catalog])
    x, y = calib.ra_dec_to_xy(cra, cdec)
    catalog_stars = list(zip(x, y, cra, cdec))
    catalog_stars = [(s[0], s[1], s[2], s[3]) for s in catalog_stars if s[0] >= -10 and s[0] <= 1930 and s[1] >= -10 and s[1] <= 1090]
    image_stars = find_stars.detect(image, mask_file)
    
    max_dist = initial_max_dist
    last_rms_res = 1000
    i = 1
    print("Cat length: %s" % len(catalog))
    while True:
        cra = np.array([cs[0] for cs in catalog])
        cdec = np.array([cs[1] for cs in catalog])
        x, y = calib.ra_dec_to_xy(cra, cdec)
        catalog_stars = list(zip(x, y, cra, cdec))
        catalog_stars = [(s[0], s[1], s[2], s[3]) for s in catalog_stars if s[0] >= -10 and s[0] <= 1930 and s[1] >= -10 and s[1] <= 1090]

        ix = [ims[0] for ims in image_stars]
        iy = [ims[1] for ims in image_stars]
        ra_data, dec_data = calib.xy_to_ra_dec(ix, iy)

        image_stars = list(zip(ix, iy, ra_data, dec_data))
        pairs = pair_stars(catalog_stars, image_stars, max_dist_px=max_dist)

        if i == max_iter:
            return pairs, catalog_stars, image_stars

        rms_res = get_rms(pairs)
        print("RMS residual after iteration %s: %s arcmin" % (i, rms_res))
        if abs(rms_res - last_rms_res) < 0.005:
            print("No RMS residual change last iteration. Skipping further iterations.")
            return pairs, catalog_stars, image_stars
        last_rms_res = rms_res

        calib.set_star_pairs(pairs)
        print("Calibrating with %s star pairs" % len(pairs))
        calib.calibrate(iter=i)

        i += 1
        # With each iteration we want to be more discriminate with how close a pair has to be to be counted
        max_dist -= max_dist_reduction

def get_rms(pairs):
    if len(pairs) == 0:
        return 0
    rms_sum = 0
    for p in pairs:
        rms_sum += (p["residual_angle"]*60)**2
    rms = math.sqrt(rms_sum/len(pairs))
    return rms

def calc_residuals(image_file, catalog, calib, mask_file="mask.png"):
    image = open_image(image_file)
    pairs, catalog_stars, image_stars = calib_iteration(image, catalog, calib, mask_file, max_iter=5, initial_max_dist=15, max_dist_reduction=2)

    print("\n== RMS: %s arcmin" % get_rms(pairs))
    plot_residuals(pairs, 1920, 1080)
    show_residuals(image, pairs, catalog_stars, image_stars, mask_file)

def read_star_catalog(catalog_file, pos):
    catalog, _, _ = StarCatalog.readStarCatalog(".", catalog_file, lim_mag=4.8)
    planets = []
    for body in [ephem.Mercury(), ephem.Venus(), ephem.Mars(), ephem.Jupiter(), ephem.Saturn()]:
        body.compute(pos)
        planets.append((body.ra/math.pi*180, body.dec/math.pi*180, body.mag))
    catalog = np.append(catalog, planets, axis=0)
    return catalog
    

if __name__ == "__main__":
    calibration = importlib.import_module("calibration_methods.rms")
    calib = calibration.Calibration("parameters.json")
    catalog = read_star_catalog("BSC5.cat")
    print("Found %s catalog stars" % len(catalog))
    calc_residuals(None, "image.png", catalog, None, calib)
