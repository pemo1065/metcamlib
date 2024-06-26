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
from rich.console import Console
from rich.table import Table
from wand.image import Image
from wand.drawing import Drawing
from wand.color import Color
from utils import ra_dec_to_alt_az


def residual_pixels(x1, y1, x2, y2):
    return math.sqrt((x1-x2)**2 + (y1-y2)**2)

def residual_angle(d1, d2, ra1, ra2):
    return math.acos(math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(ra1-ra2))

def deg_to_rad(deg):
    return deg/360*2*math.pi

catalog_stars = []

def pair_stars(catalog_stars, image_stars, calib, max_dist_px=10):
    pairings = []
    print("Image stars: %s" % len(image_stars))
    print("catalog stars: %s" % len(catalog_stars))
    lat, lon, elev, timestamp = calib.get_pos()
    no_matches = []
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
        if min_dist <= max_dist_px:
            ang_sep = degrees(angular_separation(radians(best_match[2]), radians(best_match[3]), radians(ims[2]), radians(ims[3])))
            #if ang_sep > 1:
            #    continue
            ialt, iaz = ra_dec_to_alt_az(ims[2], ims[3], lat, lon, elev, timestamp)
            calt, caz = ra_dec_to_alt_az(best_match[2], best_match[3], lat, lon, elev, timestamp)

            dists = sorted(dists)
            #if dists[1]/dists[0] < 1.5:
            #    print("Ambiguous pairing for star at (%s, %s). Found two close matches with dists %s and %s. Skipping" % (ims[0], ims[1], dists[0], dists[1]))
            #    continue

            pairings.append({
                "image_star": [*ims, ialt, iaz],
                "catalog_star": [*best_match, calt, caz],
                "residual_px": min_dist,
                "residual_x": best_match[0] - ims[0],
                "residual_y": best_match[1] - ims[1],
                "residual_alt": (ialt - calt)*60,
                "residual_az": (iaz - caz)*60,
                "residual_angle": ang_sep
            })
        else:
            no_matches.append("(%s, %s) (dist: %s)" % (ims[0], ims[1], min_dist))

    if len(no_matches) > 0:
        print("Found no match for %s stars" % len(no_matches))
    return pairings

def res_table(pairs):
    console = Console()
    table = Table(show_header=True, header_style="bold magenta")
    table.add_column("Image star (x,y)", justify="right")
    table.add_column("Image star (ra,dec)", justify="right")
    table.add_column("Catalog star (x,y)", justify="right")
    table.add_column("Catalog star (ra,dec)", justify="right")
    table.add_column("Pixel distance", justify="right")
    table.add_column("Angular distance (arcmin)", justify="right")
    for pair in pairs:
        table.add_row(
            "(%s, %s)" % (format(pair["image_star"][0], ".2f"), format(pair["image_star"][1], ".2f")),
            "(%s, %s)" % (format(pair["image_star"][2], ".2f"), format(pair["image_star"][3], ".2f")),
            "(%s, %s)" % (format(pair["catalog_star"][0], ".2f"), format(pair["catalog_star"][1], ".2f")),
            "(%s, %s)" % (format(pair["catalog_star"][2], ".2f"), format(pair["catalog_star"][3], ".2f")),
            str(format(pair["residual_px"], ".2f")),
            str(format(60*pair["residual_angle"], ".2f")),
        )
    console.print(table)

def plot_residuals(pairs, im_h, im_w):
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = pyplot.subplots(2, 4, num="Residuals")
    fig.set_figwidth(10)
    fig.set_figheight(8)
    ax1.set_xlabel("Radius")
    ax1.set_ylabel("Residual (arcmin)")
    ax2.set_xlabel("X residual (px)")
    ax2.set_ylabel("Y residual (px)")
    ax3.set_xlabel("Pixel residual")
    ax3.set_ylabel("Angle residual")
    ax4.set_xlabel("Altitude (deg)")
    ax4.set_ylabel("Altitude residual (arcmin)")
    ax5.set_xlabel("Azimuth (deg)")
    ax5.set_ylabel("Azimuth residual (arcmin)")
    ax6.set_xlabel("X coordinate")
    ax6.set_ylabel("X residual (px)")
    ax7.set_xlabel("Y coordinate")
    ax7.set_ylabel("Y residual (px)")
    center = (im_h/2+0.5, im_w/2+0.5)
    residual_angles = []
    radii = []
    res_x = []
    res_y = []
    res_alt = []
    res_az = []

    sum_x = 0
    sum_y = 0
    for pair in pairs:
        ims = pair["image_star"]
        radii.append(math.sqrt((ims[0] - center[0])**2 + (ims[1] - center[1])**2))
        residual_angles.append(pair["residual_angle"]*60)
        res_x.append(pair["residual_x"])
        res_y.append(pair["residual_y"])
        res_alt.append(pair["residual_alt"])
        res_az.append(pair["residual_az"])
        sum_x += pair["residual_x"]
        sum_y += pair["residual_y"]
    px_avg = (sum_x/len(pairs), sum_y/len(pairs))
    ax1.scatter(radii, residual_angles)
    ax2.scatter(res_x, res_y)
    ax2.scatter([px_avg[0]], [px_avg[1]], c="r", marker="x")
    ax3.scatter([p["residual_px"] for p in pairs], [p["residual_angle"] for p in pairs])
    ax4.scatter([p["image_star"][4] for p in pairs], res_alt)
    ax5.scatter([p["image_star"][5] for p in pairs], res_az)
    ax6.scatter([p["image_star"][0] for p in pairs], res_x)
    ax7.scatter([p["image_star"][1] for p in pairs], res_y)
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
    res_table(pairs)
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
    cat_stars_x = [cs[0]-0.5 for cs in catalog_stars if mask is None or mask[min(1079, max(0, int(cs[1])))][min(1919, max(0, int(cs[0])))] != 0]
    cat_stars_y = [cs[1]-0.5 for cs in catalog_stars if mask is None or mask[min(1079, max(0, int(cs[1])))][min(1919, max(0, int(cs[0])))] != 0]
    pyplot.scatter(paired_img_stars_x, paired_img_stars_y, c='magenta', marker='x', linewidths=1)
    pyplot.scatter(cat_stars_x, cat_stars_y, c='b', marker='+', linewidths=1)
    pyplot.scatter(img_stars_x, img_stars_y, c='r', marker='.', linewidths=1)

    pyplot.show()

def get_roundtrip_rms(calib):
    xs = []
    ys = []
    for x in range(0, 1921, 120):
        for y in range(0, 1081, 120):
            xs.append(x)
            ys.append(y)
    ras, decs = calib.xy_to_ra_dec(xs, ys)

    x2, y2 = calib.ra_dec_to_xy(ras, decs)
    dists_sum= 0
    for i, _ in enumerate(xs):
        dists_sum += abs(xs[i]-x2[i])**2 + abs(ys[i]-y2[i])**2

    return math.sqrt(dists_sum)/len(xs)

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

def plot_stars(image_file, pairs):
    pic = Image(filename=image_file)
    stars = pic.clone()

    with Drawing() as draw:
        draw.fill_color = Color('white')
        draw.fill_opacity=0.0
        draw.stroke_color = Color('white')
        for pair in pairs:
            img_star = pair["image_star"]
            cat_star = pair["catalog_star"]
            x_step, y_step = get_normalized_vector(img_star[0], img_star[1], cat_star[0], cat_star[1])
            x, y = cat_star[0], cat_star[1]
            ix, iy = img_star[0], img_star[1]
            draw.circle((x, y), (x+10, y))
            draw.line((x + 11*x_step, y + 11*y_step), (x + x_step*11 + pair["residual_px"]*11*x_step, y+ y_step*11 + pair["residual_px"]*11*y_step))
            # Crosshair
            draw.line((x+6, y), (x+9, y))
            draw.line((x-9, y), (x-6, y))
            draw.line((x, y+6), (x, y+9))
            draw.line((x, y-9), (x, y-6))

            draw.push()
            draw.stroke_color = Color('red')
            draw.fill_color = Color('red')
            draw.circle((ix, iy), (ix+1, iy))
            draw.pop()
        draw(stars)
        stars.save(filename='out.png')

# Finds the minimuim and maximum ra/dec coordinates in the image
def get_max_min_ra_dec(image, calib):
    height, width = image.shape
    xs = []
    ys = []
    for x in range(0, width+1, 20):
        for y in range(0, height+1, 20):
            xs.append(x)
            ys.append(y)
    ras, decs = calib.xy_to_ra_dec(xs, ys)
    max_ra = 0
    max_dec = -90
    min_ra = 360
    min_dec = 90
    for ra in ras:
        if ra > max_ra:
            max_ra = ra
        if ra < min_ra:
            min_ra = ra
    for dec in decs:
        if dec > max_dec:
            max_dec = dec
        if dec < min_dec:
            min_dec = dec
    return max_ra, max_dec, min_ra, min_dec

def calib_iteration(image, catalog, calib, mask_file="mask.png", detection_threshold=14, max_stars=60):
    (max_iter, initial_max_dist, max_dist_reduction) = calib.suggested_params()
    star_list = calib.get_star_list()
    image_stars, _ = find_stars.detect(image, mask_file, star_list)

    print("Running first calibraton with %s stars" % len(image_stars))
    calib.calibrate(np.array([(np.float64(ims[2]), np.float64(ims[3]), 0) for ims in image_stars]), np.array([(np.float64(ims[0]), np.float64(ims[1]), 0) for ims in image_stars]), curr_rms=100)

    cra = np.array([cs[0] for cs in catalog])
    cdec = np.array([cs[1] for cs in catalog])
    x, y = calib.ra_dec_to_xy(cra, cdec)
    catalog_stars = list(zip(x, y, cra, cdec))
    catalog_stars = [(s[0], s[1], s[2], s[3]) for s in catalog_stars if s[0] >= 0 and s[0] < 1920 and s[1] >= 0 and s[1] < 1080]
    pairs = pair_stars(catalog_stars, image_stars, calib)
    image_stars, gaussian_params = find_stars.detect(image, mask_file, [], detection_threshold, max_stars)
    
    max_dist = initial_max_dist
    last_rms_res = 1000
    i = 1
    print("Cat length: %s" % len(catalog))
    keep_going = True
    while True:
        cra = np.array([cs[0] for cs in catalog])
        cdec = np.array([cs[1] for cs in catalog])
        x, y = calib.ra_dec_to_xy(cra, cdec)
        catalog_stars = list(zip(x, y, cra, cdec))

        # Some calibrations seem to behave strangely outside the boundaries of the image, and ra/dec values far outside can still end up inside the image
        # Calculate the min and max celestial coordinates on the image boundaries so that we can filter out stars that fall outside of these values 
        max_ra, max_dec, min_ra, min_dec = get_max_min_ra_dec(image, calib)
        catalog_stars = [(s[0], s[1], s[2], s[3]) for s in catalog_stars if s[0] >= 0 and s[0] <= 1920 and s[1] >= 0 and s[1] <= 1080 and s[2] > min_ra and s[2] < max_ra and s[3] > min_dec and s[3] < max_dec]

        ix = [ims[0] for ims in image_stars]
        iy = [ims[1] for ims in image_stars]
        ra_data, dec_data = calib.xy_to_ra_dec(ix, iy)

        image_stars = list(zip(ix, iy, ra_data, dec_data))
        pairs = pair_stars(catalog_stars, image_stars, calib, max_dist_px=max_dist)

        if len(pairs) == 0:
            print("No stars could be paired, Calibration failed.")
            return pairs, catalog_stars, image_stars, gaussian_params

        if i == max_iter or not keep_going:
            return pairs, catalog_stars, image_stars, gaussian_params

        calib.set_star_pairs(pairs)

        ang_rms_res, px_rms_res = get_rms(pairs)
        print("RMS residual after iteration %s: %s arcmin, %s px" % (i, ang_rms_res, px_rms_res))
        if i == 5:
            return pairs, catalog_stars, image_stars, gaussian_params
        if abs(ang_rms_res - last_rms_res) < 0.005 and i > 3:
            print("No RMS residual change last iteration. Skipping further iterations.")
            return pairs, catalog_stars, image_stars, gaussian_params

        print("Running iteration %s with %s star pairs" % (i+1, len(pairs)))
        keep_going = calib.calibrate(iter=i, curr_rms=ang_rms_res, last_rms=last_rms_res)

        last_rms_res = ang_rms_res
        i += 1
        # With each iteration we want to be more discriminate with how close a pair has to be to be counted
        max_dist -= max_dist_reduction

def get_rms(pairs):
    if len(pairs) == 0:
        return 0, 0
    ang_rms_sum = 0
    px_rms_sum = 0
    for p in pairs:
        ang_rms_sum += (p["residual_angle"]*60)**2
        px_rms_sum += p["residual_px"]**2
    ang_rms = math.sqrt(ang_rms_sum/len(pairs))
    px_rms = math.sqrt(px_rms_sum/len(pairs))
    return ang_rms, px_rms

def calc_residuals(image_file, catalog, calib, mask_file="mask.png", detection_threshold=14, max_stars=60):
    image = open_image(image_file)
    pairs, catalog_stars, image_stars, gaussian_params = calib_iteration(image, catalog, calib, mask_file, detection_threshold, max_stars)
    if len(pairs) == 0:
        return

    plot_stars(image_file, pairs)

    ang_rms_res, px_rms_res = get_rms(pairs)
    if ang_rms_res < 8:
        calib.write_initial_calib_file(pairs)
        print("Wrote initial calib file")
    rt_rms = get_roundtrip_rms(calib)
    print("\n== Residual RMS: %s arcmin, %s px. Roundtrip RMS: %s" % (ang_rms_res, px_rms_res, rt_rms))
    plot_residuals(pairs, 1920, 1080)
    show_residuals(image, pairs, catalog_stars, image_stars, mask_file)

    print("\nOutput written to %s" % calib.output_file())

def read_star_catalog(catalog_file, pos, lim_mag=5.0):
    print("Detecting stars with limit magnitude %s" % lim_mag)
    catalog, _, _ = StarCatalog.readStarCatalog(".", catalog_file, lim_mag=lim_mag)
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
