import cv2 
import numpy as np 
import gaussian_fit

padding = 6

def find_max_px_coords(im, initial_guess=(padding, padding)):
    h, w = im.shape[:2]
    max_pixel = im.max()
    for y_coord in range(h):
        for x_coord in range(w):
            if im[y_coord, x_coord] == max_pixel:
                return (x_coord, y_coord)
    return initial_guess


def find_blobs(image):
    params = cv2.SimpleBlobDetector_Params()

    params.filterByArea = True
    params.maxArea = 60
    params.minArea = 8
    params.minThreshold = 5
    params.thresholdStep = 5

    detector = cv2.SimpleBlobDetector_create(params)

    return [(p.pt[0], p.pt[1]) for p in detector.detect(image)]


def detect(image, mask_file=None, star_list=[], min_px_diff=14):
    image_inv = cv2.bitwise_not(image)

    mask = None
    if mask_file is not None:
        mask = cv2.imread(mask_file, 0)

    # Don't use mask if the star list is supplied from an external source
    use_mask = False
    if star_list is None or len(star_list) == 0:
        star_list = find_blobs(image_inv)
        use_mask = True

    stars = []
    params = []
    # Remove points that are masked out and find gaussian fits for remaining points
    for point in star_list:
        x, y = point[0], point[1]
        if use_mask and mask is not None and mask[int(y)][int(x)] == 0:
            continue

        part = image[int(y)-padding:int(y)+padding, int(x)-padding:int(x)+padding]
        try:
            min = part.min()
            max = part.max()
            if max-min < min_px_diff:
                continue
        except Exception as e:
            print("Error: %s" % e)
    

        try:
            guess = find_max_px_coords(part)

            amp, x0, y0, sigma_x, sigma_y, theta = gaussian_fit.best_fit(part, guess)
            params.append((part, f'{(int(x) - padding + x0):.2f}, {(int(y) - padding + y0):.2f}', (amp, x0, y0, sigma_x, sigma_y, theta)))

            star = (int(x) - padding + x0, int(y) - padding + y0, point[2] if len(point) > 3 else 0, point[3] if len(point) > 3 else 0)
            stars.append(star)

        except Exception as e:
            print("Failed to find fit for star at [%s, %s]: %s" % (x, y, e))

    #gaussian_fit.plot_contours(params)

    return stars, params


if __name__ == "__main__":
    detect("image.png", "mask.png")
