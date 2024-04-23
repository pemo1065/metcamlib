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


def detect(image_file, mask_file=None, star_list=[]):
    image = cv2.imread(image_file, 0) 
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
            if max-min < 20:
                print("Skipping blob at %s, %s because diff is %s" % (x, y, max-min))
                continue
            else:
                print("Including blob at %s, %s because diff is %s" % (x, y, max-min))
        except Exception as e:
            print("Error: %s" % e)
    

        try:
            #min_pixel = part.min()
            #part1 = part - min_pixel
            #_, part1 = cv2.threshold(part1, 10, 255, cv2.THRESH_TOZERO)

            guess = find_max_px_coords(part)

            amp, x0, y0, sigma_x, sigma_y, theta = gaussian_fit.best_fit(part, guess)
            params.append((part, f'{(int(x) - padding + x0):.2f}, {(int(y) - padding + y0):.2f}', (amp, x0, y0, sigma_x, sigma_y, theta)))
            #if abs(int(x) - padding + x1 - 1245) < 5 and abs(int(y) - padding + y1 - 765) < 5:
            #    print("Found center at (%s, %s)" % (int(x) - padding + x1, int(y) - padding + y1))
            #    print("Fit guess (%s, %s) to (%s, %s)" % (guess[0], guess[1], int(x) - padding + x1 - 1245, int(y) - padding + y1 - 765))
            #    cv2.namedWindow("st_detail", cv2.WINDOW_NORMAL) 
            #    cv2.imshow("st_detail", np.concatenate((part, part1), axis=1) )
            #    cv2.waitKey(0)

            star = (int(x) - padding + x0, int(y) - padding + y0, point[2] if len(point) > 3 else 0, point[3] if len(point) > 3 else 0)
            stars.append(star)

        except Exception as e:
            print("Failed to find fit for star at [%s, %s]: %s" % (x, y, e))

    gaussian_fit.plot_contours(params)

    return stars


if __name__ == "__main__":
    detect("image.png", "mask.png")
