import argparse
import cv2
import datetime
import ephem
import find_stars
import importlib
import json
import os
import shutil
import subprocess
from astropy.time import Time
from calc_residuals import calc_residuals, read_star_catalog
from matplotlib import pyplot
from RMS.Formats import Platepar, StarCatalog
from pathlib import Path

OUTPUT_DIR = "output"

def get_observer(lat, lon, elev, timestamp):
  pos = ephem.Observer()
  pos.lat = str(lat)
  pos.lon = str(lon)
  pos.elevation = float(elev)
  pos.temp = 10.0
  
  dt = datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
  pos.date = dt.strftime('%Y/%m/%d %H:%M:%S')
  return pos


def copy_file(image_file, base_filename):
  Path(OUTPUT_DIR).mkdir(exist_ok=True)
  input_path = "%s/%s/%s.png" % (os.getcwd(), OUTPUT_DIR, base_filename)
  shutil.copyfile(image_file, input_path)
  return input_path

def initial_calibration(input_path):
  os.environ["AMS_HOME"] = os.getcwd() + "/ams"
  os.chdir("amscams/pipeline")
  print("Performing initial calibration. This may take a while.")
  with open('../../initial_calib.log', 'w') as out:
      subprocess.Popen(['python', 'Process.py', 'ac', input_path], stdout=out, stderr=out).wait()
  print("Initial calibration done.")
  os.chdir("../..")

def calibrate(image_file, mask_file, timestamp, calibration_method, calib_file=None, lim_mag=5.0, detection_threshold=14, max_stars=60):
  dt = datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
  base_filename = dt.strftime('%Y_%m_%d_%H_%M_%S_000_1')
  input_path = copy_file(image_file, base_filename)
  if calib_file is None:
    initial_calibration(input_path)
    calib_file = "ams/cal/freecal/%s/%s-stacked-calparams.json" % (base_filename, base_filename)

  calibration_method = importlib.import_module("calibration_methods.%s" % calibration_method)
  calib = calibration_method.Calibration(input_path, calib_file, timestamp)

  lat, lon, elev, timestamp = calib.get_pos()
  pos = get_observer(lat, lon, elev, timestamp)
  catalog = read_star_catalog("BSC5.cat", pos, lim_mag=lim_mag)

  print("Performing main calibration.")
  calc_residuals(input_path, catalog, calib, mask_file, detection_threshold, max_stars)

def open_image(image_file):
    im = cv2.imread(image_file)
    im = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    return im

def detect_stars(image_file, mask_file, detection_threshold, max_stars):
  image = open_image(image_file)
  stars, _ = find_stars.detect(image, mask_file, [], detection_threshold, max_stars)
  image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
  pyplot.imshow(image)
  pyplot.subplots_adjust(left=0.05, right=0.995, top=0.995, bottom=0.05)
  img_stars_x = []
  img_stars_y = []
  for ims in stars:
    img_stars_x.append(ims[0]-0.5)
    img_stars_y.append(ims[1]-0.5)
  pyplot.scatter(img_stars_x, img_stars_y, c='r', marker='.', linewidths=1)
  pyplot.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform a blind calibration on sky image')

    parser.add_argument('-i', '--image', dest='image', help='Input image', type=str, required=True)
    parser.add_argument('-m', '--mask', dest='mask', help='Mask image', type=str, required=False)
    parser.add_argument('-t', '--timestamp', dest='timestamp', help='UTC timestamp of capture in format "YYYY-MM-DD hh:mm:ss"', type=str, required=False, default=None)
    parser.add_argument('-c', '--calibration-method', dest='calibration_method', help='Calibration method to use', type=str, required=False, default=None)
    parser.add_argument('-l', '--lim-mag', dest='lim_mag', help='Magnitude limit for catalog stars', type=float, default=5.0, required=False)
    parser.add_argument('-f', '--initial-calib-file', dest='calib_file', help='Initial calparams file from AMS', type=str, default=None, required=False)
    parser.add_argument('-s', '--star-detection-threshold', dest='detection_threshold', help='Star detection threshold. 10-30 are reasonable values. Lower means increased sensitivity.', type=int, default=14, required=False)
    parser.add_argument('-d', '--max-stars', dest='max_stars', help='Maximum stars to return from star detector', type=int, default=60, required=False)
    parser.add_argument('function', help='Function to invoke. Possible values are "calibrate" (default) or "find_stars"', type=str, default="calibrate")

    args = parser.parse_args()

    if args.function == 'calibrate':
      if args.timestamp is None:
         print("'--timestamp' flag ('-t') is required")
         exit()
      if args.calibration_method is None:
         print("'--calibration-method' flag ('-c') is required")
         exit()

      calibrate(args.image, args.mask, args.timestamp, args.calibration_method, args.calib_file, args.lim_mag, args.detection_threshold, args.max_stars)
    elif args.function == 'find_stars':
       detect_stars(args.image, args.mask, args.detection_threshold, args.max_stars)

