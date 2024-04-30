import argparse
import datetime
import ephem
import importlib
import json
import os
import shutil
from astropy.time import Time
from calc_residuals import calc_residuals, read_star_catalog
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
  print("OBS: %s" % pos)
  return pos


def copy_file(image_file, base_filename):
  Path(OUTPUT_DIR).mkdir(exist_ok=True)
  input_path = "%s/%s/%s.png" % (os.getcwd(), OUTPUT_DIR, base_filename)
  shutil.copyfile(image_file, input_path)
  return input_path


#image_file = "/home/petter/2024_03_07_22_05_42_000_1.png"
#image_file = "/home/petter/Downloads/data/high/2024-03-07/UPP_CAM4_2024-03-07-22-05-42_avg.png"
#image_file = "/home/petter/Downloads/2023_10_08_00_00_00_000_7_avg.png"
#date = "2024-03-07-22-05-42"


def initial_calibration(input_path):
  os.environ["AMS_HOME"] = os.getcwd() + "/ams"
  os.chdir("amscams/pipeline")
  os.system("python Process.py ac %s" % input_path)
  os.chdir("../..")

#calib = calibration.Calibration(input_path, "ams/cal/freecal/2024_03_07_22_05_42_000_1/2024_03_07_22_05_42_000_1-stacked-calparams.json", timestamp.to_string())
#calib = calibration.Calibration(input_path, "ams/cal/freecal/2023_10_08_00_00_00_000_7_avg/2023_10_08_00_00_00_000_7_avg-stacked-calparams.json", timestamp.to_string())

def calibrate(image_file, mask_file, timestamp, calibration_method):
  dt = datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
  base_filename = dt.strftime('%Y_%m_%d_%H_%M_%S_000_1')
  input_path = copy_file(image_file, base_filename)
  initial_calibration(input_path)

  calib_file = "ams/cal/freecal/%s/%s-stacked-calparams.json" % (base_filename, base_filename)

  calibration_method = importlib.import_module("calibration_methods.%s" % calibration_method)
  calib = calibration_method.Calibration(input_path, calib_file, timestamp)

  lat, lon, elev, timestamp = calib.get_pos()
  pos = get_observer(lat, lon, elev, timestamp)
  catalog = read_star_catalog("BSC5.cat", pos)

  calc_residuals(input_path, catalog, calib, mask_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform a blind calibration on sky image')

    parser.add_argument('-i', '--image', dest='image', help='Input image', type=str, required=True)
    parser.add_argument('-m', '--mask', dest='mask', help='Mask image', type=str, required=False)
    parser.add_argument('-t', '--timestamp', dest='timestamp', help='UTC timestamp of capture in format "YYYY-MM-DD HH:mm:ss"', type=str, required=True)
    parser.add_argument('-c', '--calibration-method', dest='calibration_method', help='Calibration method to use', type=str, required=True)

    args = parser.parse_args()
    calibrate(args.image, args.mask, args.timestamp, args.calibration_method)
