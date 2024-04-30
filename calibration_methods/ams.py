import datetime
import json
import math
import numpy as np
import os
import subprocess
from astropy.coordinates import angular_separation
from RMS.Astrometry.Conversions import datetime2JD

class Calibration:
    def __init__(self, image_file, calparams_file, timestamp):
        self.image_file = image_file
        self.timestamp = timestamp
        self.calparams_file = calparams_file
        self.dt = datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        self.star_pairs = []
        with open(calparams_file) as cal_file:
            cal_json = json.load(cal_file)
            self.calparams = cal_json
            freecal_idx = {
                os.path.abspath(calparams_file): {
                    "cam_id": "1",
                    "total_res_px": 34.38162963092674,
                    "total_stars": len(cal_json["cat_image_stars"]),
                    "position_angle": cal_json["position_angle"],
                    "center_az": cal_json["center_az"],
                    "center_el": cal_json["center_el"],
                    "pixscale": cal_json["pixscale"],
                    "cal_date": 1
                }
            }
            with open("ams/cal/freecal_index.json", "w") as freecal_idx_json:
                json.dump(freecal_idx, freecal_idx_json, indent=4)

    def calibrate(self, c=[], i=[], iter=0):
        if len(c) > 0 and len(i) > 0:
            pairs = [{
                "image_star": (i[idx][0], i[idx][1]),
                "catalog_star": (0, 0, c[idx][0], c[idx][1])
            } for idx in range(len(c))]
            if iter > 0:
                self.set_star_pairs(pairs)
        os.chdir("amscams/pipeline")
        #out = open(os.devnull, 'wb')
        print("Calibration iteration", iter)
        with open('ams.log', 'a') as out:
            subprocess.Popen(['python', 'Process.py', 'deep_init', '1'], stdout=out, stderr=out).wait()
        os.chdir("../..")
        with open("ams/cal/multi_poly-AMSXXX-1.info") as multi_poly_file:
            self.multi_poly = json.load(multi_poly_file)
        self.calparams["x_poly"] = self.multi_poly["x_poly"]
        self.calparams["y_poly"] = self.multi_poly["y_poly"]
        self.calparams["x_poly_fwd"] = self.multi_poly["x_poly_fwd"]
        self.calparams["y_poly_fwd"] = self.multi_poly["y_poly_fwd"]
        self.calparams["ra_center"] = self.multi_poly["ra_center"]
        self.calparams["dec_center"] = self.multi_poly["dec_center"]
        self.calparams["position_angle"] = self.multi_poly["position_angle"]
        return True

    def get_pos(self):
        return float(self.calparams["device_lat"]), float(self.calparams["device_lng"]), float(self.calparams["device_alt"]), self.timestamp

    def add_additional_stars(self, pairs):
        if self.calparams["short_bright_stars"] == None:
            self.calparams["short_bright_stars"] = []
        for pair in pairs:
            #for star in self.calparams["cat_image_stars"]:
            #    name, mag, ra, dec, img_ra, img_dec, match_dist, new_x, new_y, img_az, img_el, new_cat_x, new_cat_y, imx, imy, real_res_px, star_flux = star
            #    if abs(ra -pair["catalog_star"][2]) < 0.05 and abs(dec - pair["catalog_star"][3]) < 0.05:
            #print("Adjusting (%s, %s) to (%s, %s)" % (imx, imy, pair["image_star"][0], pair["image_star"][1]))
            found = False
            for s in self.calparams["cat_image_stars"]:
                px = pair["image_star"][0]
                py = pair["image_star"][1]
                if abs(s[13]-px) < 1 and abs(s[14]-py) < 1:
                    found = True
                    break
            if not found:
                self.calparams["cat_image_stars"].append([
                    "",
                    1.0,
                    pair["catalog_star"][2],
                    pair["catalog_star"][3],
                    pair["image_star"][2],
                    pair["image_star"][3],
                    0.0,
                    pair["image_star"][0],
                    pair["image_star"][1],
                    0.0,
                    0.0,
                    pair["catalog_star"][0],
                    pair["catalog_star"][1],
                    pair["image_star"][0],
                    pair["image_star"][1],
                    0.0,
                    0.0,
                ])
                #name, name, ra, dec, mag, new_cat_x, new_cat_y, zp_cat_x, zp_cat_y, rx1,ry1,rx2,ry2
                self.calparams["short_bright_stars"].append([
                    "",
                    "",
                    pair["catalog_star"][2],
                    pair["catalog_star"][3],
                    1.0,
                    pair["catalog_star"][0],
                    pair["catalog_star"][1],
                    pair["catalog_star"][0],
                    pair["catalog_star"][1],
                    pair["catalog_star"][0] - 20,
                    pair["catalog_star"][1] - 20,
                    pair["catalog_star"][0] + 20,
                    pair["catalog_star"][1] + 20,
                ])
        with open(self.calparams_file, "w") as calparams:
            json.dump(self.calparams, calparams, indent=4)
            print("WROTE NEW STAR PAIRS")
    
    def set_star_pairs(self, pairs):
        new_cat_image_stars = []
        short_bright_stars = []
        img_x = np.array([p["image_star"][0] for p in pairs])
        img_y = np.array([p["image_star"][1] for p in pairs])
        img_ra, img_dec = self.xy_to_ra_dec(img_x, img_y)
        cat_ra = np.array([p["catalog_star"][0] for p in pairs])
        cat_dec = np.array([p["catalog_star"][1] for p in pairs])
        cat_x, cat_y = self.ra_dec_to_xy(cat_ra, cat_dec)
        for i in range(len(pairs)):
            ims = pairs[i]["image_star"]
            cs = pairs[i]["catalog_star"]
            pairs[i]["image_star"] = (ims[0], ims[1], img_ra[i], img_dec[i])
            pairs[i]["catalog_star"] = (cat_x[i], cat_y[i], cs[2], cs[3])
        
        self.add_additional_stars(pairs)

        self.star_pairs = pairs

    def get_star_list(self):
        return [(s[13], s[14], s[2], s[3]) for s in self.calparams["cat_image_stars"]]

    def ra_dec_to_xy(self, ra, dec):
        xs = []
        ys = []
        for r, d in zip(ra.tolist(), dec.tolist()):
            x, y = distort_xy(r, d, self.multi_poly["ra_center"], self.multi_poly["dec_center"],
                    self.multi_poly["x_poly"], self.multi_poly["y_poly"], self.multi_poly["imagew"],
                    self.multi_poly["imageh"], self.multi_poly["position_angle"],
                    3600/float(self.multi_poly['pixscale']))
            xs.append(x)
            ys.append(y)
        return xs, ys
    
    def xy_to_ra_dec(self, xs, ys):
        ras = []
        decs = []
        for x, y in zip(xs, ys):
            _, _, ra, dec, _, _ = XYtoRADec(x, y, self.image_file, self.calparams)
            ras.append(ra)
            decs.append(dec)
        return ras, decs
    

#####################################
### Functions from amscams repository
#####################################

def distort_xy(ra,dec,RA_center, dec_center, x_poly, y_poly, x_res, y_res, pos_angle_ref,F_scale=1):

   ra_star = np.float64(ra)
   dec_star = np.float64(dec)
   RA_center = np.float64(RA_center)
   dec_center = np.float64(dec_center)
   pos_angle_ref = np.float64(pos_angle_ref)
   F_scale = np.float64(F_scale)
   debug = 0 
   if debug == 1:
      print("DISTORT XY")
      print("STR RA/DEC:", ra_star, dec_star)
      print("CENTER RA/DEC:", RA_center, dec_center)
      print("XP:", x_poly, type(x_poly))
      print("YP:", y_poly, type(y_poly))
      print("XRES:", x_res)
      print("YRES:", y_res)
      print("POS:", pos_angle_ref)
      print("FSCALE:", F_scale)

   #F_scale = F_scale/10
   w_pix = 50*F_scale/3600
   #F_scale = 158 * 2
   #F_scale = 155
   #F_scale = 3600/16
   #F_scale = 3600/F_scale
   #F_scale = 1

   # Gnomonization of star coordinates to image coordinates
   ra1 = math.radians(float(RA_center))
   dec1 = math.radians(float(dec_center))
   ra2 = math.radians(float(ra_star))
   dec2 = math.radians(float(dec_star))
   ad = math.acos(math.sin(dec1)*math.sin(dec2) + math.cos(dec1)*math.cos(dec2)*math.cos(ra2 - ra1))
   radius = math.degrees(ad)
   
   try:
      sinA = math.cos(dec2)*math.sin(ra2 - ra1)/math.sin(ad)
      cosA = (math.sin(dec2) - math.sin(dec1)*math.cos(ad))/(math.cos(dec1)*math.sin(ad))
   except:
      sinA = 0
      cosA = 0
   theta = -math.degrees(math.atan2(sinA, cosA))
   theta = theta + pos_angle_ref - 90.0
   #theta = theta + pos_angle_ref - 90 + (1000*x_poly[12]) + (1000*y_poly[12])
   #theta = theta + pos_angle_ref - 90



   dist = np.degrees(math.acos(math.sin(dec1)*math.sin(dec2) + math.cos(dec1)*math.cos(dec2)*math.cos(ra1 - ra2)))

   # Calculate the image coordinates (scale the F_scale from CIF resolution)
   X1 = radius*math.cos(math.radians(theta))*F_scale
   Y1 = radius*math.sin(math.radians(theta))*F_scale
   # Calculate distortion in X direction
   dX = (x_poly[0]
      + x_poly[1]*X1
      + x_poly[2]*Y1
      + x_poly[3]*X1**2
      + x_poly[4]*X1*Y1
      + x_poly[5]*Y1**2
      + x_poly[6]*X1**3
      + x_poly[7]*X1**2*Y1
      + x_poly[8]*X1*Y1**2
      + x_poly[9]*Y1**3
      + x_poly[10]*X1*math.sqrt(X1**2 + Y1**2)
      + x_poly[11]*Y1*math.sqrt(X1**2 + Y1**2))

   # NEW TERMS DONT WORK WELL
   #dX += x_poly[12]*X1*math.sqrt(X1**2 + Y1**2)**3
   #dX += x_poly[13]*X1*math.sqrt(X1**2 + Y1**2)**5
   # Add the distortion correction and calculate X image coordinates
   #x_array[i] = (X1 - dX)*x_res/384.0 + x_res/2.0
   new_x = X1 - dX + x_res/2.0

   # Calculate distortion in Y direction
   dY = (y_poly[0]
      + y_poly[1]*X1
      + y_poly[2]*Y1
      + y_poly[3]*X1**2
      + y_poly[4]*X1*Y1
      + y_poly[5]*Y1**2
      + y_poly[6]*X1**3
      + y_poly[7]*X1**2*Y1
      + y_poly[8]*X1*Y1**2
      + y_poly[9]*Y1**3
      + y_poly[10]*Y1*math.sqrt(X1**2 + Y1**2)
      + y_poly[11]*X1*math.sqrt(X1**2 + Y1**2))

   # NEW TERMS DONT WORK WELL
   #dY += y_poly[12]*Y1*math.sqrt(X1**2 + Y1**2)**3
   #dY += y_poly[13]*Y1*math.sqrt(X1**2 + Y1**2)**5

   # Add the distortion correction and calculate Y image coordinates
   #y_array[i] = (Y1 - dY)*y_res/288.0 + y_res/2.0
   new_y = Y1 - dY + y_res/2.0
   return(new_x,new_y)

def convert_filename_to_date_cam(file, ms = 0):
   el = file.split("/")
   filename = el[-1]
   filename = filename.replace(".mp4" ,"")
   if "-" in filename:
      xxx = filename.split("-")
      filename = xxx[0]
   el = filename.split("_")
   if len(el) >= 8:
      fy,fm,fd,fh,fmin,fs,fms,cam = el[0], el[1], el[2], el[3], el[4], el[5], el[6], el[7]
   else:
      fy,fm,fd,fh,fmin,fs,fms,cam = "1999", "01", "01", "00", "00", "00", "000", "010001"
   if "-" in cam:
      cel = cam.split("-")
      cam = cel[0]

   #print("CAM:", cam)
   #exit()
   cam = cam.replace(".png", "")
   cam = cam.replace(".jpg", "")
   cam = cam.replace(".json", "")
   cam = cam.replace(".mp4", "")

   f_date_str = fy + "-" + fm + "-" + fd + " " + fh + ":" + fmin + ":" + fs
   f_datetime = datetime.datetime.strptime(f_date_str, "%Y-%m-%d %H:%M:%S")
   if ms == 1:
      return(f_datetime, cam, f_date_str,fy,fm,fd, fh, fmin, fs,fms)
   else:
      return(f_datetime, cam, f_date_str,fy,fm,fd, fh, fmin, fs)

def XYtoRADec(img_x,img_y,cal_file,cal_params):
   hd_datetime, hd_cam, hd_date, hd_y, hd_m, hd_d, hd_h, hd_M, hd_s = convert_filename_to_date_cam(cal_file)
   F_scale = 3600/float(cal_params['pixscale'])
   #F_scale = 24

   total_min = (int(hd_h) * 60) + int(hd_M)
   day_frac = total_min / 1440 
   hd_d = int(hd_d) + day_frac
   #jd = date_to_jd(int(hd_y),int(hd_m),float(hd_d))

   (f_datetime, cam, f_date_str,y,m,d, h, mm, s) = convert_filename_to_date_cam(cal_file)
   jd = datetime2JD(f_datetime, 0.0)
   #hour_angle = JD2HourAngle(jd)

   lat = float(cal_params['device_lat'])
   lon = float(cal_params['device_lng'])

   # Calculate the reference hour angle
   T = (jd - 2451545.0)/36525.0
   Ho = (280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*T**2 \
      - (T**3)/38710000.0)%360
   if "x_poly_fwd" in cal_params:
      x_poly_fwd = cal_params['x_poly_fwd']
      y_poly_fwd = cal_params['y_poly_fwd']
   else:
      x_poly_fwd = np.zeros(shape=(15,),dtype=np.float64)
      y_poly_fwd = np.zeros(shape=(15,),dtype=np.float64)
   
   dec_d = float(cal_params['dec_center']) 
   RA_d = float(cal_params['ra_center']) 

   dec_d = dec_d + (x_poly_fwd[13] * 100)
   dec_d = dec_d + (y_poly_fwd[13] * 100)

   RA_d = RA_d + (x_poly_fwd[14] * 100)
   RA_d = RA_d + (y_poly_fwd[14] * 100)

   pos_angle_ref = float(cal_params['position_angle']) + (1000*x_poly_fwd[12]) + (1000*y_poly_fwd[12])

   # Convert declination to radians
   dec_rad = math.radians(dec_d)

   # Precalculate some parameters
   sl = math.sin(math.radians(lat))
   cl = math.cos(math.radians(lat))

   if "imagew" not in cal_params:
      cal_params['imagew'] = 1920
      cal_params['imageh'] = 1080 

   x_det = img_x - int(cal_params['imagew'])/2
   y_det = img_y - int(cal_params['imageh'])/2

   #x = img_x
   #y = img_y
   #x0 = x_poly[0]
   #y0 = y_poly[0]

   #r = math.sqrt((x - x0)**2 + (y - y0)**2)

   dx = (x_poly_fwd[0]
      + x_poly_fwd[1]*x_det
      + x_poly_fwd[2]*y_det
      + x_poly_fwd[3]*x_det**2
      + x_poly_fwd[4]*x_det*y_det
      + x_poly_fwd[5]*y_det**2
      + x_poly_fwd[6]*x_det**3
      + x_poly_fwd[7]*x_det**2*y_det
      + x_poly_fwd[8]*x_det*y_det**2
      + x_poly_fwd[9]*y_det**3
      + x_poly_fwd[10]*x_det*math.sqrt(x_det**2 + y_det**2)
      + x_poly_fwd[11]*y_det*math.sqrt(x_det**2 + y_det**2))

   #dx += x_poly_fwd[12]*x_det*math.sqrt(x_det**2 + y_det**2)**3
   #dx += x_poly_fwd[13]*x_det*math.sqrt(x_det**2 + y_det**2)**5
   # Add the distortion correction
   x_pix = x_det + dx 

   dy = (y_poly_fwd[0]
      + y_poly_fwd[1]*x_det
      + y_poly_fwd[2]*y_det
      + y_poly_fwd[3]*x_det**2
      + y_poly_fwd[4]*x_det*y_det
      + y_poly_fwd[5]*y_det**2
      + y_poly_fwd[6]*x_det**3
      + y_poly_fwd[7]*x_det**2*y_det
      + y_poly_fwd[8]*x_det*y_det**2
      + y_poly_fwd[9]*y_det**3
      + y_poly_fwd[10]*y_det*math.sqrt(x_det**2 + y_det**2)
      + y_poly_fwd[11]*x_det*math.sqrt(x_det**2 + y_det**2))

   #dy += y_poly_fwd[12]*y_det*math.sqrt(x_det**2 + y_det**2)**3
   #dy += y_poly_fwd[13]*y_det*math.sqrt(x_det**2 + y_det**2)**5


   # Add the distortion correction
   y_pix = y_det + dy 

   x_pix = x_pix / F_scale
   y_pix = y_pix / F_scale

   ### Convert gnomonic X, Y to alt, az ###

   # Caulucate the needed parameters
   radius = math.radians(math.sqrt(x_pix**2 + y_pix**2))
   theta = math.radians((90 - pos_angle_ref + math.degrees(math.atan2(y_pix, x_pix)))%360)

   sin_t = math.sin(dec_rad)*math.cos(radius) + math.cos(dec_rad)*math.sin(radius)*math.cos(theta)
   Dec0det = math.atan2(sin_t, math.sqrt(1 - sin_t**2))

   sin_t = math.sin(theta)*math.sin(radius)/math.cos(Dec0det)
   cos_t = (math.cos(radius) - math.sin(Dec0det)*math.sin(dec_rad))/(math.cos(Dec0det)*math.cos(dec_rad))
   RA0det = (RA_d - math.degrees(math.atan2(sin_t, cos_t)))%360

   h = math.radians(Ho + lon - RA0det)
   sh = math.sin(h)
   sd = math.sin(Dec0det)
   ch = math.cos(h)
   cd = math.cos(Dec0det)

   x = -ch*cd*sl + sd*cl
   y = -sh*cd
   z = ch*cd*cl + sd*sl

   r = math.sqrt(x**2 + y**2)

   # Calculate azimuth and altitude
   azimuth = math.degrees(math.atan2(y, x))%360
   altitude = math.degrees(math.atan2(z, r))



   ### Convert alt, az to RA, Dec ###

   # Never allow the altitude to be exactly 90 deg due to numerical issues
   if altitude == 90:
      altitude = 89.9999

   # Convert altitude and azimuth to radians
   az_rad = math.radians(azimuth)
   alt_rad = math.radians(altitude)

   saz = math.sin(az_rad)
   salt = math.sin(alt_rad)
   caz = math.cos(az_rad)
   calt = math.cos(alt_rad)

   x = -saz*calt
   y = -caz*sl*calt + salt*cl
   HA = math.degrees(math.atan2(x, y))

   # Calculate the hour angle
   T = (jd - 2451545.0)/36525.0
   hour_angle = (280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*T**2 - T**3/38710000.0)%360

   RA = (hour_angle + lon - HA)%360
   dec = math.degrees(math.asin(sl*salt + cl*calt*caz))

   return(x_pix+img_x,y_pix+img_y,RA,dec,azimuth,altitude)
