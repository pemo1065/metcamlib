import datetime
import ephem
import hsi
import json
import math
import subprocess
from astropy.coordinates import angular_separation
from calibration_methods.stars import cat
from utils import ra_dec_to_alt_az, alt_az_to_ra_dec

MAX_MATCH_DIST=10

class Calibration:
    def __init__(self, image_file, ams_file, timestamp):
        self.first = True
        self.star_pairs = []
        self.timestamp = timestamp
        dt = datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        self.lens_file = "output/%s" % dt.strftime("%Y_%m_%d_%H_%M_%S.pto")
        self.image_file = image_file
        with open(ams_file) as f:
            self.data = json.load(f)

            self.width = self.data['imagew']
            self.height = self.data['imageh']

            self.center_az = self.data['center_az']
            self.center_el = self.data['center_el']
                            
            self.pos = ephem.Observer()
            self.pos.lat = self.data["device_lat"]
            self.pos.lon = self.data["device_lng"]
            self.pos.elevation = float(self.data["device_alt"])
            self.pos.temp = 10.0

            self.pos.date = dt.strftime('%Y/%m/%d %H:%M:%S')

            if 'pixel_scale' in self.data:
                self.pixel_scale = float(self.data['pixel_scale'])
            else:
                self.pixel_scale = float(self.data['pixscale'])

    def write_initial_calib_file(self, pairs):
        alt_c, az_c = self.xy_to_alt_az(self.width/2.0 - 0.5, self.height/2.0 - 0.5)
        alt_l, az_l = self.xy_to_alt_az(0, self.height/2.0 - 0.5)
        alt_r, az_r = self.xy_to_alt_az(self.width-0.5, self.height/2.0 - 0.5)
        alt_t, az_t = self.xy_to_alt_az(self.width/2.0, 0)
        alt_b, az_b = self.xy_to_alt_az(self.width/2.0, self.height - 0.5)
        print("Az/alt (l): %s/%s, Az/alt (r): %s%s" % (az_l, alt_l, az_r, alt_r))
        fov = math.degrees(angular_separation(math.radians(alt_l), math.radians(az_l), math.radians(alt_r), math.radians(az_r)))
        vfov = math.degrees(angular_separation(math.radians(alt_t), math.radians(az_t), math.radians(alt_b), math.radians(az_b)))
        print("FOV: %s" % fov)
        print("VFOV: %s" % vfov)
        ra, dec = alt_az_to_ra_dec(alt_c, az_c, self.pos.lat, self.pos.lon, self.pos.elev, self.timestamp)
        calparams = {
            "imagew": self.width,
            "imageh": self.height,
            "device_lat": math.degrees(self.pos.lat),
            "device_lng": math.degrees(self.pos.lon),
            "device_alt": self.pos.elevation,
            "center_az": az_c,
            "center_el": alt_c,
            "dec_center": dec,
            "ra_center": ra,
            "pixel_scale": self.pixel_scale,
            "fieldw": fov,
            "fieldh": vfov,
            "cat_image_stars": [],
            "short_bright_stars": [],
            "orig_pos_ang": self.data["position_angle"],
        }
        for pair in pairs:
            calparams["cat_image_stars"].append([
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
            calparams["short_bright_stars"].append([
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
        with open("output/calparams.json", "w") as cp:
            json.dump(calparams, cp)

    
    def write_pto_file(self):
        with open(self.lens_file, 'w') as output:
            print('''# hugin project file
#hugin_ptoversion 2
p f2 w36000 h18000 v360  E0 R0 n"TIFF_m c:LZW"
m g1 i0 m2 p0.00784314

# image lines
#-hugin  cropFactor=1
i w''' + str(self.width) + ' h' + str(self.height) + ' f3 v' + str(self.width * self.pixel_scale / 3600) + ' Ra0 Rb0 Rc0 Rd0 Re0 Eev0 Er1 Eb1 r0 p' + str(self.center_el) + ' y' + str(self.center_az - 180) + ''' TrX0 TrY0 TrZ0 Tpy0 Tpp0 j0 a0 b0 c0 d0 e0 g0 t0 Va1 Vb0 Vc0 Vd0 Vx0 Vy0  Vm5
i w36000 h18000 f4 v360 Ra0 Rb0 Rc0 Rd0 Re0 Eev0 Er1 Eb1 r0 p0 y0 TrX0 TrY0 TrZ0 j0 a0 b0 c0 d0 e0 g0 t0 Va1 Vb0 Vc0 Vd0 Vx0 Vy0  Vm5 n"dummy.jpg"


# specify variables that should be optimized
v v0
v r0
v p0
v y0
v a0
v b0
v c0
v d0
v e0
v

''' + ('# ' + str(self.pos.date)) + '''
# control points''', file=output)

            if len(self.star_pairs) == 0:
                for star in self.data['cat_image_stars']:
                    _,mag,ra,dec,_,_,match_dist,_,_,_,_,_,_,six,siy,_,_ = star
                    ra = ra * 24 / 360

                    if match_dist > MAX_MATCH_DIST:
                        continue

                    # Use AMS ra/dec to find the star in the NMN catalogue
                    min = 99999
                    for (ra2, pmra, dec2, pmdec, mag2, name) in cat:
                        body1 = ephem.FixedBody()
                        body1._ra, body1._pmra, body1._dec, body1._pmdec, body1._epoch = str(ra), pmra, str(dec), pmdec, ephem.J2000
                        body1.mag = mag
                        body1.compute(self.pos)
                        body2 = ephem.FixedBody()
                        body2._ra, body2._pmra, body2._dec, body2._pmdec, body2._epoch = str(ra2), pmra, str(dec2), pmdec, ephem.J2000
                        body2.mag = mag
                        body2.compute(self.pos)
                        separation = float(repr(ephem.separation(body1, body2)))
                        if (separation < min):
                            if abs(mag - mag2) > 0.3:  # Quick test whether same star
                                continue
                            min = separation
                            best = name
                            bestbody = body2
                    if min < 0.0001:
                        bestbody.compute(self.pos)
                        az = math.degrees(float(repr(bestbody.az)))
                        alt = math.degrees(float(repr(bestbody.alt)))
                        self.star_pairs.append({"image_star": (six, siy)})
                        if alt > 1:
                            print('c n0 N1 x' + str(six) + ' y' + str(siy) + ' X' + str(az*100) + ' Y' + str((90-alt)*100) + ' t0  # ' + best, file=output)
                        else:
                            print("Skipping point because alt is %s" % alt)
                    else:
                        print("Skipping star %s because separation is %s" % (star, min))
            else:
                for pair in self.star_pairs:
                    x, y = pair["image_star"][0:2]
                    ra, dec = pair["catalog_star"][2:4]
                    alt, az = ra_dec_to_alt_az(ra, dec, self.pos.lat, self.pos.lon, self.pos.elev, self.timestamp)
                    print('c n0 N1 x' + str(x) + ' y' + str(y) + ' X' + str(az*100) + ' Y' + str((90-alt)*100) + ' t0 ', file=output)
    
    def get_pos(self):
        return self.pos.lat/math.pi*180, self.pos.lon/math.pi*180, self.pos.elev, self.timestamp

    def set_star_pairs(self, pairs):
        self.star_pairs = pairs

    def output_file(self):
        return self.lens_file

    # Returns the suggested iterations, and the reduction in fitting distance with each iteration
    def suggested_params(self):
        return (6, 14, 2)

    def calibrate(self, a=None, b=None, iter=None, curr_rms=None, last_rms=None):
        self.write_pto_file()
        self.pano = hsi.Panorama()
        self.pano.ReadPTOFile(self.lens_file)
        self.pano.setOptimizeVector([('r', 'p', 'y', 'v', 'a', 'b', 'c', 'd', 'e'), ()])
        self.pano.WritePTOFile(self.lens_file)
        with open("pano.log", "a") as out:
            #subprocess.Popen(['cpclean', '-n', '1', '-o', "lens.pto", "lens.pto"], stdout=out, stderr=out).wait()
            subprocess.Popen(['autooptimiser', '-n', self.lens_file, '-o', self.lens_file], stdout=out, stderr=out).wait()
        self.pano.ReadPTOFile(self.lens_file)

        img = self.pano.getImage(0)
        self.tf = hsi.Transform()
        self.tf.createTransform(img, self.pano.getOptions())
        self.itf = hsi.Transform()
        self.itf.createInvTransform(img, self.pano.getOptions())
        return last_rms is None or abs(last_rms - curr_rms) > 0.01

    def get_star_list(self):
        return [(s["image_star"][0], s["image_star"][1]) for s in self.star_pairs]

    def xy_to_alt_az(self, x, y):
        dst = hsi.FDiff2D()
        self.itf.transformImgCoord(dst, hsi.FDiff2D(x, y))
        alt = 90 - dst.y/100
        az = (dst.x/100)%360
        return alt, az

    def alt_az_to_xy(self, alt, az):
        scale = int(self.pano.getOptions().getWidth() / self.pano.getOptions().getHFOV())
        dst = hsi.FDiff2D()
        self.tf.transformImgCoord(dst, hsi.FDiff2D(float(az*scale), float((90-alt)*scale)))
        return dst.x, dst.y

    def ra_dec_to_xy(self, ra, dec):
        scale = int(self.pano.getOptions().getWidth() / self.pano.getOptions().getHFOV())
        xs = []
        ys = []
        for r, d in zip(ra, dec): 
            alt, az = ra_dec_to_alt_az(r, d, self.pos.lat, self.pos.lon, self.pos.elevation, self.timestamp)
            dst = hsi.FDiff2D()
            self.tf.transformImgCoord(dst, hsi.FDiff2D(float(az*scale), float((90-alt)*scale)))
            xs.append(dst.x)
            ys.append(dst.y)
        
        return xs, ys
    
    def xy_to_ra_dec(self, xs, ys):
        ras = []
        decs = []
        for x, y in zip(xs, ys):
            dst = hsi.FDiff2D()
            self.itf.transformImgCoord(dst, hsi.FDiff2D(x, y))
            alt = 90 - dst.y/100
            az = (dst.x/100)%360
            ra, dec = alt_az_to_ra_dec(alt, az, self.pos.lat, self.pos.lon, self.pos.elevation, self.timestamp)
            ras.append(ra)
            decs.append(dec)
        return ras, decs
