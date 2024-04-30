import datetime
import json
import numpy as np
from astropy.time import Time
from RMS.Astrometry.ApplyAstrometry import xyToRaDecPP, raDecToXYPP
from RMS.Formats import Platepar

class Calibration:
    def __init__(self, image_file, calparams_file, timestamp):
        self.star_pairs = []
        platepar_file = self.calparams_to_platepar(calparams_file, timestamp)
        self.platepar: Platepar = Platepar.Platepar()
        self.platepar.read(platepar_file, "json")
        self.timestamp = timestamp

    
    def calparams_to_platepar(self, calparams_file, timestamp):
        with open(calparams_file) as cp:
            calparams = json.load(cp)
        dt = datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        timestamp = Time(dt.strftime('%Y-%m-%dT%H:%M:%S.000Z'), format="isot", scale="utc")
        star_list = [[timestamp.jd, s[13], s[14], 10000.0, s[2], s[3], s[1]] for s in calparams["cat_image_stars"]]
        platepar = {
            "F_scale": 1920.0/float(calparams["fieldw"]),
            "Ho": 0,
            "JD": timestamp.jd,
            "RA_d": calparams["ra_center"],
            "UT_corr": 0,
            "X_res": calparams["imagew"],
            "Y_res": calparams["imageh"],
            "alt_centre": calparams["center_el"],
            "asymmetry_corr": False,
            "auto_check_fit_refined": False,
            "auto_recalibrated": False,
            "az_centre": calparams["center_az"],
            "dec_d": calparams["dec_center"],
            "distortion_type": "radial7-odd",
            "distortion_type_list": [
                "poly3+radial",
                "poly3+radial3",
                "poly3+radial5",
                "radial3-all",
                "radial4-all",
                "radial5-all",
                "radial3-odd",
                "radial5-odd",
                "radial7-odd",
                "radial9-odd"
            ],
            "distortion_type_poly_length": [
                12,
                13,
                14,
                7,
                8,
                9,
                6,
                7,
                8,
                9
            ],
            "x_poly": calparams["x_poly"],
            "y_poly": calparams["y_poly"],
            "x_poly_fwd": calparams["x_poly_fwd"],
            "y_poly_fwd": calparams["y_poly_fwd"],
            "x_poly_rev": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "y_poly_rev": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "star_list": star_list,
            "elev": calparams["device_alt"],
            "equal_aspect": False,
            "extinction_scale": 0.6,
            "force_distortion_centre": False,
            "fov_h": float(calparams["fieldw"]),
            "fov_v": float(calparams["fieldh"]),
            "gamma": 1.0,
            "lat": float(calparams["device_lat"]),
            "lon": float(calparams["device_lng"]),
            "mag_0": -2.5,
            "mag_lev": 16.796875,
            "mag_lev_stddev": 0.770936242463839,
            "measurement_apparent_to_true_refraction": False,
            "poly_length": 14,
            "pos_angle_ref": calparams["orig_pos_ang"],
            "refraction": True,
            #"rotation_from_horiz": -1.4170180941768553,
            "station_code": "XX0001",
            "version": 2,
            "vignetting_coeff": 0.001,
            "vignetting_fixed": True,
        }

        platepar_file = "output/%s-platepar.json" % dt.strftime("%Y_%m_%d_%H_%M_%S_000_1")
        with open(platepar_file, "w") as pp:
            json.dump(platepar, pp)
        return platepar_file

    def calibrate(self, c=[], i=[], iter=0):
        if len(self.star_pairs) > 0:
            c = np.array([(np.float64(p["catalog_star"][2]), np.float64(p["catalog_star"][3]), 0) for p in self.star_pairs])
            i = np.array([(np.float64(p["image_star"][0]), np.float64(p["image_star"][1]), 0) for p in self.star_pairs])
        if iter < 1:
            self.platepar.distortion_type = "radial3-odd"
        elif iter < 3:
            self.platepar.distortion_type = "radial5-odd"
        else:
            self.platepar.distortion_type = "radial7-odd"
        self.platepar.fitAstrometry(self.platepar.JD, i, c, first_platepar_fit=True, fit_only_pointing=False, fixed_scale=False)
        self.platepar.fitAstrometry(self.platepar.JD, i, c, first_platepar_fit=False, fit_only_pointing=False, fixed_scale=False)
        return True

    def get_pos(self):
        return self.platepar.lat, self.platepar.lon, self.platepar.elev, self.timestamp
    
    def set_star_pairs(self, pairs):
        self.star_pairs = pairs

    def get_star_list(self):
        return [(s[1], s[2], s[4], s[5]) for s in self.platepar.star_list]

    def ra_dec_to_xy(self, ra, dec):
        return raDecToXYPP(ra, dec, self.platepar.JD, self.platepar)
    
    def xy_to_ra_dec(self, x, y):
        _, ras, decs, _ = xyToRaDecPP(len(x)*[self.platepar.JD], x, y, len(x)*[1], self.platepar, extinction_correction=False, jd_time=True)
        return ras, decs
