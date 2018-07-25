

import numpy as np
import logging
import os
from astropy.io import fits
from astropy.table import (Table, vstack)

logger = logging.getLogger(__name__)

sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(logging.Formatter('%(asctime)s: %(message)s'))

logger.addHandler(sh)
logger.setLevel(logging.INFO)


def symbolic_link(path, destination_folder):
    source = os.path.abspath(path)
    basename = os.path.basename(path)
    destination = os.path.abspath(os.path.join(destination_folder, basename))

    if not os.path.exists(destination_folder):
        logger.info("Creating directories {}".format(destination_folder))
        os.makedirs(destination_folder)

    if not os.path.exists(destination):
        logger.info("Creating link {} -> {}".format(source, destination))
        os.symlink(source, destination)

    else:
        logger.info("Skipping link {} -> {} because destination exists"\
                    .format(source, destination))

    return None


class ImageCollection(object):

    _parse_headers = ("object", "obstype", "obsclass", "gemprgid", "telescop",
                      "ra", "dec", "ut", "date", "detsize", "grating", 
                      "centwave", "obsepoch")

    def __init__(self, paths):
        """
        A collection of images, with some useful features.

        :param paths:
            The file system paths to the images.
        """

        rows = []
        for path in paths:
            row = dict(path=path)
            with fits.open(path) as image:
                for key in self._parse_headers:
                    row[key] = image[0].header[key]
            rows.append(row)

        self._data = Table(rows=rows)
        return None


    @property
    def data(self):
        return self._data


    def prepare(self, output_folder, folder_pattern=None, 
                group_by=("object", "date", "grating"), 
                cal_constraints=("grating", )):

        if folder_pattern is  None:
            folder_pattern = "{object}/{date}"

        logger.info("Preparing folder structure in {}".format(output_folder))
        logger.info("Sub-folder pattern will be {}".format(folder_pattern))

        subset = (self.data["obstype"] == "OBJECT") \
               * (self.data["grating"] != "MIRROR") \
               * (self.data["obsclass"] != "acq")

        associations = self.data[subset].group_by(group_by)

        logger.info("{} associations found".format(len(associations)))

        wds = dict()

        for exposures in associations.groups:

            logger.info("On association: {}".format(exposures))

            row = exposures[0]
            row = dict(zip(row.colnames, row.as_void()))

            folder = os.path.join(output_folder, folder_pattern.format(**row))

            # Create symbolic links for data.
            for path in exposures["path"]:
                symbolic_link(path, folder)

            # Get closest arcs and flat for each exposure, and ensure that is
            # linked to the same folder.
            N = len(self.data)
            idx = np.arange(N)
            science_exposures = (exposures["obsclass"] == "science") \
                              + (exposures["obsclass"] == "partnerCal")

            all_indices = []

            for exposure in exposures[science_exposures]:

                meets_constraints = np.ones(N, dtype=bool)
                for key in cal_constraints:
                    meets_constraints *= (self.data[key] == exposure[key])

                arcs = meets_constraints * (self.data["obstype"] == "ARC")
                flats = meets_constraints * (self.data["obstype"] == "FLAT")
                biases = meets_constraints * (self.data["obstype"] == "BIAS")

                # Get closest arc/flat in time to current exposure.
                dt = np.abs(exposure["obsepoch"] - self.data["obsepoch"])

                indices = np.hstack([
                    idx[arcs][np.argmin(dt[arcs])],
                    idx[flats][np.argmin(dt[flats])],
                    idx[biases] # all biases, not just closest.
                ])

                # Create symbolic links.
                for index in indices:
                    symbolic_link(self.data["path"][index], folder)

                all_indices.extend(indices)

            all_indices = np.unique(all_indices)

            wds[folder] = vstack([exposures, self.data[all_indices]])

        return wds

