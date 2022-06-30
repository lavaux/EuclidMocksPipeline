# *****************************************************************
# *                        PINOCCHIO  V4.1                        *
# *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
# *****************************************************************

# This code was written by
# Pierluigi Monaco
# Dipartimento di Fisica, Universita` di Trieste
# Copyright (C) 2016

# web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

# The original code was developed with:
# Tom Theuns           Institute for Computational Cosmology, University of Durham
# Giuliano Taffoni     INAF - Osservatorio Astronomico di Trieste

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


"""
Routine for reading Pinocchio's binary catalogs and PLCs.
Usage:

# PLC
import ReadPinocchio as rp
myplc = rp.plc("pinocchio.example.plc.out")
print(myplc.redshift)

# CATALOGS
import ReadPinocchio as rp
mycat = rp.catalog("pinocchio.0.0000.example.catalog.out")
print(mycat.Mass)

# HISTORIES
import ReadPinocchio as rp
myhist = rp.histories("pinocchio.example.histories.out")
print(myhist.name)

Written by Emiliano Munari

LAST MODIFIED:
Pierlugi Monaco 15/04/2016
Guilhem Lavaux 2022

"""


import numpy as np
import os
import sys
from struct import unpack
import logging
import contextlib
import copy

BUFFERSIZE = 1024 * 1024**2

DataBuffer = None
FileLength = 0
FirstCall = True
ByteFileCounter = 0
ByteBufferCounter = 0
BytesRead = BUFFERSIZE


_log = logging.getLogger("pinocchio")

_log.setLevel(logging.INFO)


@contextlib.contextmanager
def my_open(filename):

    global FileLength
    global ByteFileCounter
    global ByteBufferCounter
    global FirstCall

    ByteFileCounter = 0
    ByteBufferCounter = 0
    FirstCall = True

    if not os.path.exists(filename):
        _log.error("file %s does not exist", filename)
        return None

    FileLength = os.path.getsize(filename)

    _log.debug("Opened file %s of length %d", filename, FileLength)

    f = open(filename, mode="rb")
    yield f
    f.close()

    # checks that the file has been fully read
    _log.debug("Closing file, read %d bytes out of %d", ByteFileCounter, FileLength)
    if ByteFileCounter == FileLength:
        _log.debug("The file has been read to the end")
    else:
        _log.debug(
            "Warning: there are still %d bytes to read", FileLength - ByteFileCounter
        )


def my_fromfile(f, dtype, count, add_fortran_tags=False):

    global FileLength
    global ByteFileCounter
    global ByteBufferCounter
    global DataBuffer
    global FirstCall
    global BytesRead

    # load data if the buffer is void
    if FirstCall:
        BytesRead = min(BUFFERSIZE, FileLength)
        _log.debug(
            "Loading the first chunk of data: BUFFERSIZE=%d, FileLength=%d, I will load %d bytes",
            BUFFERSIZE,
            FileLength,
            BytesRead,
        )
        DataBuffer = bytearray(f.read(BytesRead))
        FirstCall = False

    my_dt = copy.deepcopy(dtype)
    if add_fortran_tags:
        my_dt.insert(0, ("fort", np.int32))
        my_dt.append(("trof", np.int32))
    my_dt = np.dtype(my_dt)
    dtsize = my_dt.itemsize
    wanted = count * dtsize
    inBuffer = BytesRead - ByteBufferCounter

    _log.debug(
        "my_fromfile called to read %d fields of length %d, total: %d",
        count,
        dtsize,
        wanted,
    )

    # check that the request does not overshoot the content of the file
    if wanted > FileLength - ByteFileCounter:
        _log.error("this request overshoots the file content")
        return None

    # if the buffer does not contain all the needed information:
    # then shift the info to the top and reload the buffer
    if wanted > inBuffer:

        # read what is present
        npresent = inBuffer // dtsize
        _log.debug(
            "the buffer is insufficient, reading %d fields of size %d", npresent, dtsize
        )

        my_data = np.frombuffer(
            DataBuffer[ByteBufferCounter:], dtype=my_dt, count=npresent
        )
        loaded = npresent * dtsize
        ByteBufferCounter += loaded
        inBuffer -= loaded
        ByteFileCounter += loaded

        # update buffer
        added = min(ByteBufferCounter, FileLength - ByteFileCounter)

        _log.debug(
                "Shifting %d bytes to the start of the buffer and updating the buffer with %d more bytes",
                inBuffer,
                added,
            )

        # then load the remaining fields
        tmp = DataBuffer[ByteBufferCounter:]
        DataBuffer[:inBuffer] = tmp[:]
        DataBuffer[inBuffer : inBuffer + added] = bytearray(f.read(added))
        del tmp

        my_data = np.concatenate(
            (my_data, np.frombuffer(DataBuffer, dtype=my_dt, count=count - npresent)),
            axis=0,
        )

        moredata = (count - npresent) * dtsize
        ByteBufferCounter = moredata
        ByteFileCounter += moredata

    else:
        # otherwise just read the data from the buffer and update the pointers
        my_data = np.frombuffer(
            DataBuffer[ByteBufferCounter:], dtype=my_dt, count=count
        )
        ByteBufferCounter += wanted
        ByteFileCounter += wanted

    return my_data


class catalog:
    def __init__(self, filename):

        # checks that the filename contains 'catalog'
        if not 'catalog' in filename:
            _log.error("Are you sure you are providing the right file name?")
            return None

        # checks that the input file ends by ".out"
        last_ext=filename.rfind('.')
        if filename[last_ext:]!='.out':

            _log.error("The catalog file should end with .out, the file number extension will be checked by the code")
            return None


        if not os.path.exists(filename):
            if not os.path.exists(f"{filename}.0"):
                _log.error("file not found: %s", filename)
                return None
            else:
                Nfiles = 1
                while os.path.exists(f"{filename}.{Nfiles}"):
                    Nfiles+=1
                _log.info("The catalog is written in {} files".format(Nfiles))
        else:
            Nfiles = 1

        _log.debug("opening file %s", filename)
        if Nfiles == 1:
            fname = filename
        else:
            fname = f"{filename}.0"

        with my_open(fname, mode="rb") as f:
            reading = np.fromfile(f, dtype=np.int32, count=10)

        NTasksPerFile = reading[1]
        _log.debug('This file has been written by {} tasks', NTasksPerFile)
        if reading[2]>10:
            newRun=True
            record_length=reading[2]
            Nslices=1
            _log.info('This is new output format, record length: {}', record_length)
        else:
            newRun=False
            record_length=reading[7]
            Nslices=reading[2]
            _log.info('This is classic output format, record length: {}', record_length)
            _log.debug('The box has been fragmented in {} slice{}', Nslices, '' if Nslices==1 else 's')
        # sets the record
        if record_length==96:

            self.cat_dtype=[ ('name',  np.int64),
                             ('Mass',  np.float64),
                             ('posin', np.float64, 3),
                             ('pos',   np.float64, 3),
                             ('vel',   np.float64, 3),
                             ('npart', np.int32),
                             ('pad'  , np.int32) ]

        elif record_length==56:

            if newRun:

                self.cat_dtype=[ ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('posin', np.float32, 3),
                                 ('npart', np.int32),
                                 ('pad'  , np.int32) ]

            else:

                self.cat_dtype=[ ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('posin', np.float32, 3),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('npart', np.int32),
                                 ('pad'  , np.int32) ]

        elif record_length==48:

            self.cat_dtype=[ ('name',  np.int64),
                             ('Mass',  np.float32),
                             ('pos',   np.float32, 3),
                             ('vel',   np.float32, 3),
                             ('posin', np.float32, 3) ]

        else:
            _log.error("sorry, I do not recognize this record length")
            return None

        # decides what files to read
        if Nfiles>1:
            if first_file is None:
                first_file=0
            elif first_file<0:
                first_file=0
            elif first_file>Nfiles:
                first_file=Nfiles
            if last_file is None:
                last_file=Nfiles
            else:
                last_file += 1
                if last_file<first_file:
                    last_file=first_file+1
                elif last_file>Nfiles:
                    last_file=Nfiles
            # this is to be used to define a pythonic range
            if not silent:
                print("I will read files in the python range from {} to {}".format(first_file,last_file))
        else:
            first_file=0
            last_file=Nfiles

        # counts the number of groups
        # loop on files
        Ngroups = np.int64(0)
        for myfile in range(first_file,last_file):

            if Nfiles==1:
                myfname = filename
            else:
                myfname = "{filename}.{myfile}"

            _log.info("reading file {}",myfname)

            with my_open(myfname, mode="rb") as f:
                Nthisfile=0

                # moves past the header
                f.seek(16,1)

                for myslice in range(Nslices):

                    _log.debug("...slice {}", myslice)

                    WritingTasks=0
                    for mytask in range(NTasksPerFile):

                        _log.debug("...task {}",mytask)

                        ngood = np.fromfile(f,dtype=np.int32,count=3)[1]
                        Ngroups+=ngood
                        Nthisfile+=ngood
                        if ngood>0:
                            WritingTasks+=1

                        _log.debug("found {} groups",ngood)
                        if newRun:
                            if ngood>0:
                                f.seek(ngood*record_length+8,1)
                        else:
                            f.seek(ngood*(record_length+8),1)

            if newRun:
                file_length = 16 + NTasksPerFile*12 + WritingTasks*8 + Nthisfile*record_length
            else:
                file_length = 16 + Nslices*NTasksPerFile*12 + Nthisfile*(record_length+8)
            if FileLength != file_length:
                _log.error("Error: inconsistent length for the file {}, I expected {} and found {}",myfname,file_length,FileLength)
                return None
            elif VERBOSE:
                _log.debug("File size of {} is as expected",myfname)

        _log.info("Number of halos in the catalog: {}", Ngroups)

        # allocates data
        self.data = np.empty(Ngroups, dtype=self.cat_dtype)

        _log.debug("Reading the catalog")

        # counts the number of groups
        # loop on files
        counter=0
        for myfile in range(first_file,last_file):

            if Nfiles==1:
                myfname = filename
            else:
                myfname = f"{filename}.{myfile}"

            _log.info("reading file {}",myfname)

            with my_open(fname, mode="rb") as f:

                # skips the header
                skip = my_fromfile(f,np.int32,4)

                for myslice in range(Nslices):

                    _log.debug("...slice {}".format(myslice))

                    for mytask in range(NTasksPerFile):
                        _log.debug("...task {}".format(mytask))

                        ngood = my_fromfile(f,np.int32,3)[1]

                        if ngood>0:
                            if newRun:
                                # reads and checks the fortran tag
                                size = my_fromfile(f,np.int32,1)[0]
                                if size != ngood * np.dtype(self.cat_dtype).itemsize:
                                    _log.error("ERROR: inconsistency in catalog file")
                                    return None

                            data = my_fromfile(f,self.cat_dtype,ngood,add_fortran_tags=not newRun)

                            if newRun:
                                # reads and checks the fortran tag
                                size = my_fromfile(f,np.int32,1)[0]
                                if size != ngood * np.dtype(self.cat_dtype).itemsize:
                                    _log.error("ERROR: inconsistency in catalog file")
                                    return None

                            self.data['name'] [counter:counter+ngood] = data['name']
                            self.data['Mass'] [counter:counter+ngood] = data['Mass']
                            self.data['pos']  [counter:counter+ngood] = data['pos']
                            self.data['vel']  [counter:counter+ngood] = data['vel']
                            self.data['posin'][counter:counter+ngood] = data['posin']
                            if record_length > 48:
                                self.data[counter:counter+ngood]['npart'] = data['npart']

                            counter+=ngood

                _log.debug(f"done with file {}", myfname)


class plc:
    def __init__(self, filename, first_file=None, last_file=None):

        # checks that the filename contains 'plc'
        if not "plc" in filename:
            print("Are you sure you are providing the right file name?")
            return None

        # checks that the input file ends by ".out"
        last_ext = filename.rfind(".")
        if filename[last_ext:] != ".out":

            print(
                "The catalog file should end with .out, the file number extension will be checked by the code"
            )
            return None

        # checks that the file exists, of that there are multiple files
        if not os.path.exists(filename):

            if not os.path.exists(filename + ".0"):

                print("file {} or {} not found:".format(filename, filename + ".0"))
                return None

            else:

                Nfiles = 1
                while os.path.exists(filename + ".{}".format(Nfiles)):
                    Nfiles += 1
                _log.info("The catalog is written in %d files", Nfiles)

        else:

            Nfiles = 1
            _log.info("The catalog is written in 1 file")

        # opens the (first) file and reads the record length
        if Nfiles == 1:
            filename_real = filename
        else:
            filename_real = f"{filename}.0"

        _log.info("opening file %s", filename_real)

        with my_open(filename_real) as f:
            reading = my_fromfile(f, dtype=np.int32, count=3)

        if reading[0] == 4:
            record_length = reading[1]
            newRun = True
            if record_length == 32:
                _log.debug(
                    "This is new light output format, record length: %d", record_length
                )
            else:
                _log.debug(
                    "This is new full output format, record length: %d", record_length
                )
        else:
            record_length = reading[0]
            newRun = False
            _log.debug(
                "This is classic output format, record length: %d", record_length
            )

        # sets the record
        if record_length == 104:

            self.cat_dtype = [
                ("name", np.uint64),
                ("truez", np.float64),
                ("pos", np.float64, 3),
                ("vel", np.float64, 3),
                ("Mass", np.float64),
                ("theta", np.float64),
                ("phi", np.float64),
                ("vlos", np.float64),
                ("obsz", np.float64),
            ]

        elif record_length == 56:

            self.cat_dtype = [
                ("name", np.uint64),
                ("truez", np.float32),
                ("pos", np.float32, 3),
                ("vel", np.float32, 3),
                ("Mass", np.float32),
                ("theta", np.float32),
                ("phi", np.float32),
                ("vlos", np.float32),
                ("obsz", np.float32),
            ]

        elif record_length == 32:

            self.cat_dtype = [
                ("name", np.uint64),
                ("truez", np.float32),
                ("Mass", np.float32),
                ("theta", np.float32),
                ("phi", np.float32),
                ("obsz", np.float32),
                ("pad", np.float32),
            ]

        else:
            print("sorry, I do not recognize this record length")
            return None

        _log.info("Detected dtype for the catalog: %s", repr(self.cat_dtype))

        # decides what files to read
        if Nfiles > 1:
            if first_file is None:
                first_file = 0
            elif first_file < 0:
                first_file = 0
            elif first_file > Nfiles:
                first_file = Nfiles
            if last_file is None:
                last_file = Nfiles
            else:
                last_file += 1
                if last_file < first_file:
                    last_file = first_file + 1
                elif last_file > Nfiles:
                    last_file = Nfiles
            # this is to be used to define a pythonic range
            _log.debug(
                "I will read files in the python range from %d to %d",
                first_file,
                last_file,
            )
        else:
            first_file = 0
            last_file = Nfiles

        # counts the number of groups
        # loop on files
        Ngroups = np.int64(0)
        for myfile in range(first_file, last_file):

            if Nfiles == 1:
                myfname = filename
            else:
                myfname = filename + ".{}".format(myfile)

            _log.info("reading file %s", myfname)

            with my_open(myfname) as f:

                Nthisfile = 0

                if newRun:

                    f.seek(12, 0)
                    size = 12
                    while size < FileLength:
                        NN = np.fromfile(f, dtype=np.int32, count=3)[1]
                        Nthisfile += NN
                        f.seek(8 + NN * record_length, 1)
                        size += 20 + NN * record_length
                        _log.debug("this chunk contains %d halos", NN)

                    if size != FileLength:
                        _log.error(
                            "Error: inconsistent length for the file %s, I expected %d and found %d",
                            myfname,
                            size,
                            FileLength,
                        )
                        return None

                    _log.debug("File %s contains %d halos", myfname, Nthisfile)

                    Ngroups += Nthisfile

                else:

                    Nthisfile = FileLength // (record_length + 8)
                    Ngroups += Nthisfile

                    _log.debug("File %s contains %d halos", myfname, Nthisfile)

                    if FileLength % (record_length + 8) > 0:
                        _log.error(
                            "Error: inconsistent length for the file %s, %d is not a multiple of %d",
                            myfname,
                            FileLength,
                            record_length + 8,
                        )
                        return None

            _log.debug("File size of %s is as expected", myfname)

        _log.info("Number of halos in the catalog: %d", Ngroups)

        # allocates data
        self.data = np.empty(Ngroups, dtype=self.cat_dtype)

        # reads data from files
        pointer = 0
        for myfile in range(first_file, last_file):

            if Nfiles == 1:
                myfname = filename
            else:
                myfname = f"{filename}.{myfile}"

            _log.info("reading file %s", myfname)

            with my_open(myfname) as f:

                if newRun:

                    my_fromfile(f, np.int32, 3)

                    size = 12
                    read = 0
                    while size < FileLength:
                        NN = my_fromfile(f, np.int32, 4)[1]
                        data = my_fromfile(
                            f, self.cat_dtype, NN, add_fortran_tags=False
                        )
                        self.data["name"][pointer : pointer + NN] = data["name"]
                        self.data["truez"][pointer : pointer + NN] = data["truez"]
                        self.data["Mass"][pointer : pointer + NN] = data["Mass"]
                        self.data["theta"][pointer : pointer + NN] = data["theta"]
                        self.data["phi"][pointer : pointer + NN] = data["phi"]
                        self.data["obsz"][pointer : pointer + NN] = data["obsz"]
                        if record_length > 32:
                            self.data["pos"][pointer : pointer + NN] = data["pos"]
                            self.data["vel"][pointer : pointer + NN] = data["vel"]
                            self.data["vlos"][pointer : pointer + NN] = data["vlos"]

                        my_fromfile(f, np.int32, 1)
                        size += 20 + NN * record_length
                        read += NN
                        pointer += NN

                else:

                    Nthisfile = FileLength // (record_length + 8)

                    _log.debug("starting to read %d halos", Nthisfile)

                    read = 0
                    while read < Nthisfile:
                        NN = min(BUFFERSIZE // (record_length + 8), Nthisfile - read)
                        data = my_fromfile(f, self.cat_dtype, NN, add_fortran_tags=True)
                        self.data["name"][pointer : pointer + NN] = data["name"]
                        self.data["truez"][pointer : pointer + NN] = data["truez"]
                        self.data["Mass"][pointer : pointer + NN] = data["Mass"]
                        self.data["theta"][pointer : pointer + NN] = data["theta"]
                        self.data["phi"][pointer : pointer + NN] = data["phi"]
                        self.data["obsz"][pointer : pointer + NN] = data["obsz"]
                        self.data["pos"][pointer : pointer + NN] = data["pos"]
                        self.data["vel"][pointer : pointer + NN] = data["vel"]
                        self.data["vlos"][pointer : pointer + NN] = data["vlos"]
                        read += NN
                        pointer += NN

            _log.debug("Read %d halos from file %s", read, myfname)

        self.name = self.data["name"]
        self.truez = self.data["truez"]
        self.Mass = self.data["Mass"]
        self.theta = self.data["theta"]
        self.phi = self.data["phi"]
        self.obsz = self.data["obsz"]
        # Compatibility field
        self.redshift = self.truez
        if "pos" in self.data.dtype.fields:
            self.pos = self.data["pos"]
            self.vel = self.data["vel"]
            self.vlos = self.data["vlos"]


class histories:
    def __init__(self, filename):
        if not os.path.exists(filename):
            print("file not found:", filename)
            sys.exit()

        f = open(filename, "rb")

        Ngroups = np.int64(0)

        header = np.fromfile(f, dtype=np.int32, count=3)
        NSlices = header[1]

        if NSlices > 10:
            lite = True
            record_length = NSlices
            NSlices = 1
            print("This is light output format, record length: %d" % record_length)

            record_dtype = np.dtype(
                [
                    ("dummy", np.int32),
                    ("name", np.uint64),
                    ("nickname", np.int32),
                    ("link", np.int32),
                    ("merged_with", np.int32),
                    ("mass_at_merger", np.int32),
                    ("mass_of_main", np.int32),
                    ("z_merging", np.float32),
                    ("z_peak", np.float32),
                    ("z_appear", np.float32),
                    ("dummy2", np.int32),
                ]
            )
        else:
            lite = False
            record_length = -1
            if NSlices == 1:
                print("The box has been fragmented in 1 slice")
            else:
                print("The box has been fragmented in %d slices" % NSlices)

            record_dtype = np.dtype(
                [
                    ("dummy", np.int32),
                    ("name", np.uint64),
                    ("nickname", np.int32),
                    ("link", np.int32),
                    ("merged_with", np.int32),
                    ("mass_at_merger", np.int32),
                    ("mass_of_main", np.int32),
                    ("z_merging", np.float64),
                    ("z_peak", np.float64),
                    ("z_appear", np.float64),
                    ("pad", np.int32),
                    ("dummy2", np.int32),
                ]
            )

        # Count the number of halos
        Total = 0
        for islice in range(NSlices):

            header = np.fromfile(f, dtype=np.int32, count=4)
            Ntrees = header[1]
            Nbranches = header[2]
            Total += Nbranches

            print(
                "Slice N. ", islice, ": ", Ntrees, " trees and ", Nbranches, " branches"
            )

            for itree in range(Ntrees):
                header = np.fromfile(f, dtype=np.int32, count=4)
                mytree = header[1]
                mynbranch = header[2]

                f.seek(mynbranch * record_dtype.itemsize, 1)

        # Go back to the starting point (NTasksPerFile already read)
        print("+++++++++++++++++++")
        f.seek(12)

        print("Total number of halos: ", Total)

        self.name = np.empty(Total, dtype=np.uint64)
        self.nickname = np.empty(Total, dtype=np.int32)
        self.link = np.empty(Total, dtype=np.int32)
        self.merged_with = np.empty(Total, dtype=np.int32)
        self.mass_at_merger = np.empty(Total, dtype=np.int32)
        self.mass_of_main = np.empty(Total, dtype=np.int32)
        self.z_merging = np.empty(Total, dtype=np.float64)
        self.z_appear = np.empty(Total, dtype=np.float64)
        self.z_peak = np.empty(Total, dtype=np.float64)

        startid = 0
        stopid = 0
        for islice in range(NSlices):
            header = np.fromfile(f, dtype=np.int32, count=4)
            Ntrees = header[1]
            Nbranches = header[2]

            for itree in range(Ntrees):
                header = np.fromfile(f, dtype=np.int32, count=4)
                mytree = header[1]
                mynbranch = header[2]

                stopid += mynbranch
                catalog = np.fromfile(f, dtype=record_dtype, count=mynbranch)

                self.name[startid:stopid] = catalog["name"]
                self.nickname[startid:stopid] = catalog["nickname"]
                self.link[startid:stopid] = catalog["link"]
                self.merged_with[startid:stopid] = catalog["merged_with"]
                self.mass_at_merger[startid:stopid] = catalog["mass_at_merger"]
                self.mass_of_main[startid:stopid] = catalog["mass_of_main"]
                self.z_merging[startid:stopid] = catalog["z_merging"]
                self.z_appear[startid:stopid] = catalog["z_appear"]
                self.z_peak[startid:stopid] = catalog["z_peak"]
                del catalog
                startid = stopid

        f.close()
