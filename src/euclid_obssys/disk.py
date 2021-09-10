from abc import ABC, abstractmethod
from astropy.io import fits
from contextlib import AbstractContextManager
import numpy as np
import numpy.typing as npt


class AbstractCatalogWrite(ABC):
    @abstractmethod
    def new_array(self, name: str, shape: int, dtype) -> npt.ArrayLike:
        raise NotImplementedError()

    @abstractmethod
    def set_array(self, name: str, array: npt.ArrayLike) -> None:
        raise NotImplementedError()


class AbstractCatalogRead(ABC):
    @abstractmethod
    def get_tag(self, name: str) -> str:
        raise NotImplementedError()

    @abstractmethod
    def get_array(self, name: str) -> npt.ArrayLike:
        raise NotImplementedError()

    @abstractmethod
    def __getitem__(self, name: str) -> npt.ArrayLike:
        raise NotImplementedError()


class CatalogFitsWrite(AbstractCatalogWrite, AbstractContextManager):
    def __init__(self, fname: str):
        self.fname = fname
        self.arrays = {}
        self.tags = {}
        self.exited = False

    def __exit__(self, exc_type, exc_value, traceback):
        self.exited = True

        primary_hdu = fits.PrimaryHDU()
        for tag in self.tags.keys():
            primary_hdu.header.append((tag, self.tags[tag]))

        primary_hdu.header.append(("OBSSYS", "1"))

        hdul = fits.HDUList([primary_hdu])

        table_id = 1
        for array_name in self.arrays.keys():
            hdul.append(fits.BinTableHDU(data=self.arrays[array_name]))
            primary_hdu.header.append((f"TN{table_id}", array_name))

        basepath, _ = os.path.split(self.fname)
        if not os.path.exists(basepath):
            os.makedirs(basepath)

        hdul.writeto(self.fname, overwrite=True)
        return None

    def new_array(self, name: str, shape: int, dtype: npt.DTypeLike) -> npt.ArrayLike:
        assert self.exited == False
        a = np.empty(shape, dtype=dtype)
        self.arrays[name] = a
        return a

    def set_array(self, name: str, array: npt.ArrayLike) -> None:
        assert self.exited == False
        self.arrays[name] = array

    def add_tag(self, tagname: str, tagdata: str):
        if type(tagdata) is not str:
            tagdata = repr(tagdata)

        self.tags[tagname] = tagdata


class CatalogFitsRead(AbstractCatalogRead, AbstractContextManager):
    def __init__(self, basepath: str):
        self.basepath = basepath
        self.arrays = None

    def __enter__(self):
        self.file = fits.open(self.basepath)
        self.header = self.file[0].header
        if "OBSSYS" in self.header and self.header["OBSSYS"] == "1":
            self.arrays = {}
            for i in range(1, len(self.file)):
                self.arrays[self.header[f"TN{i}"]] = i
        return self

    def get_tag(self, name: str) -> str:
        return self.header[name]

    def get_array(self, name: str) -> npt.ArrayLike:
        if not self.arrays is None:
            return self.file[self.arrays[name]].data[:]
        else:
            return self.file[int(name)].data[:]

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()
        return None

    def __getitem__(self, name: str) -> npt.ArrayLike:
        return self.get_array(name)


DefaultCatalogWrite = CatalogFitsWrite
DefaultCatalogRead = CatalogFitsRead
