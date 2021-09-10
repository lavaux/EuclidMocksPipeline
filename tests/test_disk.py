import numpy as np
import euclid_obssys as eo
import euclid_obssys.disk as eo_d
import pytest
import os

@pytest.fixture(scope="session")
def temp_catalog_name(tmpdir_factory):
    tmp_fitsfile = tmpdir_factory.mktemp("data").join("test.fits")
    return tmp_fitsfile

@pytest.fixture
def catalog_write(temp_catalog_name):
    with eo_d.CatalogFitsWrite(temp_catalog_name) as cat:
        yield cat

@pytest.fixture
def catalog_read(temp_catalog_name):
    with eo_d.CatalogFitsRead() as cat:
        yield cat


@pytest.fixture(autouse=True)
def test_array_write(catalog_write):
    a = catalog_write.new_array("array1", (100,), dtype=[('a',float)])
    a[:] = np.arange(100)

@pytest.fixture(autouse=True)
def test_tag_write(catalog_write):
    catalog_write.add_tag('hello','world')

def test_array_read(temp_catalog_name):
    with eo_d.CatalogFitsRead("test.fits") as cat:
      assert cat.get_tag('hello') == 'world'
      a = cat.get_array('array1')['a']
      assert (np.equal(a, np.arange(100))).all()
