################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from astropy.io import fits as astrofits
from fitsio import FITS, FITSHDR
import numpy as np
import sys
import healpy as hp
# import lxml.etree
# import lxml.builder
import filenames

if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)



# def write_xml(xml_info, out_name) :
#     '''write XML product associated to SEL-ID selected catalog
#     @par xml_info: dict with XML tags info
#     @par out_name: name and path of output xml file'''

#     def ATTR(key, value):
#         '''Utility to generate a single key-value dict to be used ad elemet attribute'''
#         return {key: value}


#     # get base product element with namespect and alias
#     PROD = lxml.builder.ElementMaker(namespace='http://euclid.esa.org/schema/dpd/le3/id/sel/out/subsampledcatalog',
#                                             nsmap={'subcat':'http://euclid.esa.org/schema/dpd/le3/id/sel/out/subsampledcatalog'}).PROD

#     # generic elemen builder
#     E = lxml.builder.ElementMaker()

#     # write all XML following structure
#     xml = PROD(
#         E.Header(
#         E.ProductId(xml_info['prod_id']),
#         E.ProductType('DpdLE3IDSELIDSubsampledCatalog'),
#         E.SoftwareName('Mock'),
#         E.SoftwareRelease('None'),
#         E.ProdSDC('SDC-IT'),
#         E.DataSetRelease(xml_info['release']),
#         E.Purpose('TEST'),
#         E.PlanId('None'),
#         E.PPOId('None'),
#         E.PipelineDefinitionId('None'),
#         E.PipelineRun('None'),
#         E.ExitStatusCode('OK'),
#         E.ManualValidationStatus('VALID'),
#         E.Curator('OU-LE3'),
#         E.CreationDate(xml_info['date'])
#         ),
#         E.Data(
#             E.SpatialCoverage(
#                 E.Polygon(
#                     E.Vertex(
#                         E.Position(
#                             E.C1('0'),
#                             E.C2('0')
#                         )
#                     )
#                 )
#             ),
#             E.CatalogDescription(
#                 E.CatalogOrigin(xml_info['catalog_XML_type']),
#                 E.CatalogType('PROXY'),  # NOT_PROXY if real observations
#                 E.PathToCatalogFile(),
#                 E.CatalogName(xml_info['cat_name'])
#             ),
#             E.Catalog(ATTR('format','le3.id.catalog.selected'), ATTR('version','0.1'),
#                 E.DataContainer(ATTR('filestatus','PROPOSED'),
#                     E.FileName(xml_info['file_name'])
#                 )
#             ),
#             E.Selection(xml_info['selection'])  
#         )
#     )

#     # finalize and write to file
#     with open(out_name, 'wb') as out :
#         #cleanup
#         clean_xml_stream = lxml.etree.tostring(xml, pretty_print=True, xml_declaration=True, encoding='utf-8')
#         # write as b-stream
#         out.write(clean_xml_stream)
# # ---------------------------------------------------------------------------- #




def write_catalog(fname, x, y, z, w, type="DATA", format="fits", coord="PSEUDO_EQUATORIAL"):
    if format=="txt":
        np.savetxt(fname, np.transpose( [x, y, z, w] ) )

    elif format=="fits":
        # current DM supported
        DM_VERSION = "2.5.0"

        # Extracted from AsciiToFits.py written by Daniele Tavagnacco
        
        if coord=="PSEUDO_EQUATORIAL":
            columns = [['SOURCE_ID', -1, 'f8'],
                       ['RA',         x, 'f8'],
                       ['DEC',        y, 'f8'],
                       ['REDSHIFT',   z, 'f8'],
                       ['WEIGHT',    -1, 'f8'],
                       ['DENSITY',    w, 'f8']]
        elif coord=="CARTESIAN":
            columns = [['SOURCE_ID', -1, 'f8'],
                       ['COMOV_X',    x, 'f8'],
                       ['COMOV_Y',    y, 'f8'],
                       ['COMOV_Z',    z, 'f8'],
                       ['WEIGHT',    -1, 'f8'],
                       ['DENSITY',    w, 'f8']]
 
        header_keywords = {"TELESCOP": "EUCLID",
                           "INSTRUME": "LE3GC-MOCKS",
                           "FILENAME": fname,
                           "CAT_ID"  : "MOCK",
                           "COORD"   : coord,
                           "ANGLE"   : "DEGREES"}

        # fname_nodir = fname[[pos for pos, char in enumerate(fname) if char == '/'][-1]+1:]
        # zshell = fname[fname.find('zshell'):-5]
        # header_keywords = {'TELESCOP' : 'EUCLID  ',
        #                    'INSTRUME' : 'LE3IDSELID',
        #                    'FILENAME' : fname_nodir,
        #                    'CAT_TYPE' : type,
        #                    'CAT_NAME' : 'MOCK-LE3GC',
        #                    'CAT_NOBJ' : x.size,
        #                    'COORD   ' : 'EQUATORIAL',
        #                    'ANGLE   ' : 'DEG     ',
        #                    'COMPLETE' : 1.0,
        #                    'PURITY  ' : 1.0,
        #                    'SELECT  ' : zshell}

        extension="CATALOG"

        xmlKeys = {"pf"    : "PK_LE3_GC_WindowMultipoles",
                   "instr" : "LE3_GC_MOCKS",
                   "id"    : "MOCK",
                   "coord" : coord}

        print("Preparing FITS structure")
        types = []
        # keep just wanted columns ...bad but..
        for c in columns :
            types.append((c[0], c[2]))

        hdr = FITSHDR()
        print("+ Add keywords")
        for k in header_keywords :
            hdr[k] = header_keywords[k]

        keep_table = {}
        for c in columns :
            # if required but not existing (-1) fill with ones
            if (c[1] is -1) :
                # add tmp column with correct name and position, bu only ones
                print(str("+  -Col '%s' filled" % c[0]))
                keep_table[c[0]] = np.ones_like(columns[1][1],dtype=np.float64)
            else :
                # keep column with  requested position in the input file
                print(str("+  -Col '%s' from '%s'" % (c[0], c[1])))
                keep_table[c[0]] = c[1].astype(np.float64)

        fullsize = len(keep_table)*len(c[1])* 8 / 1024 / 1024.
        print(str("+   ~%.2f MB in memory" % fullsize))

        # now write some data
        print(str("+ Write FITS: %s" % fname))
        fits = FITS(fname,'rw', clobber=True)
        fits.write_table(data=keep_table, header=hdr, extname=extension)
        fits.close()

        print("+ Preparing XML product")

        with open (fname.replace(input.cat4le3_format,"xml"), "w+") as f :

            f.write('''<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n''')
            if (type == "RANDOM") :
                f.write('''<p1:DpdLE3gcInputRandCat xmlns:p1="http://euclid.esa.org/schema/dpd/le3/gc/inp/catrandin">\n''')
            else:
                f.write('''<p1:DpdLE3gcInputDataCat xmlns:p1="http://euclid.esa.org/schema/dpd/le3/gc/inp/catdatain">\n''')
            f.write('''  <Header>\n''')
            f.write(str('''    <ProductId>%s</ProductId>\n''' % header_keywords["FILENAME"].split('.')[0]))
            f.write('''    <ProductType>dpdLE3gcInputRandCat</ProductType>\n''')
            f.write('''    <SoftwareName>LE3_GC_test</SoftwareName>\n''')
            f.write('''    <SoftwareRelease>1.0</SoftwareRelease>\n''')
            f.write('''    <ManualValidationStatus>UNKNOWN</ManualValidationStatus>\n''')
            f.write('''    <PipelineRun>LE3_GC_Test_Inputs</PipelineRun>\n''')
            f.write('''    <ExitStatusCode>OK</ExitStatusCode>\n''')
            f.write(str('''    <DataModelVersion>%s</DataModelVersion>\n''' % DM_VERSION))
            f.write(str('''    <MinDataModelVersion>%s</MinDataModelVersion>\n''' % DM_VERSION))
            f.write('''    <ScientificCustodian>LE3</ScientificCustodian>\n''')
            f.write('''    <AccessRights>\n''')
            f.write('''      <EuclidConsortiumRead>true</EuclidConsortiumRead>\n''')
            f.write('''      <EuclidConsortiumWrite>true</EuclidConsortiumWrite>\n''')
            f.write('''      <ScientificGroupRead>true</ScientificGroupRead>\n''')
            f.write('''      <ScientificGroupWrite>true</ScientificGroupWrite>\n''')
            f.write('''    </AccessRights>\n''')
            f.write('''    <Curator>\n''')
            f.write('''      <Name>SDC-IT</Name>\n''')
            f.write('''    </Curator>\n''')
            f.write(str('''    <Creator>%s</Creator>\n''' % xmlKeys["pf"]))
            f.write('''    <CreationDate>2019-10-31T12:12:12Z</CreationDate>\n''')
            f.write('''  </Header>\n''')
            f.write('''  <Data>\n''')
            f.write(str('''  <Instrument>%s</Instrument>\n''' % xmlKeys["instr"]))
            f.write(str('''  <Catalog_ID>%s</Catalog_ID>\n'''% xmlKeys["id"]))
            f.write(str('''  <CoordType>%s</CoordType>\n''' % xmlKeys["coord"]))
            f.write('''  <Catalog format="le3.gc.cat.test" version="0.2">\n''')
            f.write('''    <DataContainer filestatus="PROPOSED">\n''')
            f.write(str('''      <FileName>%s.fits</FileName>\n''' % fname))
            f.write('''    </DataContainer>\n''')
            f.write('''  </Catalog>\n''')
            f.write('''  </Data>\n''')
            if (type == "RANDOM") :
                f.write('''</p1:DpdLE3gcInputRandCat>\n''')
            else :
                f.write('''</p1:DpdLE3gcInputDataCat>\n''')
            f.close()

        print("files %s and %s written"%(fname,fname.replace(input.cat4le3_format,"xml")))

    else:
        
        print("ERROR: unrecognized format in write_catalog")
        sys.exit(-1)

###################################################################################


print("# Running writeCatalog4LE3.py with {}".format(sys.argv[1]))

r1=input.pinocchio_first_run
r2=input.pinocchio_last_run
if input.cat_type is 'pinocchio':
    n1=r1
    if r2 is not None:
        n2=r2
    else:
        n2=n1
else:
    n1=n2=0

fname = filenames.dndz(input)
print("# Reading dndz from {}".format(fname))
dndz = astrofits.getdata(fname)

# you can skip the writing of the random
if input.WriteLE3Random:
    print("# reading random catalog {}...".format(filenames.random(input)))
    randomcat = astrofits.getdata(filenames.random(input))

    if (not input.apply_dataselection_to_random) & (input.selection_random_tag is not None):
        fname=filenames.selection_random(input)
        print("# reading selection of random from file {}...".format(fname))
        selection = astrofits.getdata(fname)['SELECTION']
    else:
        selection = np.ones(len(randomcat), dtype=bool)

    for zshell in input.finalCatZShell:

        zmin = zshell[0]
        zmax = zshell[1]

        print("# selecting the shell at z=%f-%f..."%(zmin,zmax))
        sel = selection & (randomcat[input.redshift_key] >= zmin) & (randomcat[input.redshift_key] < zmax)
        print("# computing weights...")
        density = np.interp(randomcat[input.redshift_key][sel], dndz['z_center'], dndz['N_gal']/dndz['bin_volume'])

        print("# writing random file {}".format(filenames.LE3_random(input,zmin,zmax)))
        write_catalog(filenames.LE3_random(input,zmin,zmax),
                      randomcat['ra_gal'][sel], randomcat['dec_gal'][sel], 
                      randomcat[input.redshift_key][sel],
                      density, type = "RANDOM", format = 'fits')

    del randomcat


if input.cat_type is not 'pinocchio':
    toprocess=[None]
else:
    toprocess=range(n1,n2+1)

for myrun in toprocess:
    fname=filenames.galcat(input,myrun)
    print("# reading data catalog {}...".format(fname))
    cat = astrofits.getdata(fname)
    zused = cat[input.redshift_key]
    if input.selection_data_tag is not None:
        fname = filenames.selection_data(input,myrun)
        print("# loading selection {}...".format(fname))
        selection = astrofits.getdata(fname)['SELECTION']
    else:
        selection = np.ones(zused.size,dtype=bool)

    for zshell in input.finalCatZShell:

        zmin = zshell[0]
        zmax = zshell[1]

        print("# selecting the shell at z=%f-%f..."%(zmin,zmax))
        mysel = (zused >= zmin) & (zused < zmax) & selection
        print("# computing weights...")
        density = np.interp(zused[mysel], dndz['z_center'], dndz['N_gal']/dndz['bin_volume'])

        fname = filenames.LE3_data(input,zmin,zmax,myrun)
        print("# writing data file {}".format(fname))
        write_catalog(fname, cat['ra_gal'][mysel], cat['dec_gal'][mysel], zused[mysel], 
                      density, type = "DATA", format = 'fits')


print("# DONE!")
