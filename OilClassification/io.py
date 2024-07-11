from os import path

__all__ = ['parse_prodname_UAVSAR']

def parse_prodname_UAVSAR(filename):
    """parse_prodname_UAVSAR(filename)

    Given a filename using the UAVSAR product naming convention, returns different
    parameters if found.  If not found returns [] for that field.  As written, this
    assumes products are for 90Deg beam steering.

    The UAVSAR filename convention is:
    {site}_{lineID}_{fltID}_{datatake}_{fltdate}_L090([HV]{4}|.?*)_XX_{int}.*
    """
    ProdPolarizations = ['HHHH', 'HVHV', 'VVVV', 'HHVV', 'HVVV', 'HHHV']

    # Get the UAVSAR product information (dataName and product name)

    #prod = [join(dirname(fn), basename(fn).split('.')[0]) for fn in filename] # full path to first '.'
    filename = path.basename(filename).split('.')[0]

    assert('L090' in filename) # Assumes products with 90 deg beam steering
    pol = filename.split('L090')[1][0:4] # first 4 chars after 'L090'
    if pol in ProdPolarizations:
        dataName = filename.replace(thispol,'')
    else:
        dataName = filename
        pol      = None

    _, lineID, fltID, datatake, fltdate, *rest = dataName.split('_')

    return [dataName, filename, pol, lineID, fltID, datatake, fltdate]
