import os
import datetime
import ewok

__all__ = ["fc_file", "obs_file", "r2d2_obsfile", "r2d2_anfile", "r2d2_fixfile"]


def fc_file(fcout, step):
    fc = {}
    fc['date'] = fcout['date']
    keys = [fcout['exp'], fcout['type'], fcout['date'], ewok.jediformat(step)]
    fname = '.'.join(keys)
    fc['filename'] = os.path.join(fcout['datadir'], fname)
    return fc


def obs_file(conf):
    obsfile = conf['obsdatain']
    return obsfile


def r2d2_obsfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = ['fv3jedi', conf['source'], sdate, 'nc']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file


def _get_restarts(conf, date):
    dt_step = ewok.parse_timedelta(conf['step_cycle'])
    bdate = date - dt_step / 2
    rdate = bdate.strftime('%Y%m%d.%H')

    ntiles = conf.get('ntiles', 6)

    ftypes = conf.get('ftypes', ['fv_core.res', 'fv_srf_wnd.res', 'fv_tracer.res', 'sfc_data'])
    ftiles = [f'tile{t+1}' for t in range(0, conf.get('ntiles', 6)]

    restarts = []
    for ftype in ftypes:
        for tile in ftiles:
            #fname = f'{rdate}0000.{ftype}.{tile}.nc'
            fname = f'{ftype}.{tile}.nc'
            restarts.append(fname)
    return restarts


def _get_fv3files(conf):

    npz = conf.get('npz', 64)
    res = conf.get('resolution', 'c12')

    flist = conf.get(fv3filesList, [f'akbk{npz}.nc4', 'field_table', 'fmsmpp.nml', f'input_gfs_{res}.nml', 'inputpert_4dvar.nml'])

    fv3files = []
    for fname in flist:
        fv3files.append(fname)

    return fv3files


def _get_fieldsets(conf):

    flist = conf.get(fieldsetsList, ['dynamics.yaml'])

    fieldsets = []
    for fname in flist:
        fieldsets.append(fname)

    return fieldsets


def r2d2_anfile(conf, date):

    anlfileDict = {}
    restarts = _get_restarts(conf, date)
    fv3files = _get_fieldsets(conf)
    fieldsets = _get_fieldsets(conf)

    anlfileDict = ['restarts' : restarts,
                   'fv3files' : fv3files,
                   'fieldsets': fieldsets]

    return anlfileDict
