#!/usr/bin/env python3
"""
txt2frognc.py -- convert column-format text forcing files into netCDF files
readable by FROG's get_clim_forcing.

FROG expects a forcing file that satisfies, in spatialvars_mod.f90:

    nDims == 3            lat, lon, time
    nVars >= 4            lat, lon, time + at least one data variable
    unlimdimid != -1      time must be the UNLIMITED dimension

and reads the data variable as (time, lat, lon). The variable is looked up by
name (`name_tas_variable` in the namelist), so any CF-compliant name works as
long as the namelist matches.

The input text files are one value per line, one line per time step, with no
header. Windows line endings and UTF-8 BOMs are tolerated.

Real-world records are usually on a proleptic Gregorian calendar and so
contain 29 February. FROG has no such calendar: timer_mod knows only 360- and
365-day years, and UPDATE_climate_forcing wraps around the forcing array by
index. Feeding it a series with leap days makes the seasonal cycle drift by
one day per leap year on every wrap. Use --drop-leap-days to remove them and
--whole-years to trim any trailing partial year, giving an exactly periodic
forcing.

Typical use for the Boike/Bayelva single-point test case:

    ./txt2frognc.py \\
        --temperature Donnee/Temp_18y_day.txt \\
        --snow        Donnee/Snow_tot.txt \\
        --lat 78.9209 --lon 11.8331 \\
        --start 1998-01-01 --calendar 365_day \\
        --drop-leap-days --whole-years \\
        --outdir Input_Data/point_Bayelva

which writes:

    Input_Data/point_Bayelva/forcing_tas.nc          (variable: tas)
    Input_Data/point_Bayelva/forcing_snowthick.nc    (variable: snowthick)
    Input_Data/point_Bayelva/file_typology.nc        (grid only, no data var)

Requires: numpy, netCDF4.
"""

import argparse
import datetime as _dt
import os
import sys

import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    sys.exit("error: netCDF4 is required (pip install netCDF4)")


# --------------------------------------------------------------------------
# Variable metadata. Keys are the short names written into the netCDF file
# and referenced from FROG's namelist (name_tas_variable, etc.).
# --------------------------------------------------------------------------
VARDEFS = {
    "tas": dict(
        standard_name="air_temperature",
        long_name="Near-Surface Air Temperature",
        units="degC",
    ),
    "snowthick": dict(
        standard_name="surface_snow_thickness",
        long_name="Snow Depth",
        units="m",
    ),
    "swe": dict(
        standard_name="lwe_thickness_of_surface_snow_amount",
        long_name="Snow Water Equivalent",
        units="m",
    ),
    "tsl": dict(
        standard_name="soil_temperature",
        long_name="Soil Temperature",
        units="degC",
    ),
}

FILL = 1.0e20


def read_column(path):
    """Read a one-value-per-line text file into a 1-D float array.

    Tolerates CRLF endings, UTF-8 BOMs, blank lines and trailing whitespace,
    all of which occur in the files under Donnee/.
    """
    values = []
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for lineno, raw in enumerate(fh, start=1):
            token = raw.strip()
            if not token:
                continue
            try:
                values.append(float(token))
            except ValueError:
                sys.exit(f"error: {path}:{lineno}: cannot parse {token!r} as a number")
    if not values:
        sys.exit(f"error: {path}: no numeric data found")
    return np.asarray(values, dtype=np.float32)


def leap_filter(ntime, start_date, drop_leap, whole_years, dpy):
    """Return the indices to keep, plus a human-readable description.

    Input series are assumed to be daily and contiguous from start_date on a
    real (proleptic Gregorian) calendar. Leap days are removed first, then
    any trailing partial year is trimmed, so that the result is an exact
    multiple of dpy and wraps cleanly.
    """
    dates = [start_date + _dt.timedelta(days=i) for i in range(ntime)]
    keep = list(range(ntime))
    notes = []

    if drop_leap:
        leaps = [i for i in keep if dates[i].month == 2 and dates[i].day == 29]
        keep = [i for i in keep if i not in set(leaps)]
        if leaps:
            listed = ", ".join(str(dates[i]) for i in leaps)
            notes.append(f"removed {len(leaps)} leap day(s): {listed}")
        else:
            notes.append("no leap days found in the series")

    if whole_years:
        nfull = (len(keep) // dpy) * dpy
        if nfull == 0:
            sys.exit(f"error: fewer than {dpy} steps remain -- cannot form a whole year")
        if nfull != len(keep):
            trimmed = len(keep) - nfull
            last = dates[keep[nfull - 1]]
            notes.append(f"trimmed {trimmed} trailing step(s), series now ends {last}")
            keep = keep[:nfull]

    return keep, dates, notes


def days_per_year(calendar):
    return 360 if calendar == "360_day" else 365


def make_dataset(path, lat, lon, ntime, time_units, calendar, title, warnings=()):
    """Create a NETCDF3_CLASSIC file with the lat/lon/time skeleton FROG wants.

    NETCDF3_CLASSIC is used deliberately: it matches the existing
    Boikedata_lonlat.nc, and avoids requiring HDF5 support in the netCDF
    Fortran library FROG is linked against.
    """
    ds = nc.Dataset(path, "w", format="NETCDF3_CLASSIC")

    ds.createDimension("lat", 1)
    ds.createDimension("lon", 1)
    ds.createDimension("time", None)  # UNLIMITED -- required by get_clim_forcing

    v = ds.createVariable("lat", "f4", ("lat",))
    v.standard_name = "latitude"
    v.long_name = "Latitude"
    v.axis = "Y"
    v.units = "degrees_north"
    v[:] = np.float32(lat)

    v = ds.createVariable("lon", "f4", ("lon",))
    v.standard_name = "longitude"
    v.long_name = "Longitude"
    v.axis = "X"
    v.units = "degrees_east"
    v[:] = np.float32(lon)

    v = ds.createVariable("time", "f8", ("time",))
    v.standard_name = "time"
    v.long_name = "time"
    v.axis = "T"
    v.units = time_units
    v.calendar = calendar
    v[:] = np.arange(1, ntime + 1, dtype=np.float64)

    ds.Conventions = "CF-1.8"
    ds.title = title
    ds.history = "{}: created by txt2frognc.py".format(
        _dt.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    )
    for i, text in enumerate(warnings):
        setattr(ds, "WARNING" if i == 0 else f"WARNING_{i + 1}", text)
    return ds


def write_forcing(path, name, data, lat, lon, time_units, calendar, title,
                  warnings=()):
    ds = make_dataset(path, lat, lon, len(data), time_units, calendar, title,
                      warnings)

    meta = VARDEFS.get(name)
    if meta is None:
        # Unknown short name: still write it, but without CF metadata we
        # cannot invent. Warn so the user can add it to VARDEFS.
        print(f"  warning: no CF metadata known for variable {name!r}", file=sys.stderr)
        meta = dict(long_name=name, units="1")

    v = ds.createVariable(name, "f4", ("time", "lat", "lon"), fill_value=np.float32(FILL))
    for attr, value in meta.items():
        setattr(v, attr, value)
    v.missing_value = np.float32(FILL)
    v[:, 0, 0] = data

    ds.close()
    print(f"  wrote {path}  ({name}, {len(data)} steps)")


def write_typology(path, lat, lon, time_units, calendar, title, warnings=()):
    """Grid-definition file: same skeleton, one time step, no data variable.

    Mirrors file_typology-OnePointBoike.nc, which was produced with
    `ncks -d time,0,0 -x -v tas`.
    """
    ds = make_dataset(path, lat, lon, 1, time_units, calendar, title, warnings)
    ds.close()
    print(f"  wrote {path}  (grid definition only)")


def main():
    p = argparse.ArgumentParser(
        description="Convert FROG text forcing files to netCDF.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Typical use")[1].rsplit("Requires:", 1)[0].strip(),
    )
    p.add_argument("--temperature", metavar="FILE",
                   help="text file of near-surface air temperature [degC]")
    p.add_argument("--snow", metavar="FILE",
                   help="text file of snow thickness [m]")
    p.add_argument("--swe", metavar="FILE",
                   help="text file of snow water equivalent [m]")
    p.add_argument("--extra", metavar="NAME=FILE", action="append", default=[],
                   help="any further variable, e.g. --extra tsl=Donnee/Temp_soil.txt")
    p.add_argument("--lat", type=float, required=True, help="latitude [degrees_north]")
    p.add_argument("--lon", type=float, required=True, help="longitude [degrees_east]")
    p.add_argument("--start", default="1901-01-01", metavar="YYYY-MM-DD",
                   help="reference date for the time axis (default: 1901-01-01)")
    p.add_argument("--calendar", default="365_day", choices=["365_day", "360_day"],
                   help="must match YearType in the FROG namelist (default: 365_day)")
    p.add_argument("--drop-leap-days", action="store_true",
                   help="remove 29 February from the series (input assumed daily "
                        "and contiguous from --start on a real calendar)")
    p.add_argument("--whole-years", action="store_true",
                   help="trim any trailing partial year so the series is an exact "
                        "multiple of the calendar year length")
    p.add_argument("--warning", metavar="TEXT", action="append", default=[],
                   help="add a caveat as a global attribute; repeatable")
    p.add_argument("--outdir", default=".", help="output directory (created if needed)")
    p.add_argument("--title", default="FROG forcing, single point",
                   help="value of the netCDF global attribute 'title'")
    p.add_argument("--no-typology", action="store_true",
                   help="skip writing file_typology.nc")
    args = p.parse_args()

    # Collect requested variables -----------------------------------------
    jobs = []
    if args.temperature:
        jobs.append(("tas", args.temperature, "forcing_tas.nc"))
    if args.snow:
        jobs.append(("snowthick", args.snow, "forcing_snowthick.nc"))
    if args.swe:
        jobs.append(("swe", args.swe, "forcing_swe.nc"))
    for spec in args.extra:
        if "=" not in spec:
            sys.exit(f"error: --extra expects NAME=FILE, got {spec!r}")
        name, path = spec.split("=", 1)
        jobs.append((name, path, f"forcing_{name}.nc"))

    if not jobs:
        sys.exit("error: nothing to do -- give at least one of "
                 "--temperature / --snow / --swe / --extra")

    # Read and cross-check -------------------------------------------------
    series = {}
    for name, path, _ in jobs:
        if not os.path.exists(path):
            sys.exit(f"error: {path}: no such file")
        series[name] = read_column(path)

    lengths = {name: len(a) for name, a in series.items()}
    if len(set(lengths.values())) > 1:
        print("error: input files have different lengths:", file=sys.stderr)
        for name, n in lengths.items():
            print(f"         {name}: {n}", file=sys.stderr)
        sys.exit("       forcing variables must share a common time axis")

    ntime = next(iter(lengths.values()))
    dpy = days_per_year(args.calendar)

    # Calendar conditioning ------------------------------------------------
    try:
        start_date = _dt.date.fromisoformat(args.start)
    except ValueError:
        sys.exit(f"error: --start {args.start!r} is not a valid YYYY-MM-DD date")

    if args.drop_leap_days or args.whole_years:
        keep, dates, notes = leap_filter(ntime, start_date, args.drop_leap_days,
                                         args.whole_years, dpy)
        for note in notes:
            print(f"  {note}")
        if len(keep) != ntime:
            idx = np.asarray(keep, dtype=int)
            series = {name: arr[idx] for name, arr in series.items()}
            print(f"  {ntime} -> {len(keep)} steps "
                  f"({dates[keep[0]]} to {dates[keep[-1]]})")
            ntime = len(keep)
    nyears = ntime / dpy
    if abs(nyears - round(nyears)) > 1e-6:
        print(f"warning: {ntime} steps is {nyears:.3f} years on a {dpy}-day "
              f"calendar -- not a whole number of years. FROG wraps around the "
              f"forcing when a run outlasts it, so a partial year shifts the "
              f"seasonal cycle at every wrap. Consider --drop-leap-days "
              f"--whole-years.", file=sys.stderr)

    # Sanity checks on physical ranges ------------------------------------
    if "tas" in series:
        t = series["tas"]
        if t.min() > 100.0:
            print("warning: temperature looks like Kelvin, not degC. FROG's "
                  "fix_Kelvin_or_Celsius will convert it, but the units "
                  "attribute written here says degC.", file=sys.stderr)
    for name in ("snowthick", "swe"):
        if name in series and series[name].min() < 0.0:
            print(f"warning: {name} contains negative values "
                  f"(min = {series[name].min():.4g})", file=sys.stderr)

    # Write ----------------------------------------------------------------
    os.makedirs(args.outdir, exist_ok=True)
    time_units = f"days since {args.start} 00:00:00"

    print(f"{ntime} time steps, {nyears:.2f} years on a {dpy}-day calendar")
    print(f"point at {args.lat:.4f} N, {args.lon:.4f} E")

    caveats = list(args.warning)
    if args.drop_leap_days:
        caveats.append(
            "Leap days (29 February) were removed from the original record to "
            "produce an exactly periodic series on a {}-day calendar. "
            "Do not compare this time axis against real dates without "
            "accounting for that.".format(dpy))

    for name, _, fname in jobs:
        write_forcing(os.path.join(args.outdir, fname), name, series[name],
                      args.lat, args.lon, time_units, args.calendar, args.title,
                      caveats)

    if not args.no_typology:
        write_typology(os.path.join(args.outdir, "file_typology.nc"),
                       args.lat, args.lon, time_units, args.calendar,
                       args.title + " (grid definition)", caveats)

    # Namelist snippet -----------------------------------------------------
    print("\nCorresponding namelist entries:\n")
    print("&inputFiles")
    for name, _, fname in jobs:
        if name == "tas":
            print(f'  forc_tas_file       = "{os.path.join(args.outdir, fname)}"')
            print(f'  name_tas_variable   = "tas"')
    print("  ...")
    print("/\n")
    print("&inputsGrid")
    print(f'  typology_file = "{os.path.join(args.outdir, "file_typology.nc")}"')
    print(f'  mask_file     = "{os.path.join(args.outdir, jobs[0][2])}"')
    print("  ...")
    print("/")

    snow_jobs = [j for j in jobs if j[0] == "snowthick"]
    if snow_jobs:
        print("\nNote: as of commit 9919138, FROG does not read snow forcing from "
              "file in offline mode -- spatialvars_init sets "
              "forcing_snow_thick(:,:) = 1.5 unconditionally under "
              "#if ( SNOW_EFFECT == 1 ). The snow file written here is only "
              "usable once that path is implemented.", file=sys.stderr)


if __name__ == "__main__":
    main()
