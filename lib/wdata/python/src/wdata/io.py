"""Input and output of states to various data formats.

Binary data is stored in one of the following formats denoted by their
extension:

wdata:
   Filename ends in ``.wdat``.  This is a simple binary format
   contiguous in C ordering.  In this format there is no information
   about the size, datatype, or shape of the data.  Thus, a
   corresponding metadata file is needed.
npy:
   Filename ends in ``.wdat``  This is the NPY format.  This includes
   the size and datatype, so one can infer the type of data without a
   metadata file.  We assume that the shape is one of:

   * ``(Nt,) + Nxyz``: Scalar data with `Nt` blocks each of shape ``Nxyz``.
   * ``(Nt, Nv) + Nxyz``: Vector data with ``Nt`` blocks, and ``Nv <= 3``
     components each.  Note: we assume ``Nx > 3`` so that these two cases
     can be distinguished by looking at ``shape[1]``.

Note: Both of these formats are simple binary representations, and hence can be loaded
as memory mapped files.  We do this by default in the current code so users can extract
slices as needed without loading the entire dataset.  Be aware that the data can be big,
so if you do anything that makes a copy of the array (like computing `np.sin(...)`),
this will load the entire array into memory.

See Also
--------
* https://numpy.org/doc/stable/reference/generated/numpy.memmap.html

"""
import collections.abc
from collections import OrderedDict
import os.path
from warnings import warn

from datetime import datetime
import numpy as np

from pytz import timezone
import tzlocal

from zope.interface import implementer, Attribute, Interface
from zope.interface.common.collections import IMapping


class IVar(Interface):
    """Interface for a variable descriptor.

    Instances may actually store the data, or may represent data on
    disk that has not yet been loaded.

    If data is provided, then it must have the shape:

    * Abscissa data: ``(N,)``
    * Scalar data: ``(Nt,) + Nxyz``
    * Vector data: ``(Nt, 3) + Nxyz``

    Assumes that ``Nx > 3`` so that if ``shape[1] <= 3``, the data is vector.
    """

    # Required attributes
    name = Attribute("Name of variable")
    data = Attribute("Array-like object with access to the data.")
    description = Attribute("Single-line description of the data.")
    ext = Attribute("Extension for data file format.")
    unit = Attribute("String representation of units for the data.")
    filename = Attribute("Full path of the data file.")

    # Derived attributes: can be computed from data if present.
    descr = Attribute("NumPy data descriptor dtype.descr")
    shape = Attribute("Shape of the data so that wdat files can be read.")

    def __init__(
        name=None,
        data=None,
        description=None,
        ext=None,
        unit="none",
        filename=None,
        descr=None,
        shape=None,
        **kwargs,
    ):
        """Initialize the object.

        Arguments
        ---------
        name : str
           Name of the variable.  If not provided, then the data
           should be set using kwarg syntax ``__init__(name=data)``.
        kwargs : {name: data}
           If no name is provided, then kwargs should have a single
           entry which is the data and the key is the name.
        filename : str
           Name of file (including path).

        description : str
           Single-line comment for variable.  Defaults to the name.
        ext : str
           Extension for data file format - will override default from
           IWData object if provided.
        unit : str
           Unit of data for metadata file.
        shape : tuple, None
           If no data is provided, then this is needed so that wdat
           files can be loaded.
        """

    def write_data(force=False, ext="wdat"):
        """Write the data to disk.

        Arguments
        ---------
        force : bool
           If True, overwrite existing files, otherwise raise IOError
           if file exists.
        ext : str
           File format.  Currently supported values are ``'wdat'`` (pure
           binary) and ``'npy'`` for the numpy NPY format.
        """


class IWData(IMapping):
    """Interface of a complete dataset.

    Also used to define abscissa as follows::

        x = np.arange(Nx)*dx + x0
        t = np.arange(Nt)*dt + t0

    If not provided, the default values are ``x0 = -Lx/2 = -Nx*dx/2`` so
    that the spatial abscissa are centered on the orgin and ``t0 = 0``.
    """

    prefix = Attribute("Prefix for all data files")
    description = Attribute("Description of the data.")
    data_dir = Attribute("Location of data.")
    ext = Attribute("Default extension for data file format.")

    variables = Attribute("List of IVar instances for the variables")

    aliases = Attribute("Dictionary of aliases.")
    constants = Attribute("Dictionary of constants.")

    xyz = Attribute("Abscissa (x, y, z) broadcast appropriately.")
    dim = Attribute("Dimensions. 1: Nx, 2: Nx*Ny, 3: Nx*Ny*Nz")

    Nxyz = Attribute("Shape of data (Nx, Ny, Nz) and lattice shape.")
    xyz0 = Attribute("Offsets (x0, y0, z0).")

    dxyz = Attribute("Spacing (dx, dy, dz).  `NotImplemented` if not uniform.")

    t = Attribute("Times for each frame")

    t0 = Attribute("Time of initial frame.")
    dt = Attribute("Time steps between frames.  None if not uniform.")

    infofile = Attribute("Name of infofile")

    def __init__(
        prefix="tmp",
        description="",
        data_dir=".",
        ext="wdat",
        Nxyz=None,
        dxyz=(1, 1, 1),
        xyz0=None,
        xyz=None,
        Nt=None,
        t0=0,
        dt=1,
        t=None,
        variables=None,
        aliases=None,
        constants=None,
    ):
        """Constructor.

        Arguments
        ---------
        prefix : str
           Prefix for files.  Default is ``'tmp'`` allowing this class to
           be used to generate abscissa.
        description : str
           User-friendly description of the data.  Will be inserted as
           a comment in the metadata file.
        data_dir : str
           Path to the directory containing the data.
        ext : str
           Default extension for data file format used if variables do
           not specify something else.  Also used for abscissa.

        Nxyz : (int,)*dim

        Nxyz, dxyz, xyz0 : (int,)*dim, (float,)*dim, (float,)*dim or None
           If these are provided, then the abscissa ``xyz`` are computed
           with equal spacings ``x = np.arange(Nx)*dx + x0``.  Default
           offeset (if ``xyz0 == None``) is centered ``x0 = -Lx/2 = -Nx*dx/2``.
        xyz : (array_like,)*dim
           Alternatively, the abscissa can be provided and the
           previous properties will be computed (if defined).
        Nt, t0, dt : int, float, float
           If these are provided, then the times are equally spaced.
           ``t = np.arange(Nt)*dt + t0``.
        t : array_like
           Alternatively, the times can be provided and the previous
           properties will be computed.

        variables : [IVar]
           List of Variables.
        """

    def get_metadata(header=None):
        """Return a string representation of the metadata file.

        Arguments
        ---------
        header : str
           Descriptive header to be added as a comment at the top of
           the metadata file before self.description.
        """

    def save(force=False):
        """Save data to disk.

        Arguments
        ---------
        force : bool
           If True, create needed directories and overwrite existing
           data files.  Otherwise, raise an IOError.
        """

    def load(infofile=None, full_prefix=None):
        """Load data from disk.

        Arguments
        ---------
        infofile : str
           If provided, then read this file and use the information
           contained within it to load the data.
        full_prefix : str
           Full prefix for data files, including directory if not in
           the current directory.  I.e. ``prefix='data_dir/run1'`` will
           look for data in the form of ``'data_dir/run1_*.*'`` and an
           infofile of the form ``'data_dir/run1.wtxt'``.

           The full_prefix will be split at the final path
           separation and what follows will be the ``prefix``.


        **No infofile option**

        Data can be provided without metadata if the following files
        are present::

           <full_prefix>_x.npy
           <full_prefix>_y.npy
           <full_prefix>_z.npy
           <full_prefix>_t.npy
           <full_prefix>_*.npy

        In his case, the abscissa will be determined by loading the
        first four files, and the remaining files will defined the
        variables.  No constants or links will be defined in this
        case.

        One other option is provided when ``full_prefix`` is a
        directory: i.e. ends with a path separator.  In this case, the
        directory will be assumed to contain a single set of data with
        a prefix determined by the filenames if the files above exist.
        """

    def __getattr__(key):
        """Convenience method for variable access.

        Returns the data of the named variable or the named
        constant, following aliases if defined."""


@implementer(IVar)
class Var(object):
    def __init__(
        self,
        name=None,
        data=None,
        description=None,
        filename=None,
        ext=None,
        unit="none",
        descr=None,
        shape=None,
        **kwargs,
    ):
        if name is None:
            if data is not None:
                raise ValueError("Got data but no name.")
            if not len(kwargs) == 1:
                raise ValueError(
                    f"Must provide `name` or data as a kwarg: got {kwargs}"
                )
            name, data = kwargs.popitem()
        self.name = name
        self.description = description
        self.filename = filename
        self.unit = unit
        self.ext = ext
        self._descr = descr
        self._data = data
        self.shape = shape
        self.init()

    def init(self):
        """Perform all initialization and checks.

        Called by the constructor, and before saving.
        """
        if self._data is not None:
            # Convert data and check some properties
            args = {}
            if self._descr is not None:
                args.update(astype=np.dtype(self._descr))
            self._data = np.ascontiguousarray(self._data, **args)

        if self.description is None:
            self.description = self.name

    @property
    def data(self):
        """Return or load data."""
        if self._data is None:
            self.load_data()
        return self._data

    @data.setter
    def data(self, data):
        self._data = data
        self.init()

    @property
    def descr(self):
        if self._data is None:
            return self._descr
        else:
            descr = self._data.dtype.descr
            assert len(descr) == 1
            return descr[0][-1]

    @property
    def vector(self):
        return self.shape[1] <= 3

    @property
    def shape(self):
        if self._data is None:
            return self._shape
        else:
            return self._data.shape

    @shape.setter
    def shape(self, shape):
        if self._data is not None and shape is not None and shape != self._data.shape:
            raise ValueError(
                f"Property shape={shape} incompatible "
                + f"with data.shape={self._data.shape}"
            )
        self._shape = shape

    def write_data(self, filename=None, force=False, ext="wdat"):
        """Write self.data to the specified file."""
        self.init()

        if filename is None:
            filename = self.filename

        if filename is None:
            raise ValueError("No filename specified in Var.")

        if self._data is None:
            raise IOError(f"Missing data for '{self.name}'!")

        if os.path.exists(filename) and not force:
            raise IOError("File '{filename}' already exists!")

        A = self._data
        if self.ext is not None:
            ext = self.ext

        if ext == "wdat":
            with open(".".join([filename, ext]), "wb") as fd:
                fd.write(A.tobytes())
        elif ext == "npy":
            np.save(filename, A)
        else:
            raise NotImplementedError("Unsupported ext={self.ext}")

    def load_data(self):
        """Load the data from file."""
        if self.ext == "npy":
            self._data = np.load(self.filename, mmap_mode="r")
        elif self.ext == "wdat":
            self._data = np.memmap(self.filename, dtype=np.dtype(self.descr)).reshape(
                self.shape
            )
        else:
            raise NotImplementedError(f"Data format ext={self.ext} not supported.")


@implementer(IWData)
class WData(collections.abc.Mapping):
    """Base implementation."""

    # This is the extension used for infofiles
    _infofile_extension = "wtxt"

    def __init__(
        self,
        prefix="tmp",
        description="",
        data_dir=".",
        ext="wdat",
        Nxyz=None,
        dxyz=(1, 1, 1),
        xyz0=None,
        xyz=None,
        Nt=None,
        t0=0,
        dt=1,
        t=None,
        variables=None,
        aliases=None,
        constants=None,
    ):
        self.prefix = prefix
        self.description = description
        self.data_dir = data_dir
        self.ext = ext
        self.variables = variables
        self.aliases = aliases
        self.constants = constants
        self.xyz, self.Nxyz, self.dxyz, self.xyz0 = xyz, Nxyz, dxyz, xyz0
        self.t, self.Nt, self.dt, self.t0 = t, Nt, dt, t0
        self.init()

    def init(self):
        """Perform all initialization and checks.

        Called by the constructor, and before saving.
        """
        # Abscissa
        if self.xyz is not None:
            xyz = list(map(np.ravel, self.xyz))
            Nxyz, dxyz, xyz0 = [], [], []
            for x in xyz:
                dxs = np.diff(x)
                dx = dxs.mean()
                if not np.allclose(dxs, dx):
                    dx = NotImplemented

                Nxyz.append(len(x))
                dxyz.append(dx)
                xyz0.append(x[0])
        else:
            if self.Nxyz is None:
                raise ValueError("Must provide one of xyz or Nxyz")

            Nxyz, dxyz, xyz0 = self.Nxyz, self.dxyz, self.xyz0
            _xyz0 = -np.multiply(Nxyz, dxyz) / 2
            if xyz0 is None:
                xyz0 = _xyz0
            else:
                # Allow for individual None values or NaN values.
                xyz0 = tuple(
                    _x0 if (x0 is None or x0 != x0) else x0
                    for x0, _x0 in zip(xyz0, _xyz0)
                )

            xyz = [np.arange(_N) * _dx + _x0 for _N, _dx, _x0 in zip(Nxyz, dxyz, xyz0)]

        # Make sure abscissa are appropriately broadcast.
        self.xyz = np.meshgrid(*xyz, indexing="ij", sparse=True)
        self.Nxyz, self.dxyz, self.xyz0 = map(tuple, (Nxyz, dxyz, xyz0))

        # Times
        if self.t is not None:
            t = np.ravel(self.t)
            Nt = len(t)
            dts = np.diff(t)
            dt = dts.mean()
            if not np.allclose(dts, dt):
                dt = NotImplemented
            t0 = t[0]
        else:
            Nt, dt, t0 = self.Nt, self.dt, self.t0
            if Nt is None:
                if self.variables is not None:
                    for var in self.variables:
                        if var.data is not None:
                            Nt = var.data.shape[0]
                            break
                if Nt is None:
                    raise ValueError("Must provide t, Nt, or a variable with data")

            t = np.arange(Nt) * dt + t0

        self.t, self.Nt, self.dt, self.t0 = t, Nt, dt, t0

        # Check variables
        if self.variables is not None:
            dim = self.dim
            Nt = self.Nt
            for var in self.variables:
                if var.data is not None:
                    name, data = var.name, var.data
                    if Nt != var.data.shape[0]:
                        raise ValueError(
                            f"Variable '{name}'has incompatible Nt={Nt}:"
                            + f" shape[0] = {data.shape[0]}"
                        )
                    if var.data.shape[-self.dim :] != self.Nxyz:
                        raise ValueError(
                            f"Variable '{name}' has incompatible Nxyz={Nxyz}:"
                            + f" shape[-{self.dim}:] = {data.shape[-self.dim:]}"
                        )
                    if (var.vector and len(var.shape) - 2 != dim) or (
                        not var.vector and len(var.shape) - 1 != dim
                    ):
                        raise ValueError(
                            f"Variable '{name}' has incompatible dim={dim}:"
                            + f" data.shape = {data.shape}"
                        )

    @property
    def dim(self):
        return len(self.xyz)

    def get_metadata(self, header=None):
        # Pad these with 1's for backwards compatibility with
        # Gabriel's existing code
        Nxyz = self.Nxyz + (1,) * (3 - len(self.Nxyz))
        dxyz = self.dxyz + (1,) * (3 - len(self.dxyz))
        if "None" in dxyz:
            dxyz = tuple(_dx if _dx is not None else "varying" for _dx in dxyz)

        descriptors = [
            ("NX", Nxyz[0], "Lattice size in X direction"),
            ("NY", Nxyz[1], "            ... Y ..."),
            ("NZ", Nxyz[2], "            ... Z ..."),
            ("DX", dxyz[0], "Spacing in X direction"),
            ("DY", dxyz[1], "       ... Y ..."),
            ("DZ", dxyz[2], "       ... Z ..."),
            ("prefix", self.prefix, "datafile prefix: <prefix>_<var>.<format>"),
            ("datadim", self.dim, "Block size: 1:NX, 2:NX*NY, 3:NX*NY*NZ"),
            ("cycles", self.Nt, "Number Nt of frames/cycles per dataset"),
            ("t0", self.t0, "Time value of first frame"),
            (
                "dt",
                self.dt if self.dt is not None else "varying",
                "Time interval between frames",
            ),
        ]

        # Add X0, Y0, Z0 if not default
        for x0, dx, Nx, X in zip(self.xyz0, self.dxyz, self.Nxyz, "XYZ"):
            if dx is not None and np.allclose(x0, -Nx * dx / 2):
                continue
            descriptors.append((f"{X}0", x0, f"First point in {X} lattice"))

        # Add comments here
        lines = []

        if header is not None:
            lines.extend(["# " + _v for _v in header.splitlines()])

        if self.description is not None:
            lines.extend(["# " + _v for _v in self.description.splitlines()])

        if lines:
            lines.extend([""])

        # Process text and pad for descriptors
        lines.extend(
            pad_and_justify(
                [
                    (_name, str(_value), "# " + _comment)
                    for (_name, _value, _comment) in descriptors
                ]
            )
        )

        # Add variables here
        if self.variables:
            lines.extend(["", "# variables"])
            lines.extend(
                pad_and_justify(
                    [("# tag", "name", "type", "unit", "format", "# description")]
                    + [
                        (
                            "var",
                            _v.name,
                            self._get_type(_v),
                            _v.unit,
                            _v.ext if _v.ext is not None else self.ext,
                            f"# {_v.description}",
                        )
                        for _v in self.variables
                    ]
                )
            )

        # Add aliases here
        if self.aliases:
            lines.extend(["", "# links"])
            lines.extend(
                pad_and_justify(
                    [("# tag", "name", "link-to")]
                    + [("link", _k, self.aliases[_k]) for _k in self.aliases]
                )
            )

        # Add constants
        if self.constants:
            lines.extend(["", "# consts"])
            lines.extend(
                pad_and_justify(
                    [("# tag", "name", "value")]
                    + [("const", _k, repr(self.constants[_k])) for _k in self.constants]
                )
            )
        return "\n".join(lines)

    @property
    def infofile(self):
        return os.path.join(
            self.data_dir, ".".join([self.prefix, self._infofile_extension])
        )

    def save(self, force=False):
        t1, t2 = current_time()
        metadata = self.get_metadata(header=f"Generated by wdata.io: [{t1} = {t2}]")

        data_dir = self.data_dir
        if not os.path.exists(data_dir):
            if force:
                os.makedirs(data_dir, exist_ok=True)
            else:
                raise IOError(f"Directory data_dir={data_dir} does not exist.")

        infofile = self.infofile
        if os.path.exists(infofile) and not force:
            raise IOError(f"File {infofile} already exists!")

        with open(infofile, "w") as f:
            f.write(metadata)

        variables = list(self.variables)

        if self.xyz is not None:
            for _x, _n in zip(self.xyz, "xyz"):
                variables.append(Var(**{_n: np.ravel(_x)}))

        if self.t is not None:
            variables.append(Var(t=np.ravel(self.t)))

        for var in self.variables:
            filename = os.path.join(data_dir, "{}_{}".format(self.prefix, var.name))
            var.write_data(filename=filename, force=force, ext=self.ext)

    @classmethod
    def load(cls, infofile=None, full_prefix=None):
        if infofile is not None:
            # Load from infofile
            if full_prefix is None:
                return cls.load_from_infofile(infofile=infofile)
            raise ValueError(
                f"Got both infofile={infofile} and" + f" full_prefix={full_prefix}."
            )
        elif full_prefix is None:
            raise ValueError("Must provide either infofile or full_prefix.")

        # Check if full_prefix shows an infofile.
        infofile = ".".join([full_prefix, cls._infofile_extension])
        if os.path.exists(infofile):
            return cls.load_from_infofile(infofile=infofile)

        # No infofile option.
        raise NotImplementedError()

    @classmethod
    def load_from_infofile(cls, infofile):
        """Load data from specified infofile."""
        with open(infofile, "r") as f:
            lines = f.readlines()

        data_dir = os.path.dirname(infofile)

        # Extract header - initial comments (skipping blank lines)
        header = []
        while lines:
            if not lines[0].strip():
                pass
            elif lines[0].startswith("#"):
                header.append(lines[0][1:].strip())
            else:
                break
            lines.pop(0)

        description = "\n".join(header)

        # Pairs of [line, comment] or [line]
        lines = [_l.split("#", 1) for _l in lines]

        # Pairs of ([terms], [comment]) or ([terms], []) if no comment
        lines = [([_w.strip() for _w in _l[0].split()], _l[1:]) for _l in lines if _l]
        lines = [_l for _l in lines if _l[0]]  # Skip blank lines

        parameters = OrderedDict()
        variables = []
        aliases = OrderedDict()
        constants = OrderedDict()

        for line in lines:
            terms, comments = line
            if terms[0] == "var":
                name, type, unit, ext = terms[1:]

                # Dummy shape for now.
                if type == "vector":
                    shape = (None, 3, 100)
                else:
                    shape = (None, 100)

                variables.append(
                    Var(
                        name=name,
                        unit=unit,
                        ext=ext,
                        shape=shape,
                        descr=cls._get_descr(type=type),
                        description=" ".join(comments),
                    )
                )
            elif terms[0] == "link":
                name, link = terms[1:]
                aliases[name] = link
            elif terms[0] == "const":
                name, value = terms[1:]
                constants[name] = eval(value, constants)
            else:
                name, value = terms
                parameters[name] = value

        # Process parameters
        if "prefix" not in parameters:
            prefix = os.path.basename(infofile)
            if prefix.endswith("." + cls._infofile_ext):
                prefix = prefix[: -1 + len(cls._infofile_ext)]
            else:
                prefix = prefix.rsplit(".", 1)[0]

            warn("No prefix specified in {infofile}: assuming prefix={prefix}")
            parameters["prefix"] = prefix

        prefix = parameters["prefix"]

        Nxyz = tuple(int(parameters.pop(f"N{X}", 1)) for X in "XYZ")
        dim = int(parameters.pop("datadim", len(Nxyz)))
        Nxyz = Nxyz[:dim]

        Nt = int(parameters.pop("cycles", 0))

        # Add filenames and shapes.  Defer loading until the user needs the data
        for var in variables:
            var.filename = os.path.join(data_dir, f"{prefix}_{var.name}.{var.ext}")
            if var.vector:
                var.shape = (Nt, dim) + Nxyz
            else:
                var.shape = (Nt,) + Nxyz

        args = dict(
            prefix=prefix,
            description=description,
            data_dir=data_dir,
            Nxyz=Nxyz,
            dxyz=tuple(float(parameters.pop(f"D{X}", 1.0)) for X in "XYZ")[:dim],
            xyz0=tuple(float(parameters.pop(f"{X}0", np.nan)) for X in "XYZ")[:dim],
            Nt=Nt,
            t0=float(parameters.pop("t0", 0)),
            dt=float(parameters.pop("dt", 1)),
            variables=variables,
            aliases=aliases,
            constants=constants,
        )

        wdata = cls(**args)
        return wdata

    def load_data(self, *names):
        """Load the specified data."""
        raise NotImplementedError

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            return super().__getattribute__(key)

    ######################################################################
    # Methods for the IMapping interface.  Used by collections.abc.Mapping
    def __getitem__(self, key):
        key = self.aliases.get(key, key)
        for var in self.variables:
            if var.name == key:
                if key in self.constants:
                    warn(f"Variable {key} hides constant of the same name")
                return var.data

        if key in self.constants:
            return self.constants[key]

        raise KeyError(key)

    def keys(self):
        keys = set([_var.name for _var in self.variables])
        if self.aliases:
            keys.update(self.aliases)
        if self.constants:
            keys.update(self.constants)
        return sorted(keys)

    def __iter__(self):
        return self.keys()

    def __len__(self):
        return len(self.keys())

    ######################################################################
    ####################
    # Helpers
    @staticmethod
    def _get_type(var):
        """Return the type string for backward compatibility.

        Arguments
        ---------
        var : IVar
           Variable.
        """
        if len(var.data.shape) == 1:
            return "abscissa"
        elif var.vector:
            assert var.descr == "<f8"
            return "vector"
        elif var.descr == "<c16":
            return "complex"
        elif var.descr == "<f8":
            return "real"
        else:
            return var.descr

    @staticmethod
    def _get_descr(type):
        """Return the descr from type string for backward compatibility.

        Arguments
        ---------
        type : str
           The type field in the original WDAT format.
        """
        if type in ("vector", "real"):
            return "<f8"

        if type == ("complex"):
            return "<c16"

        return type


def load_wdata(prefix=None, infofile=None):
    if infofile is not None:
        return WData.load(infofile=infofile)


######################################################################
# Utilities and helpers
def current_time(format="%Y-%m-%d %H:%M:%S %Z%z"):
    """Return the date and time as strings (utc, local)."""
    now_utc = datetime.now(timezone("UTC"))
    # Convert to local time zone
    now_local = now_utc.astimezone(tzlocal.get_localzone())
    return (now_utc.strftime(format), now_local.strftime(format))


def pad_and_justify(lines, padding="    "):
    """Return list of padded aligned strings.

    First string is left-justified, the remaining are
    right-justified except comments.
    """
    # Add padding and justify
    justify = [str.ljust] + [str.rjust] * len(lines[0])
    for _line in lines:
        for _n in range(1, len(_line)):
            if "#" in _line[_n]:
                # Left-justify comment columns
                justify[_n] = str.ljust

    lens = np.max([[len(_v) for _v in _d] for _d in lines], axis=0)
    return [
        padding.join([justify[_n](_v, lens[_n]) for _n, _v in enumerate(_line)])
        for _line in lines
    ]
