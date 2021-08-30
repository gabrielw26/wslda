"""Test IO Routines"""
import glob
import os.path
import stat
import tempfile

import numpy as np

from zope.interface.verify import verifyClass, verifyObject
import pytest

from wdata import io


@pytest.fixture
def data_dir():
    with tempfile.TemporaryDirectory() as data_dir:
        yield data_dir


@pytest.fixture(params=["wdat", "npy"])
def ext(request):
    yield request.param


@pytest.fixture(params=[1, 2, 3])
def dim(request):
    yield request.param


@pytest.fixture
def infofile(data_dir, ext, dim):
    """Reasonable datasets for testing."""
    Nt = 4
    Nxyz = (4, 8, 16)
    dxyz = (0.1, 0.2, 0.3)

    variables = [
        io.Var(density=np.random.random((Nt,) + Nxyz[:dim])),
        io.Var(current1=np.random.random((Nt, 1) + Nxyz[:dim])),
        io.Var(current2=np.random.random((Nt, 2) + Nxyz[:dim])),
        io.Var(current3=np.random.random((Nt, 3) + Nxyz[:dim])),
    ]

    data = io.WData(
        prefix="tmp",
        data_dir=data_dir,
        ext=ext,
        dim=dim,
        Nxyz=Nxyz,
        dxyz=dxyz,
        variables=variables,
        Nt=Nt,
        aliases={"n": "density"},
        constants=dict(hbar=1.23),
    )

    data.save()
    infofile = data.infofile
    del data

    yield infofile


class TestIO:
    def test_interfaces(self, data_dir):
        assert verifyClass(io.IVar, io.Var)
        assert verifyClass(io.IWData, io.WData)

        Nt = 1
        Nxyz = (4, 8, 16)
        dxyz = (0.1, 0.2, 0.3)

        var = io.Var(density=np.random.random((Nt,) + Nxyz))
        data = io.WData(
            prefix="tmp",
            data_dir=data_dir,
            Nxyz=Nxyz,
            dxyz=dxyz,
            variables=[var],
            Nt=Nt,
        )

        assert verifyObject(io.IVar, var)
        assert verifyObject(io.IWData, data)

    def test_wdata1(self, data_dir, ext):
        Nxyz = (4, 8, 16)
        dxyz = (0.1, 0.2, 0.3)
        prefix = "tmp"
        full_prefix = os.path.join(data_dir, prefix)

        data = io.WData(
            prefix=prefix, data_dir=data_dir, Nxyz=Nxyz, dxyz=dxyz, Nt=1, ext=ext
        )

        xyz = data.xyz

        np.random.seed(2)

        psi = np.random.random(Nxyz + (2,)) - 0.5
        psi = psi.view(dtype=complex)[..., 0]
        assert psi.shape == Nxyz

        psis = [psi]
        densities = [abs(psi) ** 2]

        data = io.WData(
            data_dir=data_dir,
            xyz=xyz,
            variables=[io.Var(density=densities), io.Var(delta=psis)],
            ext=ext,
        )
        data.save()

        res = io.WData.load(full_prefix=full_prefix)
        assert np.allclose(res.delta, psi)
        assert np.allclose(res.density, densities)

    def test_wdata_backward_compatible(self, data_dir):
        """Test loading from wdata using description in the docs."""
        prefix = "test"
        full_prefix = os.path.join(data_dir, f"{prefix}")
        infofile = f"{full_prefix}.wtxt"
        with open(infofile, "w") as f:
            f.write(
                """
# Comments with additional  info about data set
# Comments are ignored when reading by parser

NX       24 # lattice
NY       28 # lattice
NZ       32 # lattice
DX        1 # spacing
DY        1 # spacing
DZ        1 # spacing
datadim   3 # dimension of block size: 1=NX, 2=NX*NY, 3=NX*NY*NZ
prefix test # prefix for files ... files have names prefix_variable.format
cycles   10 # number of cycles (measurements)
t0        0 # time value for the first cycle
dt        1 # time interval between cycles

# variables
# tag           name     type      unit      format
var        density_a     real      none        wdat
var            delta  complex      none        wdat
var        current_a   vector      none        wdat

# links
# tag           name    link-to
link       density_b  density_a
link       current_b  current_a

# consts
# tag       name       value
const         eF         0.5
const         kF           1
"""
            )

        Nxyz = (24, 28, 32)
        dxyz = (1, 1, 1)
        Nt = cycles = 10
        t0 = 0
        dt = 1
        xyz = np.meshgrid(
            *((np.arange(_N) - _N / 2) * _dx for _N, _dx in zip(Nxyz, dxyz)),
            indexing="ij",
            sparse=True,
        )
        ts = np.arange(Nt) * dt + t0

        # Construct some data
        x, y, z = xyz
        Lxyz = np.multiply(Nxyz, dxyz)
        wx, wy, wz = 2 * np.pi / Lxyz
        cos, sin = np.cos, np.sin
        density = np.array(
            [cos(wx * t * x) * cos(wy * t * y) * cos(wz * t * z) for t in ts]
        )
        delta = 1 + 1j * density
        gradient = [
            [
                -wx * t * sin(wx * t * x) * cos(wy * t * y) * cos(wz * t * z),
                -wy * t * cos(wx * t * x) * sin(wy * t * y) * cos(wz * t * z),
                -wz * t * cos(wx * t * x) * cos(wy * t * y) * sin(wz * t * z),
            ]
            for t in ts
        ]
        vars = dict(density_a=density, delta=delta, current_a=gradient)
        for var in vars:
            with open(f"{full_prefix}_{var}.wdat", "wb") as f:
                f.write(np.ascontiguousarray(vars[var]).tobytes())

        wdata = io.WData.load(infofile=infofile)
        assert wdata.prefix == "test"
        assert wdata.description == "\n".join(
            [
                "Comments with additional  info about data set",
                "Comments are ignored when reading by parser",
            ]
        )
        assert wdata.data_dir == data_dir
        assert wdata.ext == "wdat"
        assert dict(wdata.aliases) == {
            "density_b": "density_a",
            "current_b": "current_a",
        }
        assert dict(wdata.constants) == {"eF": 0.5, "kF": 1}

        assert all(np.allclose(_x, __x) for (_x, __x) in zip(wdata.xyz, xyz))
        assert wdata.dim == 3
        assert wdata.Nxyz == (24, 28, 32)
        assert wdata.xyz0 == (-24 / 2, -28 / 2, -32 / 2)
        assert wdata.dxyz == (1, 1, 1)
        assert wdata.eF == 0.5
        assert wdata.kF == 1.0
        assert sorted(wdata) == [
            "current_a",
            "current_b",
            "delta",
            "density_a",
            "density_b",
            "eF",
            "kF",
        ]

        assert np.allclose(wdata.t, ts)
        assert np.allclose(wdata.t0, 0)
        assert np.allclose(wdata.dt, 1)

        density_a, delta_, current_a = wdata.variables
        assert wdata.density_a is density_a.data
        assert density_a.name == "density_a"
        assert density_a.description == ""
        assert density_a.filename.endswith(".wdat")
        assert density_a.unit == "none"
        assert density_a.filename == f"{full_prefix}_density_a.wdat"
        assert density_a.descr == "<f8"
        assert not density_a.vector
        assert np.allclose(density_a.data, density)

        assert delta_.name == "delta"
        assert delta_.description == ""
        assert delta_.filename.endswith(".wdat")
        assert delta_.unit == "none"
        assert delta_.filename == f"{full_prefix}_delta.wdat"
        assert delta_.descr == "<c16"
        assert not delta_.vector
        assert np.allclose(delta.data, delta)

        assert current_a.name == "current_a"
        assert current_a.description == ""
        assert current_a.filename.endswith(".wdat")
        assert current_a.unit == "none"
        assert current_a.filename == f"{full_prefix}_current_a.wdat"
        assert current_a.descr == "<f8"
        assert current_a.vector
        assert np.allclose(current_a.data, gradient)

    ######################################################################
    # Regression tests
    def test_issue_2(self, data_dir):
        """Test datadim=2 issue."""
        prefix = "test"
        full_prefix = os.path.join(data_dir, f"{prefix}")
        infofile = f"{full_prefix}.wtxt"
        with open(infofile, "w") as f:
            f.write(
                """
# Generated by td-wslda-2d [12/13/20-16:51:30]

NX       24 # lattice
NY       28 # lattice
NZ       32 # lattice
DX        1 # spacing
DY        1 # spacing
DZ        1 # spacing
datadim   2 # dimension of block size: 1=NX, 2=NX*NY, 3=NX*NY*NZ
prefix test # prefix for files ... files have names prefix_variable.format
cycles    3 # number of cycles (measurements)
t0        0 # time value for the first cycle
dt        1 # time interval between cycles

# variables
# tag           name     type      unit      format
var        density       real      none        wdat
var        current2      vector    none        wdat
var        current3      vector    none        wdat
"""
            )

        Nxyz = (24, 28)
        dxyz = (1, 1)
        Nt = cycles = 3
        t0 = 0
        dt = 1
        xyz = np.meshgrid(
            *((np.arange(_N) - _N / 2) * _dx for _N, _dx in zip(Nxyz, dxyz)),
            indexing="ij",
            sparse=True,
        )
        ts = np.arange(Nt) * dt + t0

        # Construct some data
        x, y = xyz
        Lxyz = np.multiply(Nxyz, dxyz)
        wx, wy = 2 * np.pi / Lxyz
        cos, sin = np.cos, np.sin
        density = np.array([cos(wx * t * x) * cos(wy * t * y) for t in ts])
        current2 = np.array(
            [
                [cos(wx * t * x) * cos(wx * t * y), sin(wx * t * x) * sin(wx * t * y)]
                for t in ts
            ]
        )
        current3 = np.array(
            [
                [
                    cos(wx * t * x) * cos(wx * t * y),
                    sin(wx * t * x) * sin(wx * t * y),
                    0 * t * x * y,
                ]
                for t in ts
            ]
        )
        vars = dict(density=density, current2=current2, current3=current3)
        for var in vars:
            with open(f"{full_prefix}_{var}.wdat", "wb") as f:
                f.write(np.ascontiguousarray(vars[var]).tobytes())

        wdata = io.WData.load(infofile=infofile)
        assert wdata.prefix == "test"
        assert wdata.description == "Generated by td-wslda-2d [12/13/20-16:51:30]"
        assert wdata.data_dir == data_dir
        assert wdata.ext == "wdat"

        assert all(np.allclose(_x, __x) for (_x, __x) in zip(wdata.xyz, xyz))
        assert wdata.dim == 2
        assert wdata.Nxyz[: wdata.dim] == (24, 28)
        assert wdata.xyz0[: wdata.dim] == (-24 / 2, -28 / 2)
        assert wdata.dxyz[: wdata.dim] == (1, 1)

        assert np.allclose(wdata.t, ts)
        assert np.allclose(wdata.t0, 0)
        assert np.allclose(wdata.dt, 1)

        (density_, current2_, current3_) = wdata.variables

        assert density_.name == "density"
        assert density_.description == ""
        assert density_.filename.endswith(".wdat")
        assert density_.unit == "none"
        assert density_.filename == f"{full_prefix}_density.wdat"
        assert density_.descr == "<f8"
        assert not density_.vector
        assert np.allclose(density_.data, density)

        assert current2_.name == "current2"
        assert current2_.description == ""
        assert current2_.filename.endswith(".wdat")
        assert current2_.unit == "none"
        assert current2_.filename == f"{full_prefix}_current2.wdat"
        assert current2_.descr == "<f8"
        assert current2_.vector
        assert np.allclose(current2_.data, current2)

        assert current3_.name == "current3"
        assert current3_.description == ""
        assert current3_.filename.endswith(".wdat")
        assert current3_.unit == "none"
        assert current3_.filename == f"{full_prefix}_current3.wdat"
        assert current3_.descr == "<f8"
        assert current3_.vector
        assert np.allclose(current3_.data, current3)

    @pytest.mark.filterwarnings("error")
    def test_issue5(self, data_dir, ext):
        x = np.array([1, 2, 4, 5])
        y = np.array([1, 2, 3, 4, 5])
        xyz = [x[:, np.newaxis], y[np.newaxis, :]]
        t = [0]
        prefix = "tmp"
        full_prefix = os.path.join(data_dir, prefix)

        data = io.WData(prefix=prefix, data_dir=data_dir, xyz=xyz, t=t, ext=ext)
        assert np.isnan(data.dt)
        assert np.isnan(data.dxyz[0])
        assert np.allclose(1.0, data.dxyz[1])
        data.save()

        infofile = f"{full_prefix}.wtxt"
        with open(infofile, "r") as f:
            found = set()
            for line in f.readlines():
                if line.startswith("DX"):
                    found.add("DX")
                    assert line == "DX         varying    # Spacing in X direction\n"
                if line.startswith("DY"):
                    found.add("DY")
                    assert line == "DY             1.0    #        ... Y ...\n"
                if line.startswith("dt"):
                    found.add("dt")
                    assert (
                        line == "dt         varying    # Time interval between frames\n"
                    )
            assert len(found) == 3

        wdata = io.WData.load(infofile=infofile)

        assert np.isnan(wdata.dt)
        assert np.isnan(wdata.dxyz[0])
        assert np.allclose(1.0, wdata.dxyz[1])

    def test_metadata(self, data_dir, ext):
        x = np.array([1, 2, 3, 5])
        y = np.array([1, 2, 3, 4, 5])
        xyz = [x[:, np.newaxis], y[np.newaxis, :]]
        t = [0, 1.2]
        prefix = "tmp"
        full_prefix = os.path.join(data_dir, prefix)

        Nxyz = sum(xyz).shape
        x = io.Var(x=np.ones((len(t),) + Nxyz, dtype=complex))
        data = io.WData(
            prefix=prefix,
            data_dir=data_dir,
            xyz=xyz,
            t=t,
            ext=ext,
            variables=[x],
            aliases={"y": "x"},
            constants=dict(hbar=1.23),
        )

        data.save()

        infofile = f"{full_prefix}.wtxt"
        with open(infofile, "r") as f:
            found = set()
            for line in f.readlines():
                if line.startswith("DX"):
                    found.add("DX")
                    assert line == "DX         varying    # Spacing in X direction\n"
                if line.startswith("DY"):
                    found.add("DY")
                    assert line == "DY             1.0    #        ... Y ...\n"
                if line.startswith("dt"):
                    found.add("dt")
                    assert (
                        line == "dt             1.2    # Time interval between frames\n"
                    )
                if line.startswith("var"):
                    found.add("var")
                    assert (
                        line
                        == f"var         x    complex    none      {ext:>4}    # x\n"
                    )
                if line.startswith("link"):
                    found.add("link")
                    assert line == "link        y          x\n"
                if line.startswith("const"):
                    found.add("const")
                    assert line == "const    hbar     1.23"
            assert len(found) == 6

    def test_missing_prefix(self, data_dir):
        prefix = "test"
        full_prefix = os.path.join(data_dir, f"{prefix}")
        info_contents = """
# Generated by td-wslda-2d [12/13/20-16:51:30]

NX        4 # lattice
NY        5 # lattice
DX        1 # spacing
DY        1 # spacing
datadim   2 # dimension of block size: 1=NX, 2=NX*NY, 3=NX*NY*NZ
cycles    3 # number of cycles (measurements)
t0        0 # time value for the first cycle
dt        1 # time interval between cycles

# variables
# tag           name     type      unit      format
var        density       real      none         npy
"""
        density = np.ones((3, 4, 5), dtype=float)
        filename = f"{full_prefix}_density.npy"
        np.save(filename, density)

        for info_ext in ["wtxt", "custom"]:
            infofile = f"{full_prefix}.{info_ext}"
            with open(infofile, "w") as f:
                f.write(info_contents)
            data = io.WData.load(infofile)
            assert data.prefix == prefix

    def test_save(self, infofile):
        data = io.WData.load(infofile)
        data_dir = os.path.dirname(infofile)
        new_data_dir = os.path.join(data_dir, "new")
        data.data_dir = new_data_dir
        data.save(force=True)

    def test_dts(self):
        data = io.WData(Nxyz=(4, 5, 6), t=[1, 2, 3, 4, 5])
        assert data.dt == 1
        data = io.WData(Nxyz=(4, 5, 6), t=[1, 2, 3, 5])
        assert np.isnan(data.dt)

    def test__get_ext(self, ext):
        x = io.Var(x=np.zeros((1, 4, 5)), filename=f"x.{ext}")

        data = io.WData(Nxyz=x.shape[1:], variables=[x])
        info = data.get_metadata()
        assert ext in info

    def test_readonly(self, data_dir, ext):
        """Issue 13: fail on read-only file systems."""
        Nxyz = (4, 8)
        dxyz = (0.1, 0.2)
        prefix = "tmp"
        full_prefix = os.path.join(data_dir, prefix)

        data = io.WData(
            prefix=prefix, data_dir=data_dir, Nxyz=Nxyz, dxyz=dxyz, Nt=1, ext=ext
        )

        xyz = data.xyz

        np.random.seed(2)

        psi = np.random.random(Nxyz + (2,)) - 0.5
        psi = psi.view(dtype=complex)[..., 0]
        assert psi.shape == Nxyz

        psis = [psi]
        densities = [abs(psi) ** 2]

        data = io.WData(
            data_dir=data_dir,
            xyz=xyz,
            variables=[io.Var(density=densities), io.Var(delta=psis)],
            ext=ext,
        )
        data.save()

        # Make files read-only
        for (dirpath, dirnames, filenames) in os.walk(data_dir):
            for file in filenames:
                os.chmod(os.path.join(dirpath, file), stat.S_IRUSR)

        # Try loading
        res = io.WData.load(full_prefix=full_prefix)
        assert np.allclose(res.delta, psi)
        assert np.allclose(res.density, densities)


class TestVar:
    def test_descr(self):
        """Test setting descr."""
        v = io.Var(x=np.ones(10), descr=complex)
        assert v.data.dtype == complex
        assert io.WData._get_type(v) == "abscissa"

        v = io.Var(x=np.ones((3, 10), dtype="<f4"))
        assert io.WData._get_type(v) == "<f4"
        assert io.WData._get_descr("<f4") == "<f4"

    def test_set(self):
        """Test setting attributes."""
        v = io.Var(x=np.ones(10))
        assert v.data.dtype == float
        assert v.shape == (10,)

        v.data = 1j * np.arange(12)
        assert v.data.dtype == complex
        assert v.shape == (12,)

        v.shape = (3, 4)
        assert v.data.shape == (3, 4)

        with pytest.raises(ValueError) as excinfo:
            v.shape = (3, 5)
        assert (
            str(excinfo.value)
            == "Property shape=(3, 5) incompatible with data.shape=(3, 4)"
        )

    def test_write_data(self, data_dir):
        filename = os.path.join(data_dir, "tmp.wdat")

        v = io.Var(x=np.arange(10))
        with pytest.raises(ValueError) as excinfo:
            v.write_data()
        assert str(excinfo.value) == "No filename specified in Var."

        v.write_data(filename)
        assert os.path.exists(f"{filename}")

        with pytest.raises(IOError) as excinfo:
            v.write_data(filename)
        assert str(excinfo.value) == f"File '{filename}' already exists!"

        with pytest.raises(NotImplementedError) as excinfo:
            v.write_data(filename="x.unknown")
        assert str(excinfo.value) == "Unsupported extension for 'x.unknown'"

        v = io.Var(name="x", filename=filename)
        with pytest.raises(ValueError) as excinfo:
            v.write_data()
        assert str(excinfo.value) == "Missing data for 'x'!"

    def test_getattr(self):
        data = io.WData(Nxyz=(4, 5), t=[1, 2])
        with pytest.raises(AttributeError) as excinfo:
            data.unknown
        assert str(excinfo.value) == "'WData' object has no attribute 'unknown'"


class TestErrors:
    """Test coverage and errors."""

    def test_var_errors(self):
        with pytest.raises(ValueError) as excinfo:
            io.Var()
        assert str(excinfo.value) == "Must provide `name` or data as a kwarg: got {}"

        with pytest.raises(ValueError) as excinfo:
            io.Var(data=np.ones((4, 4)))
        assert str(excinfo.value) == "Got data but no name."

        v = io.Var(name="x", filename="tmp.unknown")
        with pytest.raises(NotImplementedError) as excinfo:
            v.load_data()
        assert str(excinfo.value) == "Data format of 'tmp.unknown' not supported."

    def test_wdata_errors(self):
        with pytest.raises(ValueError) as excinfo:
            io.WData()
        assert str(excinfo.value) == "Must provide one of xyz or Nxyz"

        with pytest.raises(ValueError) as excinfo:
            io.WData(Nxyz=(4, 5))
        assert str(excinfo.value) == "Must provide t, Nt, or a variable with data"

        with pytest.raises(ValueError) as excinfo:
            io.WData(Nxyz=(3, 4, 5), Nt=1)
        assert str(excinfo.value) == "First dimension of Nxyz==(3, 4, 5) must be > 3."

        with pytest.raises(ValueError) as excinfo:
            x = io.Var(x=np.zeros((2, 4, 5)))
            io.WData(Nxyz=(4, 5), Nt=1, variables=[x])
        assert (
            str(excinfo.value)
            == "Variable 'x' has incompatible Nt=1: data.shape[0] = 2"
        )

        with pytest.raises(ValueError) as excinfo:
            x = io.Var(x=np.zeros((2, 4)))
            io.WData(Nxyz=(5,), Nt=2, variables=[x])
        assert (
            str(excinfo.value)
            == "Variable 'x' has incompatible Nxyz=(5,): data.shape[-1:] = (4,)"
        )

        with pytest.raises(ValueError) as excinfo:
            x = io.Var(x=np.zeros((2, 4, 4, 6)))
            io.WData(Nxyz=(4, 6), dim=2, Nt=2, variables=[x])
        assert (
            str(excinfo.value)
            == "Variable 'x' has incompatible dim=2: data.shape = (2, 4, 4, 6)"
        )

    @pytest.mark.filterwarnings("error")
    def test_getitem(self):
        x = io.Var(x=np.zeros((2, 4, 5)))
        data = io.WData(Nxyz=(4, 5), t=[1, 2], variables=[x], constants=dict(x=3))
        with pytest.warns(
            UserWarning, match="Variable x hides constant of the same name"
        ):
            data.x

        data = io.WData(Nxyz=(4, 5), t=[1, 2], variables=[x], constants=dict(y=3))
        with pytest.raises(AttributeError) as excinfo:
            data.q
        assert str(excinfo.value) == "'WData' object has no attribute 'q'"

        data.keys()

    def test__dir__(self):
        x = io.Var(x=np.zeros((2, 4, 5)))
        data = io.WData(Nxyz=(4, 5), t=[1, 2], variables=[x], constants=dict(x=3))
        assert dir(data) == ["x"]

    def test_load_errors(self, data_dir):
        with pytest.raises(ValueError) as excinfo:
            io.WData.load()
        assert str(excinfo.value) == "Must provide either infofile or full_prefix."

        with pytest.raises(ValueError) as excinfo:
            io.WData.load(infofile="info", full_prefix="/")
        assert str(excinfo.value) == "Got both infofile=info and full_prefix=/."

        with pytest.raises(NotImplementedError) as excinfo:
            io.WData.load(full_prefix=data_dir)

    def test_save_errors(self, infofile):
        with pytest.raises(NotImplementedError) as excinfo:
            data = io.load_wdata(infofile)

        data = io.load_wdata(infofile=infofile)

        with pytest.raises(IOError) as excinfo:
            data.save()
        assert str(excinfo.value) == f"File '{infofile}' already exists!"

        data_dir = os.path.dirname(infofile)
        new_data_dir = os.path.join(data_dir, "new")
        data.data_dir = new_data_dir
        with pytest.raises(IOError) as excinfo:
            data.save()
        assert (
            str(excinfo.value) == f"Directory data_dir={new_data_dir} does not exist."
        )

    def test_invalid_vector(self, infofile):
        data = io.WData.load(infofile)

        if data.dim == 1:
            for var in data.variables:
                if var.name == "current1":
                    break

            if var.filename.endswith(".wdat"):
                # Break data by changing size... make sure we copy it first though since
                # we are using mem-mapped files!
                A = np.array(var.data.ravel()[:-1])
                filename = var.filename
                assert os.path.exists(filename)
                with open(filename, "wb") as fd:
                    fd.write(A.tobytes())

                _data = np.memmap(filename, dtype=np.dtype(A.dtype))
                assert _data.shape == A.shape

                # This should work...
                data_ = io.WData.load(infofile, check_data=False)

                # but this should fail
                with pytest.raises(ValueError) as excinfo:
                    data_.current1
                msg = str(excinfo.value)
                assert "Shape of data" in msg
                assert "Nv=0.9375 must be an integer." in msg
