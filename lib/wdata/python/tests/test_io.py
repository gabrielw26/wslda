"""Test IO Routines"""
import os.path
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
const         kF           1"""
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

        assert np.allclose(wdata.t, ts)
        assert np.allclose(wdata.t0, 0)
        assert np.allclose(wdata.dt, 1)

        density_a, delta_, current_a = wdata.variables
        assert density_a.name == "density_a"
        assert density_a.description == ""
        assert density_a.ext == "wdat"
        assert density_a.unit == "none"
        assert density_a.filename == f"{full_prefix}_density_a.wdat"
        assert density_a.descr == "<f8"
        assert not density_a.vector
        assert np.allclose(density_a.data, density)

        assert delta_.name == "delta"
        assert delta_.description == ""
        assert delta_.ext == "wdat"
        assert delta_.unit == "none"
        assert delta_.filename == f"{full_prefix}_delta.wdat"
        assert delta_.descr == "<c16"
        assert not delta_.vector
        assert np.allclose(delta.data, delta)

        assert current_a.name == "current_a"
        assert current_a.description == ""
        assert current_a.ext == "wdat"
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
var        current       vector    none        wdat
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
        current = np.array(
            [
                [cos(wx * t * x) * cos(wx * t * y), sin(wx * t * x) * sin(wx * t * y)]
                for t in ts
            ]
        )
        vars = dict(density=density, current=current)
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
        assert wdata.Nxyz == (24, 28)
        assert wdata.xyz0 == (-24 / 2, -28 / 2)
        assert wdata.dxyz == (1, 1)

        assert np.allclose(wdata.t, ts)
        assert np.allclose(wdata.t0, 0)
        assert np.allclose(wdata.dt, 1)

        (density_, current_) = wdata.variables

        assert density_.name == "density"
        assert density_.description == ""
        assert density_.ext == "wdat"
        assert density_.unit == "none"
        assert density_.filename == f"{full_prefix}_density.wdat"
        assert density_.descr == "<f8"
        assert not density_.vector
        assert np.allclose(density_.data, density)

        assert current_.name == "current"
        assert current_.description == ""
        assert current_.ext == "wdat"
        assert current_.unit == "none"
        assert current_.filename == f"{full_prefix}_current.wdat"
        assert current_.descr == "<f8"
        assert current_.vector
        assert np.allclose(current_.data, current)
