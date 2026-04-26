"""Python-only 2DAlphabet-style pass/fail validation fitter.

This module mirrors the core 2DAlphabet model algebra without ROOT/Combine:

    qcd_fail[x, y] = one free positive bin parameter per bin
    qcd_pass[x, y] = qcd_fail[x, y] * rpf(x, y; theta)

The transfer function ``rpf`` is a positive parametric function evaluated at
bin centers mapped to [0, 1], matching ``TwoDAlphabet.alphawrap.ParametricFunction``.
For speed, the free fail-bin parameters are analytically profiled for each
transfer-function parameter point.

Inputs are intentionally limited to ``hist.Hist`` objects or NumPy arrays.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from types import SimpleNamespace
from typing import Mapping, Sequence

import numpy as np


_EPS = 1.0e-9


@dataclass(frozen=True)
class Hist2D:
    """A small NumPy representation of a 2D histogram."""

    values: np.ndarray
    x_edges: np.ndarray
    y_edges: np.ndarray
    variances: np.ndarray | None = None
    x_name: str = "x"
    y_name: str = "y"

    @classmethod
    def from_hist(cls, histo) -> "Hist2D":
        """Build from a 2D scikit-hep ``hist.Hist`` object."""

        if len(histo.axes) != 2:
            raise ValueError(f"Expected a 2D hist.Hist, got {len(histo.axes)} axes")

        variances = histo.variances()
        return cls(
            values=np.asarray(histo.values(), dtype=float),
            variances=None if variances is None else np.asarray(variances, dtype=float),
            x_edges=np.asarray(histo.axes[0].edges, dtype=float),
            y_edges=np.asarray(histo.axes[1].edges, dtype=float),
            x_name=histo.axes[0].name,
            y_name=histo.axes[1].name,
        )

    @classmethod
    def from_numpy(
        cls,
        values: np.ndarray,
        x_edges: Sequence[float],
        y_edges: Sequence[float],
        variances: np.ndarray | None = None,
        x_name: str = "x",
        y_name: str = "y",
    ) -> "Hist2D":
        """Build from a raw NumPy array and bin edges."""

        return cls(
            values=np.asarray(values, dtype=float),
            variances=None if variances is None else np.asarray(variances, dtype=float),
            x_edges=np.asarray(x_edges, dtype=float),
            y_edges=np.asarray(y_edges, dtype=float),
            x_name=x_name,
            y_name=y_name,
        )

    def checked(self, name: str) -> "Hist2D":
        """Validate dimensional consistency and return ``self``."""

        if self.values.ndim != 2:
            raise ValueError(f"{name} must be 2D, got shape {self.values.shape}")
        expected = (len(self.x_edges) - 1, len(self.y_edges) - 1)
        if self.values.shape != expected:
            raise ValueError(f"{name} shape {self.values.shape} does not match edges {expected}")
        if self.variances is not None and self.variances.shape != self.values.shape:
            raise ValueError(f"{name} variances shape {self.variances.shape} does not match values")
        return self


@dataclass(frozen=True)
class ParameterSpec:
    """Configuration for one transfer-function parameter."""

    name: str
    nominal: float = 0.1
    minimum: float = -1000.0
    maximum: float = 1000.0
    error: float = 0.1
    constraint: str = "flatParam"

    @classmethod
    def from_2dalphabet(
        cls,
        name: str,
        spec: Mapping[str, float | str] | None = None,
    ) -> "ParameterSpec":
        """Build from a 2DAlphabet-like constraint dictionary."""

        spec = spec or {}
        return cls(
            name=name,
            nominal=float(spec.get("NOM", 0.1)),
            minimum=float(spec.get("MIN", -1000.0)),
            maximum=float(spec.get("MAX", 1000.0)),
            error=float(spec.get("ERROR", 0.1)),
            constraint=str(spec.get("constraint", "flatParam")),
        )


class FormulaTransferFunction:
    """2DAlphabet-style parametric transfer function.

    The formula syntax follows the useful subset of RooFit formulas used by
    2DAlphabet: parameters are ``@0``, ``@1``, ..., and mapped bin-center
    coordinates are ``x`` and ``y``.
    """

    _allowed_names = {
        "abs": np.abs,
        "acos": np.arccos,
        "asin": np.arcsin,
        "atan": np.arctan,
        "cos": np.cos,
        "cosh": np.cosh,
        "exp": np.exp,
        "log": np.log,
        "log10": np.log10,
        "max": np.maximum,
        "min": np.minimum,
        "pow": np.power,
        "sin": np.sin,
        "sinh": np.sinh,
        "sqrt": np.sqrt,
        "tan": np.tan,
        "tanh": np.tanh,
        "np": np,
    }

    def __init__(
        self,
        formula: str,
        constraints: Mapping[int, Mapping[str, float | str]] | None = None,
        force_positive: bool = True,
    ):
        self.formula = formula.replace(" ", "")
        self.force_positive = force_positive
        self.n_params = self._count_params(self.formula)
        constraints = constraints or {}
        self.parameters = [
            ParameterSpec.from_2dalphabet(f"rpf_par{i}", constraints.get(i))
            for i in range(self.n_params)
        ]
        self._python_expr = self._to_python_expr(self.formula)

    @staticmethod
    def _count_params(formula: str) -> int:
        matches = [int(match) for match in re.findall(r"@(\d+)", formula)]
        return max(matches) + 1 if matches else 0

    @staticmethod
    def _to_python_expr(formula: str) -> str:
        expr = formula.replace("^", "**")
        expr = re.sub(r"@(\d+)", r"p[\1]", expr)
        return expr

    @property
    def initial(self) -> np.ndarray:
        return np.asarray([p.nominal for p in self.parameters], dtype=float)

    @property
    def bounds(self) -> list[tuple[float, float]]:
        return [(p.minimum, p.maximum) for p in self.parameters]

    def evaluate(self, params: Sequence[float], x: np.ndarray, y: np.ndarray) -> np.ndarray:
        p = np.asarray(params, dtype=float)
        if len(p) != self.n_params:
            raise ValueError(f"Expected {self.n_params} parameters, got {len(p)}")
        value = eval(self._python_expr, {"__builtins__": {}}, {**self._allowed_names, "p": p, "x": x, "y": y})
        value = np.asarray(value, dtype=float)
        if value.shape == ():
            value = np.full_like(x, float(value), dtype=float)
        if self.force_positive:
            value = np.maximum(_EPS, value)
        return value

    def constraint_nll(self, params: Sequence[float]) -> float:
        out = 0.0
        for value, spec in zip(params, self.parameters):
            if spec.constraint == "flatParam":
                continue
            pieces = spec.constraint.split()
            if len(pieces) == 3 and pieces[0] == "param":
                mean = float(pieces[1])
                sigma = float(pieces[2])
                if sigma > 0:
                    out += 0.5 * ((value - mean) / sigma) ** 2
        return out


def polynomial_rpf(order_x: int, order_y: int, scale: float = 0.1) -> FormulaTransferFunction:
    """Create an additive 2D polynomial transfer function.

    Example: ``order_x=1, order_y=1`` gives
    ``scale * (@0 + @1*x + @2*y + @3*x*y)``.
    """

    terms = []
    ipar = 0
    for ix in range(order_x + 1):
        for iy in range(order_y + 1):
            factor = f"@{ipar}"
            if ix:
                factor += f"*x**{ix}" if ix > 1 else "*x"
            if iy:
                factor += f"*y**{iy}" if iy > 1 else "*y"
            terms.append(factor)
            ipar += 1
    return FormulaTransferFunction(f"{scale}*(" + "+".join(terms) + ")")


@dataclass(frozen=True)
class PassFailModelInput:
    """Inputs to the pass/fail 2DAlphabet-style likelihood."""

    data_fail: Hist2D
    data_pass: Hist2D
    bkg_fail: Hist2D | None = None
    bkg_pass: Hist2D | None = None
    signal_fail: Hist2D | None = None
    signal_pass: Hist2D | None = None

    def checked(self) -> "PassFailModelInput":
        self.data_fail.checked("data_fail")
        self.data_pass.checked("data_pass")
        _check_compatible("data_pass", self.data_fail, self.data_pass)
        for name in ("bkg_fail", "bkg_pass", "signal_fail", "signal_pass"):
            histo = getattr(self, name)
            if histo is not None:
                histo.checked(name)
                _check_compatible(name, self.data_fail, histo)
        return self


@dataclass(frozen=True)
class ABCDEFRegions:
    """Masks for the 2DAlphabet ABCDEF layout.

    The X axis is the top-candidate mass axis. LOW, SIG, and HIGH slices form
    the three alphabet columns:

        fail: A C E
        pass: B D F

    In the usual blinded background fit, pass-SIG (D) is excluded from the
    likelihood while fail-SIG (C) remains included to anchor the QCD template.
    """

    sig_start: float
    sig_end: float
    x_edges: np.ndarray
    low: np.ndarray
    sig: np.ndarray
    high: np.ndarray
    fail_fit_mask: np.ndarray
    pass_fit_mask: np.ndarray

    @classmethod
    def from_edges(
        cls,
        x_edges: Sequence[float],
        y_edges: Sequence[float],
        sig_start: float = 105.0,
        sig_end: float = 210.0,
        blind_pass_signal: bool = True,
    ) -> "ABCDEFRegions":
        """Create ABCDEF masks from 2D histogram bin edges."""

        x_edges = np.asarray(x_edges, dtype=float)
        y_edges = np.asarray(y_edges, dtype=float)
        if sig_start not in x_edges or sig_end not in x_edges:
            raise ValueError(
                "sig_start and sig_end must coincide with X-axis bin edges. "
                f"Got {sig_start}, {sig_end}; available edges are {x_edges.tolist()}"
            )

        x_low = x_edges[:-1]
        x_high = x_edges[1:]
        low_1d = x_high <= sig_start
        sig_1d = (x_low >= sig_start) & (x_high <= sig_end)
        high_1d = x_low >= sig_end
        if not np.all(low_1d | sig_1d | high_1d):
            raise ValueError("X-axis bins must not straddle ABCDEF region boundaries")

        shape = (len(x_edges) - 1, len(y_edges) - 1)
        low = np.broadcast_to(low_1d[:, None], shape).copy()
        sig = np.broadcast_to(sig_1d[:, None], shape).copy()
        high = np.broadcast_to(high_1d[:, None], shape).copy()
        fail_fit_mask = np.ones(shape, dtype=bool)
        pass_fit_mask = np.ones(shape, dtype=bool)
        if blind_pass_signal:
            pass_fit_mask &= ~sig

        return cls(
            sig_start=sig_start,
            sig_end=sig_end,
            x_edges=x_edges,
            low=low,
            sig=sig,
            high=high,
            fail_fit_mask=fail_fit_mask,
            pass_fit_mask=pass_fit_mask,
        )

    @property
    def A(self) -> np.ndarray:
        return self.low

    @property
    def B(self) -> np.ndarray:
        return self.low

    @property
    def C(self) -> np.ndarray:
        return self.sig

    @property
    def D(self) -> np.ndarray:
        return self.sig

    @property
    def E(self) -> np.ndarray:
        return self.high

    @property
    def F(self) -> np.ndarray:
        return self.high

    def yields(self, fail_values: np.ndarray, pass_values: np.ndarray) -> dict[str, float]:
        """Return ABCDEF yields for two full fail/pass templates."""

        return {
            "A": float(np.sum(fail_values[self.A])),
            "B": float(np.sum(pass_values[self.B])),
            "C": float(np.sum(fail_values[self.C])),
            "D": float(np.sum(pass_values[self.D])),
            "E": float(np.sum(fail_values[self.E])),
            "F": float(np.sum(pass_values[self.F])),
        }


@dataclass(frozen=True)
class FitResult:
    """Fit output with profiled QCD templates and expectations."""

    success: bool
    message: str
    nll: float
    params: dict[str, float]
    qcd_fail: np.ndarray
    qcd_pass: np.ndarray
    expected_fail: np.ndarray
    expected_pass: np.ndarray
    rpf: np.ndarray
    x_mapped: np.ndarray
    y_mapped: np.ndarray
    fail_fit_mask: np.ndarray
    pass_fit_mask: np.ndarray
    abcdef_yields: dict[str, dict[str, float]] | None
    scipy_result: object


def _check_compatible(name: str, reference: Hist2D, other: Hist2D) -> None:
    if other.values.shape != reference.values.shape:
        raise ValueError(f"{name} shape {other.values.shape} does not match {reference.values.shape}")
    if not np.allclose(other.x_edges, reference.x_edges):
        raise ValueError(f"{name} x edges do not match data_fail")
    if not np.allclose(other.y_edges, reference.y_edges):
        raise ValueError(f"{name} y edges do not match data_fail")


def _zero_like(reference: Hist2D) -> np.ndarray:
    return np.zeros_like(reference.values, dtype=float)


def _checked_mask(mask: np.ndarray | None, shape: tuple[int, int], name: str) -> np.ndarray:
    if mask is None:
        return np.ones(shape, dtype=bool)
    out = np.asarray(mask, dtype=bool)
    if out.shape != shape:
        raise ValueError(f"{name} shape {out.shape} does not match histogram shape {shape}")
    return out


def _mapped_centers(x_edges: np.ndarray, y_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    x = (x_centers - x_edges[0]) / (x_edges[-1] - x_edges[0])
    y = (y_centers - y_edges[0]) / (y_edges[-1] - y_edges[0])
    return np.meshgrid(x, y, indexing="ij")


def _poisson_nll(observed: np.ndarray, expected: np.ndarray, mask: np.ndarray | None = None) -> float:
    if mask is not None:
        observed = observed[mask]
        expected = expected[mask]
    expected = np.maximum(expected, _EPS)
    return float(np.sum(expected - observed * np.log(expected)))


def _profile_qcd_fail(
    data_fail: np.ndarray,
    data_pass: np.ndarray,
    bkg_fail: np.ndarray,
    bkg_pass: np.ndarray,
    rpf: np.ndarray,
    fail_mask: np.ndarray | None = None,
    pass_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Profile one positive fail-QCD bin parameter per bin.

    This solves the same free-bin likelihood maximization that Combine performs
    numerically for the 2DAlphabet ``BinnedDistribution`` fail object.
    """

    n_fail = np.maximum(data_fail, 0.0)
    n_pass = np.maximum(data_pass, 0.0)
    b_fail = np.maximum(bkg_fail, 0.0)
    b_pass = np.maximum(bkg_pass, 0.0)
    r = np.maximum(rpf, _EPS)
    fail_mask = np.ones_like(r, dtype=bool) if fail_mask is None else np.asarray(fail_mask, dtype=bool)
    pass_mask = np.ones_like(r, dtype=bool) if pass_mask is None else np.asarray(pass_mask, dtype=bool)

    if fail_mask.shape != r.shape or pass_mask.shape != r.shape:
        raise ValueError("fit masks must match histogram shape")

    q = np.full_like(r, _EPS, dtype=float)

    fail_only = fail_mask & ~pass_mask
    pass_only = pass_mask & ~fail_mask
    both = fail_mask & pass_mask
    neither = ~(fail_mask | pass_mask)

    q[fail_only] = np.maximum(n_fail[fail_only] - b_fail[fail_only], _EPS)
    q[pass_only] = np.maximum((n_pass[pass_only] - b_pass[pass_only]) / r[pass_only], _EPS)

    one_plus_r = 1.0 + r

    a = one_plus_r * r
    b = one_plus_r * (b_pass + r * b_fail) - r * (n_fail + n_pass)
    c = one_plus_r * b_fail * b_pass - n_fail * b_pass - n_pass * r * b_fail

    disc = np.maximum(b * b - 4.0 * a * c, 0.0)
    q[both] = (-b[both] + np.sqrt(disc[both])) / (2.0 * a[both])
    q[neither] = np.maximum(n_fail[neither] - b_fail[neither], _EPS)
    return np.maximum(q, _EPS)


class PassFail2DFitter:
    """Binned Poisson pass/fail fitter using the 2DAlphabet QCD model."""

    def __init__(
        self,
        inputs: PassFailModelInput,
        transfer_function: FormulaTransferFunction,
        fixed_signal_strength: float = 0.0,
        float_signal: bool = False,
        signal_bounds: tuple[float, float] = (-10.0, 10.0),
        signal_initial: float = 0.0,
        fail_fit_mask: np.ndarray | None = None,
        pass_fit_mask: np.ndarray | None = None,
        regions: ABCDEFRegions | None = None,
    ):
        self.inputs = inputs.checked()
        self.transfer_function = transfer_function
        self.fixed_signal_strength = fixed_signal_strength
        self.float_signal = float_signal
        self.signal_bounds = signal_bounds
        self.signal_initial = signal_initial
        shape = self.inputs.data_fail.values.shape
        self.fail_fit_mask = _checked_mask(fail_fit_mask, shape, "fail_fit_mask")
        self.pass_fit_mask = _checked_mask(pass_fit_mask, shape, "pass_fit_mask")
        self.regions = regions
        self.x_mapped, self.y_mapped = _mapped_centers(
            self.inputs.data_fail.x_edges, self.inputs.data_fail.y_edges
        )

    @classmethod
    def abcdef(
        cls,
        inputs: PassFailModelInput,
        transfer_function: FormulaTransferFunction,
        sig_start: float = 105.0,
        sig_end: float = 210.0,
        blind_pass_signal: bool = True,
        **kwargs,
    ) -> "PassFail2DFitter":
        """Construct the standard ABCDEF sideband fitter."""

        checked_inputs = inputs.checked()
        regions = ABCDEFRegions.from_edges(
            checked_inputs.data_fail.x_edges,
            checked_inputs.data_fail.y_edges,
            sig_start=sig_start,
            sig_end=sig_end,
            blind_pass_signal=blind_pass_signal,
        )
        return cls(
            checked_inputs,
            transfer_function,
            fail_fit_mask=regions.fail_fit_mask,
            pass_fit_mask=regions.pass_fit_mask,
            regions=regions,
            **kwargs,
        )

    @property
    def initial(self) -> np.ndarray:
        if not self.float_signal:
            return self.transfer_function.initial
        return np.r_[self.transfer_function.initial, self.signal_initial]

    @property
    def bounds(self) -> list[tuple[float, float]]:
        if not self.float_signal:
            return self.transfer_function.bounds
        return self.transfer_function.bounds + [self.signal_bounds]

    def _split_params(self, params: Sequence[float]) -> tuple[np.ndarray, float]:
        params = np.asarray(params, dtype=float)
        if not self.float_signal:
            return params, self.fixed_signal_strength
        return params[:-1], float(params[-1])

    def templates_for(self, params: Sequence[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        tf_params, signal_strength = self._split_params(params)
        rpf = self.transfer_function.evaluate(tf_params, self.x_mapped, self.y_mapped)

        bkg_fail = _zero_like(self.inputs.data_fail) if self.inputs.bkg_fail is None else self.inputs.bkg_fail.values
        bkg_pass = _zero_like(self.inputs.data_pass) if self.inputs.bkg_pass is None else self.inputs.bkg_pass.values
        if self.inputs.signal_fail is not None:
            bkg_fail = bkg_fail + signal_strength * self.inputs.signal_fail.values
        if self.inputs.signal_pass is not None:
            bkg_pass = bkg_pass + signal_strength * self.inputs.signal_pass.values

        qcd_fail = _profile_qcd_fail(
            self.inputs.data_fail.values,
            self.inputs.data_pass.values,
            bkg_fail,
            bkg_pass,
            rpf,
            self.fail_fit_mask,
            self.pass_fit_mask,
        )
        qcd_pass = qcd_fail * rpf
        expected_fail = qcd_fail + bkg_fail
        expected_pass = qcd_pass + bkg_pass
        return qcd_fail, qcd_pass, expected_fail, expected_pass, rpf

    def nll(self, params: Sequence[float]) -> float:
        tf_params, _ = self._split_params(params)
        _, _, expected_fail, expected_pass, _ = self.templates_for(params)
        if not np.all(np.isfinite(expected_fail)) or not np.all(np.isfinite(expected_pass)):
            return math.inf

        out = _poisson_nll(self.inputs.data_fail.values, expected_fail, self.fail_fit_mask)
        out += _poisson_nll(self.inputs.data_pass.values, expected_pass, self.pass_fit_mask)
        out += self.transfer_function.constraint_nll(tf_params)
        return out

    def fit(self, **minimize_kwargs) -> FitResult:
        """Fit and return profiled templates.

        Uses SciPy's L-BFGS-B when available. If SciPy is not installed, falls
        back to a small bounded pattern search that is adequate for low-order
        transfer-function validation scans.
        """

        options = {"maxiter": 2000, "ftol": 1.0e-8}
        options.update(minimize_kwargs.pop("options", {}))
        try:
            from scipy.optimize import minimize

            result = minimize(
                self.nll,
                self.initial,
                method=minimize_kwargs.pop("method", "L-BFGS-B"),
                bounds=self.bounds,
                options=options,
                **minimize_kwargs,
            )
        except ModuleNotFoundError:
            if minimize_kwargs:
                unknown = ", ".join(minimize_kwargs)
                raise TypeError(f"Pure NumPy fallback does not accept minimize kwargs: {unknown}")
            result = self._fit_pattern_search(options)

        qcd_fail, qcd_pass, expected_fail, expected_pass, rpf = self.templates_for(result.x)
        tf_params, signal_strength = self._split_params(result.x)

        params = {
            spec.name: float(value)
            for spec, value in zip(self.transfer_function.parameters, tf_params)
        }
        if self.float_signal:
            params["signal_strength"] = signal_strength
        abcdef_yields = None
        if self.regions is not None:
            abcdef_yields = {
                "data": self.regions.yields(self.inputs.data_fail.values, self.inputs.data_pass.values),
                "qcd": self.regions.yields(qcd_fail, qcd_pass),
                "expected": self.regions.yields(expected_fail, expected_pass),
            }

        return FitResult(
            success=bool(result.success),
            message=str(result.message),
            nll=float(result.fun),
            params=params,
            qcd_fail=qcd_fail,
            qcd_pass=qcd_pass,
            expected_fail=expected_fail,
            expected_pass=expected_pass,
            rpf=rpf,
            x_mapped=self.x_mapped,
            y_mapped=self.y_mapped,
            fail_fit_mask=self.fail_fit_mask,
            pass_fit_mask=self.pass_fit_mask,
            abcdef_yields=abcdef_yields,
            scipy_result=result,
        )

    def _fit_pattern_search(self, options: Mapping[str, float | int]) -> SimpleNamespace:
        maxiter = int(options.get("maxiter", 2000))
        ftol = float(options.get("ftol", 1.0e-8))
        x = np.asarray(self.initial, dtype=float)
        bounds = self.bounds

        widths = np.asarray(
            [
                high - low if np.isfinite(high - low) else max(abs(start), 1.0)
                for start, (low, high) in zip(x, bounds)
            ],
            dtype=float,
        )
        param_errors = np.asarray(
            [spec.error for spec in self.transfer_function.parameters],
            dtype=float,
        )
        if self.float_signal:
            param_errors = np.r_[param_errors, 0.2]
        steps = np.maximum(param_errors, 0.05 * widths)
        steps = np.where(np.isfinite(steps), steps, 1.0)

        def clip_to_bounds(values):
            out = np.array(values, copy=True, dtype=float)
            for i, (low, high) in enumerate(bounds):
                out[i] = min(max(out[i], low), high)
            return out

        x = clip_to_bounds(x)
        best = self.nll(x)
        nfev = 1

        for iteration in range(maxiter):
            improved = False
            start_best = best
            for ipar in range(len(x)):
                for direction in (1.0, -1.0):
                    trial = np.array(x, copy=True)
                    trial[ipar] += direction * steps[ipar]
                    trial = clip_to_bounds(trial)
                    value = self.nll(trial)
                    nfev += 1
                    if value + ftol < best:
                        x = trial
                        best = value
                        improved = True

            if not improved:
                steps *= 0.5
            if abs(start_best - best) < ftol and np.max(steps) < 1.0e-4:
                return SimpleNamespace(
                    x=x,
                    fun=best,
                    success=True,
                    message="Converged with pure NumPy bounded pattern search",
                    nit=iteration + 1,
                    nfev=nfev,
                )

        return SimpleNamespace(
            x=x,
            fun=best,
            success=False,
            message="Reached maxiter in pure NumPy bounded pattern search",
            nit=maxiter,
            nfev=nfev,
        )


def subtract_background(data: Hist2D, backgrounds: Sequence[Hist2D]) -> Hist2D:
    """Return the 2DAlphabet ``InitQCDHists`` equivalent: data minus backgrounds."""

    values = np.array(data.values, copy=True, dtype=float)
    variances = None if data.variances is None else np.array(data.variances, copy=True, dtype=float)
    for i, background in enumerate(backgrounds):
        _check_compatible(f"background[{i}]", data, background)
        values -= background.values
        if variances is not None and background.variances is not None:
            variances += background.variances
    return Hist2D(
        values=values,
        variances=variances,
        x_edges=data.x_edges,
        y_edges=data.y_edges,
        x_name=data.x_name,
        y_name=data.y_name,
    )
