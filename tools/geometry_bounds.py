#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path


INF = float("inf")
ATOL = 1.0e-12


@dataclass(frozen=True)
class BBox:
    lower: tuple[float, float, float]
    upper: tuple[float, float, float]

    @classmethod
    def infinite(cls) -> BBox:
        return cls((-INF, -INF, -INF), (INF, INF, INF))

    @classmethod
    def empty(cls) -> BBox:
        return cls((INF, INF, INF), (-INF, -INF, -INF))

    def intersect(self, other: BBox) -> BBox:
        return BBox(
            tuple(max(a, b) for a, b in zip(self.lower, other.lower)),
            tuple(min(a, b) for a, b in zip(self.upper, other.upper)),
        )

    def union(self, other: BBox) -> BBox:
        return BBox(
            tuple(min(a, b) for a, b in zip(self.lower, other.lower)),
            tuple(max(a, b) for a, b in zip(self.upper, other.upper)),
        )

    def expand(self, padding: float) -> BBox:
        return BBox(
            tuple(x - padding for x in self.lower),
            tuple(x + padding for x in self.upper),
        )

    def is_finite(self) -> bool:
        return all(math.isfinite(x) for x in (*self.lower, *self.upper))


@dataclass(frozen=True)
class Surface:
    surface_id: int
    surface_type: str
    coeffs: tuple[float, ...]

    def halfspace_bbox(self, side: str) -> BBox:
        if self.surface_type in {"x-plane", "y-plane", "z-plane"}:
            axis = {"x-plane": 0, "y-plane": 1, "z-plane": 2}[self.surface_type]
            return _axis_plane_bbox(axis, 1.0, self.coeffs[0], side)

        if self.surface_type == "plane":
            return _general_plane_bbox(self.coeffs, side)

        if side == "+":
            return BBox.infinite()

        if self.surface_type == "x-cylinder":
            y0, z0, radius = self.coeffs
            return BBox((-INF, y0 - radius, z0 - radius), (INF, y0 + radius, z0 + radius))

        if self.surface_type == "y-cylinder":
            x0, z0, radius = self.coeffs
            return BBox((x0 - radius, -INF, z0 - radius), (x0 + radius, INF, z0 + radius))

        if self.surface_type == "z-cylinder":
            x0, y0, radius = self.coeffs
            return BBox((x0 - radius, y0 - radius, -INF), (x0 + radius, y0 + radius, INF))

        if self.surface_type == "sphere":
            x0, y0, z0, radius = self.coeffs
            return BBox(
                (x0 - radius, y0 - radius, z0 - radius),
                (x0 + radius, y0 + radius, z0 + radius),
            )

        # OpenMC also returns an infinite halfspace bounding box for cones,
        # generic quadrics, and non-axis-aligned surfaces.
        return BBox.infinite()

    def contextual_halfspace_bbox(self, side: str, context: BBox) -> BBox | None:
        if side != "-":
            return None

        if self.surface_type == "x-cone" and _axis_is_finite(context, 0):
            x0, y0, z0, r2 = self.coeffs
            radius = math.sqrt(r2) * max(abs(context.lower[0] - x0), abs(context.upper[0] - x0))
            return BBox((-INF, y0 - radius, z0 - radius), (INF, y0 + radius, z0 + radius))

        if self.surface_type == "y-cone" and _axis_is_finite(context, 1):
            x0, y0, z0, r2 = self.coeffs
            radius = math.sqrt(r2) * max(abs(context.lower[1] - y0), abs(context.upper[1] - y0))
            return BBox((x0 - radius, -INF, z0 - radius), (x0 + radius, INF, z0 + radius))

        if self.surface_type == "z-cone" and _axis_is_finite(context, 2):
            x0, y0, z0, r2 = self.coeffs
            radius = math.sqrt(r2) * max(abs(context.lower[2] - z0), abs(context.upper[2] - z0))
            return BBox((x0 - radius, y0 - radius, -INF), (x0 + radius, y0 + radius, INF))

        return None


def _axis_plane_bbox(axis: int, coefficient: float, value: float, side: str) -> BBox:
    lower = [-INF, -INF, -INF]
    upper = [INF, INF, INF]
    if side == "-":
        if coefficient > 0.0:
            upper[axis] = value
        else:
            lower[axis] = value
    else:
        if coefficient > 0.0:
            lower[axis] = value
        else:
            upper[axis] = value
    return BBox(tuple(lower), tuple(upper))


def _general_plane_bbox(coeffs: tuple[float, ...], side: str) -> BBox:
    if len(coeffs) != 4:
        return BBox.infinite()

    a, b, c, d = coeffs
    nonzero = [(axis, value) for axis, value in enumerate((a, b, c)) if abs(value) > ATOL]
    if len(nonzero) != 1:
        return BBox.infinite()

    axis, coefficient = nonzero[0]
    return _axis_plane_bbox(axis, coefficient, d / coefficient, side)


def _axis_is_finite(box: BBox, axis: int) -> bool:
    return math.isfinite(box.lower[axis]) and math.isfinite(box.upper[axis])


class Region:
    def bbox(self, surfaces: dict[int, Surface]) -> BBox:
        raise NotImplementedError

    def complement(self) -> Region:
        raise NotImplementedError


@dataclass(frozen=True)
class Halfspace(Region):
    surface_id: int
    side: str

    def bbox(self, surfaces: dict[int, Surface]) -> BBox:
        try:
            surface = surfaces[self.surface_id]
        except KeyError as exc:
            raise ValueError(f"region references missing surface id {self.surface_id}") from exc
        return surface.halfspace_bbox(self.side)

    def complement(self) -> Region:
        return Halfspace(self.surface_id, "+" if self.side == "-" else "-")


@dataclass(frozen=True)
class Intersection(Region):
    nodes: tuple[Region, ...]

    def bbox(self, surfaces: dict[int, Surface]) -> BBox:
        box = BBox.infinite()
        for node in self.nodes:
            box = box.intersect(node.bbox(surfaces))
        for node in self.nodes:
            if not isinstance(node, Halfspace):
                continue
            surface = surfaces[node.surface_id]
            contextual_box = surface.contextual_halfspace_bbox(node.side, box)
            if contextual_box is not None:
                box = box.intersect(contextual_box)
        return box

    def complement(self) -> Region:
        return Union(tuple(node.complement() for node in self.nodes))


@dataclass(frozen=True)
class Union(Region):
    nodes: tuple[Region, ...]

    def bbox(self, surfaces: dict[int, Surface]) -> BBox:
        box = BBox.empty()
        for node in self.nodes:
            box = box.union(node.bbox(surfaces))
        return box

    def complement(self) -> Region:
        return Intersection(tuple(node.complement() for node in self.nodes))


@dataclass(frozen=True)
class Complement(Region):
    node: Region

    def bbox(self, surfaces: dict[int, Surface]) -> BBox:
        return self.node.complement().bbox(surfaces)

    def complement(self) -> Region:
        return self.node


def _is_integer_token(token: str) -> bool:
    return bool(re.fullmatch(r"[+-]?\d+", token))


def _ends_region(token: str) -> bool:
    return _is_integer_token(token) or token == ")"


def _starts_region(token: str) -> bool:
    return _is_integer_token(token) or token in {"(", "~"}


def _tokenize_region(expression: str) -> list[str]:
    raw_tokens = re.findall(r"[()]|\||~|[+-]?\d+", expression)
    compact_expr = re.sub(r"\s+", "", expression)
    if "".join(raw_tokens) != compact_expr:
        raise ValueError(f"unsupported character in region expression: {expression!r}")

    tokens: list[str] = []
    previous: str | None = None
    for token in raw_tokens:
        if previous is not None and _ends_region(previous) and _starts_region(token):
            tokens.append(" ")
        tokens.append(token)
        previous = token
    return tokens


def parse_region(expression: str) -> Region:
    output: list[Region] = []
    operators: list[str] = []
    precedence = {"|": 1, " ": 2, "~": 3}
    associativity = {"|": "left", " ": "left", "~": "right"}

    def make_halfspace(token: str) -> Halfspace:
        value = int(token)
        return Halfspace(abs(value), "-" if value < 0 else "+")

    def apply_operator(operator: str) -> None:
        if operator == "~":
            if not output:
                raise ValueError(f"missing operand for '~' in region expression: {expression!r}")
            output.append(Complement(output.pop()))
            return

        if len(output) < 2:
            raise ValueError(f"missing operand for '{operator}' in region expression: {expression!r}")
        right = output.pop()
        left = output.pop()
        if operator == " ":
            if isinstance(left, Intersection):
                output.append(Intersection((*left.nodes, right)))
            elif isinstance(right, Intersection):
                output.append(Intersection((left, *right.nodes)))
            else:
                output.append(Intersection((left, right)))
        elif operator == "|":
            if isinstance(left, Union):
                output.append(Union((*left.nodes, right)))
            elif isinstance(right, Union):
                output.append(Union((left, *right.nodes)))
            else:
                output.append(Union((left, right)))
        else:
            raise ValueError(f"unknown operator {operator!r}")

    for token in _tokenize_region(expression):
        if token in precedence:
            while operators:
                operator = operators[-1]
                if operator in {"(", ")"}:
                    break
                should_apply = (
                    associativity[token] == "right" and precedence[token] < precedence[operator]
                ) or (
                    associativity[token] == "left" and precedence[token] <= precedence[operator]
                )
                if not should_apply:
                    break
                apply_operator(operators.pop())
            operators.append(token)
        elif token == "(":
            operators.append(token)
        elif token == ")":
            while operators and operators[-1] != "(":
                apply_operator(operators.pop())
            if not operators:
                raise ValueError(f"mismatched parentheses in region expression: {expression!r}")
            operators.pop()
        elif _is_integer_token(token):
            output.append(make_halfspace(token))
        else:
            raise ValueError(f"bad token {token!r} in region expression: {expression!r}")

    while operators:
        operator = operators.pop()
        if operator in {"(", ")"}:
            raise ValueError(f"mismatched parentheses in region expression: {expression!r}")
        apply_operator(operator)

    if len(output) != 1:
        raise ValueError(f"could not parse region expression: {expression!r}")
    return output[0]


def get_xml_text(elem: ET.Element, name: str, default: str | None = None) -> str | None:
    value = elem.get(name)
    if value is not None:
        return value
    child = elem.find(name)
    if child is not None and child.text is not None:
        return child.text
    return default


def parse_surfaces(root: ET.Element) -> dict[int, Surface]:
    surfaces: dict[int, Surface] = {}
    for elem in root.findall("surface"):
        surface_id = int(get_xml_text(elem, "id", "0"))
        surface_type = get_xml_text(elem, "type")
        coeff_text = get_xml_text(elem, "coeffs", "")
        if surface_type is None:
            raise ValueError(f"surface {surface_id} is missing type")
        coeffs = tuple(float(value) for value in coeff_text.split())
        surfaces[surface_id] = Surface(surface_id, surface_type, coeffs)
    return surfaces


@dataclass(frozen=True)
class Cell:
    universe_id: int
    fill_id: int | None
    region: str | None


def parse_cells(root: ET.Element) -> list[Cell]:
    cells: list[Cell] = []
    for elem in root.findall("cell"):
        universe_id = int(get_xml_text(elem, "universe", "0"))
        has_material = get_xml_text(elem, "material") is not None
        fill_text = get_xml_text(elem, "fill")
        fill_id = None if has_material or fill_text is None else int(fill_text)
        cells.append(Cell(universe_id, fill_id, get_xml_text(elem, "region")))
    return cells


def infer_root_universe(cells: list[Cell]) -> int:
    universes = {cell.universe_id for cell in cells}
    child_universes = {cell.fill_id for cell in cells if cell.fill_id is not None}
    candidates = sorted(universes - child_universes)
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ValueError("could not infer root universe; pass --root-universe")
    joined = ", ".join(str(value) for value in candidates)
    raise ValueError(f"ambiguous root universe candidates ({joined}); pass --root-universe")


def geometry_bbox(path: Path, root_universe: int | None) -> tuple[int, BBox]:
    tree = ET.parse(path)
    root = tree.getroot()
    if root.tag != "geometry":
        raise ValueError(f"{path} does not look like an OpenMC geometry.xml")

    surfaces = parse_surfaces(root)
    cells = parse_cells(root)
    if root_universe is None:
        root_universe = infer_root_universe(cells)

    box = BBox.empty()
    matched_regions = 0
    for cell in cells:
        if cell.universe_id != root_universe or cell.region is None:
            continue
        matched_regions += 1
        box = box.union(parse_region(cell.region).bbox(surfaces))

    if matched_regions == 0:
        raise ValueError(f"root universe {root_universe} has no bounded cell regions")

    return root_universe, box


def format_float(value: float, digits: int) -> str:
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    text = f"{value:.{digits}g}"
    if "e" not in text.lower() and "." not in text:
        text += ".0"
    return text


def format_tuple(values: tuple[float, float, float], digits: int) -> str:
    return "(" + ", ".join(format_float(value, digits) for value in values) + ")"


def format_xml_vec(values: tuple[float, float, float], digits: int) -> str:
    return " ".join(format_float(value, digits) for value in values)


def print_python(box: BBox, pitch: float, digits: int) -> None:
    print("kinetics_mesh={")
    print(f"    'pitch': {format_float(pitch, digits)},")
    print("    'auto_bounds': False,")
    print(f"    'lower_left': {format_tuple(box.lower, digits)},")
    print(f"    'upper_right': {format_tuple(box.upper, digits)},")
    print("}")


def print_xml(box: BBox, pitch: float, digits: int) -> None:
    print("<kinetics_mesh>")
    print(f"  <pitch>{format_float(pitch, digits)}</pitch>")
    print("  <auto_bounds>false</auto_bounds>")
    print(f"  <lower_left>{format_xml_vec(box.lower, digits)}</lower_left>")
    print(f"  <upper_right>{format_xml_vec(box.upper, digits)}</upper_right>")
    print("</kinetics_mesh>")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read an OpenMC geometry.xml and print CLUTCH kinetics mesh bounds."
    )
    parser.add_argument(
        "geometry_xml",
        nargs="?",
        default="geometry.xml",
        type=Path,
        help="Path to geometry.xml. Defaults to ./geometry.xml.",
    )
    parser.add_argument(
        "--root-universe",
        type=int,
        help="Root universe id. If omitted, it is inferred from fill relationships.",
    )
    parser.add_argument(
        "--pitch",
        type=float,
        default=1.0,
        help="Pitch value to include in the printed kinetics_mesh snippet.",
    )
    parser.add_argument(
        "--padding",
        type=float,
        default=0.0,
        help="Optional padding added to every lower/upper bound.",
    )
    parser.add_argument(
        "--format",
        choices=("python", "xml", "both"),
        default="both",
        help="Output format. Defaults to both.",
    )
    parser.add_argument(
        "--digits",
        type=int,
        default=12,
        help="Significant digits used when printing floats.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        root_universe, box = geometry_bbox(args.geometry_xml, args.root_universe)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.padding:
        if args.padding < 0.0:
            print("error: --padding must be non-negative", file=sys.stderr)
            return 2
        box = box.expand(args.padding)

    print(f"# geometry_xml: {args.geometry_xml}")
    print(f"# root_universe: {root_universe}")
    if not box.is_finite():
        print(
            "# warning: calculated bounds contain infinity; provide manual bounds "
            "or check for missing finite boundary surfaces",
            file=sys.stderr,
        )

    if args.format in {"python", "both"}:
        print_python(box, args.pitch, args.digits)
    if args.format == "both":
        print()
    if args.format in {"xml", "both"}:
        print_xml(box, args.pitch, args.digits)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
