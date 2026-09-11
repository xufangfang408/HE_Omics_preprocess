#!/usr/bin/env python3
"""Batch H&E masks with TIAToolbox; Python >=3.11.

Install in a clean environment: pip install tiatoolbox==2.1.3
Batch: python he_tissue_mask.py --HE-dir HE_dir --mask-dir mask_dir
Coarse outline (interior holes removed): add --mask-mode coarse
One image: add --only section01.tif --set sensitivity=1.2 --overwrite
Export settings: python he_tissue_mask.py --init-config config.json
Persistent overrides: edit config.json, then add --config config.json.

Outputs under --mask-dir, named after the input stem (section01.tif ->):
    mask/section01.tif        single-channel binary TIFF (0/255)
    preview/section01.jpg     H&E | mask | overlay panel
    json/section01.json       parameters, sizes, coordinate scale, diagnostics
    json/section01.error.json only when that image fails
Stems must be unique per subdirectory; section01.tif and section01.png collide.
Default mask is at working resolution, NOT necessarily original resolution.
--full-size exports nearest-neighbour original-size mask (bounded by max pixels).
No GPU/weights required. Enhanced mode is a custom heuristic built around
TIAToolbox Otsu, not an official TIAToolbox model or validated universal method.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi

DEFAULTS = dict(
    mode="enhanced", mask_mode="fine", max_dim=3000, work_mpp=None, input_mpp=None,
    sensitivity=1.0, gray_threshold=None, background_rgb=None,
    background_correct=True, od_floor=0.025, saturation_floor=0.025,
    blur_sigma=0.6, min_area=150, hole_area=100, close_radius=1,
    dilate_radius=0, keep_largest=0,
    coarse_close_radius=6, coarse_open_radius=2, coarse_min_area=2000,
    coarse_fill_holes=True, coarse_max_hole_frac=0.5, edge_artifact_ratio=0.8,
)
EXTENSIONS = {".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff",
              ".svs", ".ndpi", ".mrxs", ".scn", ".vms", ".vmu"}
RASTER = {".png", ".jpg", ".jpeg", ".bmp"}
INTEGERS = ("max_dim", "min_area", "hole_area", "close_radius", "dilate_radius",
            "keep_largest", "coarse_close_radius", "coarse_open_radius", "coarse_min_area")


def validate(p):
    if set(p) != set(DEFAULTS):
        raise ValueError(f"Unknown parameters: {set(p) - set(DEFAULTS)}")
    if p["mode"] not in {"enhanced", "otsu"}:
        raise ValueError("mode must be enhanced or otsu")
    if p["mask_mode"] not in {"fine", "coarse"}:
        raise ValueError("mask_mode must be fine or coarse")
    if not 0.25 <= p["sensitivity"] <= 3:
        raise ValueError("sensitivity must be 0.25..3")
    for k in ("background_correct", "coarse_fill_holes"):
        if not isinstance(p[k], bool):
            raise ValueError(f"{k} must be true or false")
    for k in INTEGERS:
        if isinstance(p[k], bool) or not isinstance(p[k], int) or p[k] < (16 if k == "max_dim" else 0):
            raise ValueError(f"Invalid nonnegative integer: {k}")
    for k in ("work_mpp", "input_mpp"):
        if p[k] is not None and (not np.isfinite(p[k]) or p[k] <= 0):
            raise ValueError(f"{k} must be positive")
    for k in ("blur_sigma", "od_floor", "saturation_floor"):
        if not np.isfinite(p[k]) or p[k] < 0:
            raise ValueError(f"Invalid {k}")
    for k in ("coarse_max_hole_frac", "edge_artifact_ratio"):
        if not np.isfinite(p[k]) or not 0 <= p[k] <= 1:
            raise ValueError(f"{k} must be 0..1")
    if p["gray_threshold"] is not None and not 0 <= p["gray_threshold"] <= 255:
        raise ValueError("gray_threshold must be 0..255")
    if p["background_rgb"] is not None:
        bg = np.asarray(p["background_rgb"], dtype=float)
        if bg.shape != (3,) or not np.isfinite(bg).all() or np.any((bg <= 0) | (bg > 255)):
            raise ValueError("background_rgb must be three values in (0,255]")


def rgb8(a):
    a = np.asarray(a)
    if a.dtype != np.uint8 or a.ndim != 3 or a.shape[-1] not in (3, 4):
        raise ValueError("Expected 8-bit RGB/RGBA H&E; explicitly convert other bit depths first")
    if a.shape[-1] == 4:
        alpha = a[..., 3:4].astype(np.float32) / 255
        a = np.rint(a[..., :3] * alpha + 255 * (1 - alpha)).astype(np.uint8)
    return a


def repair_edge(rgb, ratio, notes):
    """Thumbnail readers blend a dark 1-px frame where the filter runs off the slide.

    That frame is darker than any real background, so Otsu keeps it and the mask gains
    a closed ring along the image border; coarse mode would then fill the whole image.
    """
    if not ratio or min(rgb.shape[:2]) < 4:
        return rgb
    edges = {"top": (rgb[0], rgb[1]), "bottom": (rgb[-1], rgb[-2]),
             "left": (rgb[:, 0], rgb[:, 1]), "right": (rgb[:, -1], rgb[:, -2])}
    dark = [k for k, (edge, inner) in edges.items() if edge.mean() < ratio * inner.mean()]
    if not dark:
        return rgb
    rgb = rgb.copy()
    if "top" in dark:
        rgb[0] = rgb[1]
    if "bottom" in dark:
        rgb[-1] = rgb[-2]
    if "left" in dark:
        rgb[:, 0] = rgb[:, 1]
    if "right" in dark:
        rgb[:, -1] = rgb[:, -2]
    notes.append("Replaced dark 1-px reader border (" + ", ".join(sorted(dark)) + ")")
    return rgb


def read_image(path, p):
    notes = []
    if path.suffix.lower() in RASTER:
        with Image.open(path) as im:
            if im.mode not in ("RGB", "RGBA"):
                raise ValueError(f"Expected RGB/RGBA, got {im.mode}")
            wh = im.size
            mpp = None if p["input_mpp"] is None else [p["input_mpp"]] * 2
            scale = min(1., p["max_dim"] / max(wh))
            if p["work_mpp"] is not None:
                if mpp is None:
                    raise ValueError("Raster work_mpp requires input_mpp")
                scale = min(scale, mpp[0] / p["work_mpp"])
            size = tuple(max(1, round(v * scale)) for v in wh)
            # Composite before resizing, so transparent black does not bleed in.
            if im.mode == "RGBA":
                bg = Image.new("RGBA", im.size, "white")
                im = Image.alpha_composite(bg, im).convert("RGB")
            im = im.resize(size, Image.Resampling.LANCZOS)
            rgb = np.array(im)
    else:
        from tiatoolbox.wsicore.wsireader import WSIReader
        kwargs = {} if p["input_mpp"] is None else {"mpp": p["input_mpp"]}
        reader = WSIReader.open(path, **kwargs)
        try:
            wh = tuple(int(v) for v in reader.info.slide_dimensions)
            mpp = reader.info.mpp
            mpp = None if mpp is None else np.asarray(mpp).reshape(-1).tolist()
            if mpp is not None and len(mpp) == 1:
                mpp *= 2
            scale = min(1., p["max_dim"] / max(wh))
            if p["work_mpp"] is not None:
                if mpp is None:
                    raise ValueError("Missing MPP: set input_mpp, or leave work_mpp null")
                scale = min(scale, min(mpp) / p["work_mpp"])
            rgb = rgb8(reader.slide_thumbnail(resolution=scale, units="baseline"))
        finally:
            close = getattr(reader, "close", None)
            if callable(close):
                close()
    h, w = rgb.shape[:2]
    actual_mpp = None if mpp is None else [mpp[0] * wh[0] / w, mpp[1] * wh[1] / h]
    if p["work_mpp"] is not None and max(actual_mpp) > p["work_mpp"] * 1.05:
        notes.append("max_dim or native resolution limits requested work_mpp")
    return repair_edge(rgb8(rgb), p["edge_artifact_ratio"], notes), wh, actual_mpp, notes


def disk(radius):
    y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
    return x*x + y*y <= radius*radius


def remove_small(mask, area):
    labels, _ = ndi.label(mask, structure=np.ones((3, 3)))
    sizes = np.bincount(labels.ravel())
    keep = sizes >= area
    keep[0] = False
    return keep[labels]


def morph(mask, radius, closing):
    """Closing/opening padded with edge values, so tissue at the border survives."""
    if not radius:
        return mask
    op = ndi.binary_closing if closing else ndi.binary_opening
    padded = np.pad(mask, radius, mode="edge")
    return op(padded, structure=disk(radius))[radius:-radius, radius:-radius]


def cleanup(mask, p):
    mask = morph(mask, p["close_radius"], closing=True)
    mask = remove_small(mask, p["min_area"])
    if p["hole_area"]:
        labels, _ = ndi.label(~mask, structure=np.ones((3, 3)))
        sizes = np.bincount(labels.ravel())
        fill = sizes <= p["hole_area"]
        # Never fill background connected to the image boundary.
        border = np.unique(np.concatenate([labels[0], labels[-1], labels[:, 0], labels[:, -1]]))
        fill[border] = False
        fill[0] = False
        mask |= fill[labels]
    if p["keep_largest"]:
        labels, n = ndi.label(mask, structure=np.ones((3, 3)))
        sizes = np.bincount(labels.ravel())[1:]
        ids = np.argsort(sizes)[-min(n, p["keep_largest"]):] + 1 if n else []
        mask = np.isin(labels, ids)
    if p["dilate_radius"]:
        mask = ndi.binary_dilation(mask, structure=disk(p["dilate_radius"]))
    return mask


def fill_holes(mask, max_frac, notes):
    """Fill enclosed background, but refuse holes so large they are a ring artifact.

    Any foreground ring around the field of view encloses the whole background, and
    filling it would report the entire image as tissue instead of failing loudly.
    """
    filled = ndi.binary_fill_holes(mask)
    if max_frac >= 1:
        return filled
    labels, n = ndi.label(filled & ~mask, structure=np.ones((3, 3)))
    if not n:
        return filled
    oversized = np.bincount(labels.ravel()) > max_frac * mask.size
    oversized[0] = False
    if oversized.any():
        note = ("Skipped filling a hole larger than coarse_max_hole_frac: "
                "tissue may be ringed by an artifact")
        if note not in notes:
            notes.append(note)
        filled &= ~oversized[labels]
    return filled


def coarsen(mask, p, notes):
    """Outline-only mask: bridge intra-tissue gaps and drop every interior hole."""
    mask = morph(mask, p["coarse_close_radius"], closing=True)
    if p["coarse_fill_holes"]:
        mask = fill_holes(mask, p["coarse_max_hole_frac"], notes)
    # Opening runs after filling, so it smooths the outline instead of thin webbing.
    mask = morph(mask, p["coarse_open_radius"], closing=False)
    if p["coarse_fill_holes"]:
        mask = fill_holes(mask, p["coarse_max_hole_frac"], notes)
    if p["coarse_min_area"]:
        mask = remove_small(mask, p["coarse_min_area"])
    return mask


def extract_mask(rgb, p):
    from tiatoolbox.tools.tissuemask import OtsuTissueMasker
    from skimage.filters import threshold_otsu
    notes = []
    arr = rgb.astype(np.float32)
    # Bright, low-chroma candidates over the whole image, not just its border.
    sample = arr[::max(1, arr.shape[0] // 512), ::max(1, arr.shape[1] // 512)].reshape(-1, 3)
    bright = sample.mean(1)
    chroma = np.ptp(sample, axis=1)
    candidates = sample[(bright >= np.percentile(bright, 85)) & (chroma < 25) & (bright > 180)]
    if p["background_rgb"] is not None:
        bg = np.array(p["background_rgb"], dtype=np.float32)
    elif len(candidates) >= max(20, int(len(sample) * 0.005)):
        bg = np.median(candidates, axis=0)
    else:
        bg = np.array([255., 255., 255.])
        notes.append("No reliable bright background: inspect or set background_rgb")
    if p["background_correct"]:
        arr = np.clip(arr / bg * 255, 0, 255)
    if p["blur_sigma"]:
        arr = ndi.gaussian_filter(arr, sigma=(p["blur_sigma"], p["blur_sigma"], 0))
    corrected = np.rint(arr).astype(np.uint8)
    masker = OtsuTissueMasker()
    masker.fit(corrected[None])  # Fit independently for every section.
    original_threshold = float(masker.threshold)
    threshold = p["gray_threshold"]
    if threshold is None:
        threshold = 255 - (255 - original_threshold) / p["sensitivity"]
    masker.threshold = float(threshold)
    base = masker.transform(corrected[None])[0].astype(bool)
    intensity = arr / 255
    od = -np.log(np.clip(intensity, 1/255, 1)).mean(-1)
    sat = (intensity.max(-1) - intensity.min(-1)) / np.maximum(intensity.max(-1), 1/255)
    od_t = max(p["od_floor"], float(threshold_otsu(od))) / p["sensitivity"]
    sat_t = max(p["saturation_floor"], float(threshold_otsu(sat))) / p["sensitivity"]
    if p["mode"] == "enhanced":
        # Floors prevent nearly uniform blank slides being divided by Otsu noise.
        signal = (od >= p["od_floor"] / p["sensitivity"]) | (sat >= p["saturation_floor"] / p["sensitivity"])
        mask = (base | (od > od_t) | (sat > sat_t)) & signal
    else:
        mask = base
    mask = cleanup(mask, p)
    fine_coverage = float(mask.mean())
    if p["mask_mode"] == "coarse":
        mask = coarsen(mask, p, notes)
    coverage = float(mask.mean())
    if coverage < 0.01 or coverage > 0.95:
        notes.append("Extreme tissue fraction: review preview (may be legitimate)")
    return mask, dict(background_rgb=bg.tolist(), otsu_threshold=original_threshold,
                      gray_threshold=float(threshold), od_threshold=od_t,
                      saturation_threshold=sat_t, tissue_fraction=coverage,
                      fine_tissue_fraction=fine_coverage,
                      warnings=notes)


def save_preview(rgb, mask, path):
    im = Image.fromarray(rgb)
    im.thumbnail((1200, 1200))
    m = np.asarray(Image.fromarray(mask).resize(im.size, Image.Resampling.NEAREST))
    overlay = np.array(im).copy()
    overlay[m] = (0.65 * overlay[m] + 0.35 * np.array([0, 230, 90])).astype(np.uint8)
    panels = [im, Image.fromarray(m.astype(np.uint8)*255).convert("RGB"), Image.fromarray(overlay)]
    canvas = Image.new("RGB", (im.width*3, im.height+28), "white")
    draw = ImageDraw.Draw(canvas)
    for i, (panel, title) in enumerate(zip(panels, ["H&E", "Mask", "Overlay"])):
        canvas.paste(panel, (i*im.width, 28))
        draw.text((i*im.width+5, 6), title, fill="black")
    canvas.save(path, format="JPEG", quality=90)


def output_paths(out, rel):
    """json/<stem>.json, mask/<stem>.tif, preview/<stem>.jpg, mirroring subdirectories."""
    return {name: out / name / rel.parent / (rel.stem + suffix)
            for name, suffix in (("json", ".json"), ("mask", ".tif"), ("preview", ".jpg"))}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--HE-dir", type=Path)
    parser.add_argument("--mask-dir", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--init-config", type=Path)
    parser.add_argument("--only", help="Exact input-relative filename, e.g. batch1/section01.tif")
    parser.add_argument("--mask-mode", choices=["fine", "coarse"],
                        help="fine keeps interior holes; coarse keeps tissue outlines only")
    parser.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--full-size", action="store_true")
    parser.add_argument("--max-full-pixels", type=int, default=100_000_000)
    args = parser.parse_args()
    if args.init_config:
        if args.init_config.exists():
            parser.error(f"{args.init_config} already exists; delete or rename it first")
        with args.init_config.open("x", encoding="utf-8") as f:
            json.dump({"defaults": DEFAULTS, "overrides": {"section01.tif": {"sensitivity": 1.2}}}, f, indent=2)
        return 0
    if not args.HE_dir or not args.mask_dir:
        parser.error("--HE-dir and --mask-dir are required")
    root, out = args.HE_dir.resolve(), args.mask_dir.resolve()
    if not root.is_dir():
        parser.error("HE-dir is not a directory")
    if root == out or root in out.parents:
        parser.error("mask-dir must be outside HE-dir to avoid reading generated outputs")
    config = json.loads(args.config.read_text(encoding="utf-8")) if args.config else {}
    if set(config) - {"defaults", "overrides"}:
        parser.error("Config keys must be defaults and overrides")
    cli = {}
    if args.mask_mode:
        cli["mask_mode"] = args.mask_mode
    for entry in args.set:
        key, value = entry.split("=", 1)
        try:
            cli[key] = json.loads(value)
        except json.JSONDecodeError:
            cli[key] = value
    files = sorted(p for p in (root.rglob("*") if args.recursive else root.iterdir())
                   if p.is_file() and p.suffix.lower() in EXTENSIONS)
    # Outputs are named by stem, so two inputs sharing one stem would overwrite each other.
    stems = {}
    for path in files:
        rel = path.relative_to(root)
        stems.setdefault((rel.parent, rel.stem), []).append(rel.as_posix())
    clashes = sorted(v for v in stems.values() if len(v) > 1)
    if clashes:
        parser.error("Inputs share an output name; rename them first: "
                     + "; ".join(" vs ".join(c) for c in clashes))
    if args.only:
        files = [p for p in files if p.relative_to(root).as_posix() == args.only]
    if not files:
        parser.error("No matching images; check --only path and --recursive")
    failures = 0
    for path in files:
        rel = path.relative_to(root)
        paths = output_paths(out, rel)
        for target in paths.values():
            target.parent.mkdir(parents=True, exist_ok=True)
        mask_path, meta_path, preview_path = paths["mask"], paths["json"], paths["preview"]
        error_path = meta_path.with_name(rel.stem + ".error.json")
        if not args.overwrite and all(x.exists() for x in paths.values()):
            print(f"SKIP {rel}; use --overwrite to regenerate")
            continue
        try:
            p = DEFAULTS | config.get("defaults", {}) | config.get("overrides", {}).get(rel.as_posix(), {}) | cli
            validate(p)
            rgb, wh, mpp, notes = read_image(path, p)
            if args.full_size and math.prod(wh) > args.max_full_pixels:
                raise ValueError("Full-size mask exceeds max-full-pixels; use working-resolution mask")
            mask, info = extract_mask(rgb, p)
            image = Image.fromarray(mask.astype(np.uint8) * 255)
            if args.full_size:
                image = image.resize(wh, Image.Resampling.NEAREST)
            # Commit metadata last; it describes only successfully saved outputs.
            meta_path.unlink(missing_ok=True)
            image.save(mask_path, format="TIFF", compression="tiff_deflate")
            save_preview(rgb, mask, preview_path)
            info["warnings"] += notes
            meta = dict(input=str(path), parameters=p, diagnostics=info,
                        tiatoolbox_version=importlib.metadata.version("tiatoolbox"),
                        original_wh=list(wh), working_wh=[mask.shape[1], mask.shape[0]],
                        mask_wh=list(image.size), working_mpp_xy=mpp,
                        baseline_pixels_per_mask_pixel_xy=[wh[0]/image.width, wh[1]/image.height],
                        mask_values={"background": 0, "tissue": 255},
                        outputs={k: str(v) for k, v in paths.items()},
                        coordinates="Origin top-left; x=column, y=row; no crop or rotation")
            meta_path.write_text(json.dumps(meta, indent=2, allow_nan=False), encoding="utf-8")
            error_path.unlink(missing_ok=True)
            print(f"OK {rel} [{p['mask_mode']}]: tissue={info['tissue_fraction']:.1%}; "
                  + "; ".join(info["warnings"]))
        except Exception as exc:
            failures += 1
            error_path.parent.mkdir(parents=True, exist_ok=True)
            error_path.write_text(json.dumps({"input": str(path), "error": str(exc)}, indent=2), encoding="utf-8")
            print(f"ERROR {rel}: {exc}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
