#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from Bio import SeqIO
from PIL import Image, ImageDraw, ImageFont

# Base -> code (small integers)
# 0:A 1:T 2:G 3:C 4:N 5:- 6:other
BASE_TO_CODE = np.full(256, 6, dtype=np.uint8)
for ch, code in [
    (b"A", 0),
    (b"a", 0),
    (b"T", 1),
    (b"t", 1),
    (b"G", 2),
    (b"g", 2),
    (b"C", 3),
    (b"c", 3),
    (b"N", 4),
    (b"n", 4),
    (b"-", 5),
]:
    BASE_TO_CODE[ch[0]] = code
# Set1-inspired palette (A,T,G,C,N,-,other)
PALETTE_RGB = np.array(
    [
        [228, 26, 28],  # A red
        [77, 175, 74],  # T green
        [255, 127, 0],  # G orange
        [55, 126, 184],  # C blue
        [153, 153, 153],  # N grey
        [245, 245, 245],  # - near-white
        [200, 200, 200],  # other light grey
    ],
    dtype=np.uint8,
)
GAP = 5


def parse_args():
    p = argparse.ArgumentParser(
        prog="msa_overview.py",
        description="Render a whole-alignment FASTA overview (downsampled) as PNG or PDF.",
    )
    p.add_argument("alignment", help="Aligned FASTA file (all sequences same length).")
    p.add_argument(
        "-o",
        "--out",
        default=None,
        help="Output filename. Default: <basename>.overview.<format>",
    )
    p.add_argument(
        "--format", choices=["png", "pdf"], default="png", help="Output format."
    )
    p.add_argument(
        "--width",
        type=int,
        default=4000,
        help="Output image width in pixels (alignment area).",
    )
    p.add_argument(
        "--height",
        type=int,
        default=2000,
        help="Output image height in pixels (alignment area).",
    )
    p.add_argument(
        "--max-seqs",
        type=int,
        default=None,
        help="Optional cap on sequences to read (for huge files).",
    )
    # Left labels / tracking
    p.add_argument(
        "--label-margin",
        type=int,
        default=120,
        help="Left margin (pixels) to draw row numbers/ticks. Set 0 to disable.",
    )
    p.add_argument(
        "--label-step",
        type=int,
        default=50,
        help="Label every Nth sequence row (1-based).",
    )
    p.add_argument(
        "--label-map",
        default=None,
        help="Write TSV mapping of rows/metrics. Default: <basename>.rowmap.tsv",
    )
    # Top scale / bases
    p.add_argument(
        "--top-margin",
        type=int,
        default=60,
        help="Top margin (pixels) to draw base-position scale. Set 0 to disable.",
    )
    p.add_argument(
        "--x-step",
        type=int,
        default=None,
        help="Tick step in bases (e.g. 1000). Default: auto based on alignment length.",
    )
    # Font
    p.add_argument("--font-size", type=int, default=12, help="Font size for labels.")
    p.add_argument(
        "--font",
        default=None,
        help="Path to a .ttf font file (optional). If not set, uses PIL default.",
    )
    # Sorting
    p.add_argument(
        "--sort",
        choices=["none", "ungapped", "consensus"],
        default="none",
        help="Sort sequences by ungapped length or by identity to consensus.",
    )
    p.add_argument(
        "--consensus-ignore-gaps",
        action="store_true",
        help="When building consensus, ignore '-' in each column (recommended).",
    )
    p.add_argument(
        "--identity-ignore-gaps",
        choices=["both", "ref", "query", "none"],
        default="both",
        help="How to treat gaps when computing identity to consensus. Default: both.",
    )
    return p.parse_args()


def read_alignment_as_codes_and_ids(fasta_path, max_seqs=None):
    seq_codes = []
    ids = []
    ncol = None
    for i, rec in enumerate(SeqIO.parse(fasta_path, "fasta")):
        if max_seqs is not None and i >= max_seqs:
            break

        s = str(rec.seq)
        if ncol is None:
            ncol = len(s)
            if ncol == 0:
                raise ValueError("First sequence has length 0.")

        elif len(s) != ncol:
            raise ValueError(
                f"Not a rectangular alignment: sequence {rec.id} has length {len(s)} vs {ncol}."
            )

        b = s.encode("ascii", "replace")
        arr = np.frombuffer(b, dtype=np.uint8)
        codes = BASE_TO_CODE[arr]
        seq_codes.append(codes)
        ids.append(rec.id)
    if not seq_codes:
        raise ValueError("No sequences found in FASTA.")

    return np.stack(seq_codes, axis=0), ids


def ungapped_lengths(codes):
    return (codes != GAP).sum(axis=1)


def consensus_reference(codes, ignore_gaps=True, n_codes=7):
    n_seqs, n_cols = codes.shape
    ref = np.empty(n_cols, dtype=np.uint8)
    for j in range(n_cols):
        col = codes[:, j]
        if ignore_gaps:
            col = col[col != GAP]
            if col.size == 0:
                ref[j] = GAP
                continue

        ref[j] = np.bincount(col, minlength=n_codes).argmax()
    return ref


def identity_to_ref(codes, ref_codes, ignore_gaps="both"):
    ref = ref_codes[None,:]
    q = codes
    if ignore_gaps == "both":
        mask = (q != GAP) & (ref != GAP)
    elif ignore_gaps == "ref":
        mask = (ref != GAP)
    elif ignore_gaps == "query":
        mask = (q != GAP)
    elif ignore_gaps == "none":
        mask = np.ones_like(q, dtype=bool)
    else:
        raise ValueError("ignore_gaps must be one of: both, ref, query, none")

    denom = mask.sum(axis=1)
    matches = ((q == ref) & mask).sum(axis=1)
    out = np.full(q.shape[0], np.nan, dtype=float)
    ok = denom > 0
    out[ok] = matches[ok] / denom[ok]
    return out


def sort_codes_and_ids(
    codes, ids, sort_mode, consensus_ignore_gaps, identity_ignore_gaps
):
    n = codes.shape[0]
    ung = ungapped_lengths(codes)
    ref = consensus_reference(codes, ignore_gaps=consensus_ignore_gaps)
    ident = identity_to_ref(codes, ref, ignore_gaps=identity_ignore_gaps)
    if sort_mode == "none":
        order = np.arange(n)
    elif sort_mode == "ungapped":
        order = np.argsort(-ung)
    elif sort_mode == "consensus":
        ident_key = np.nan_to_num(ident, nan=-1.0)
        order = np.lexsort((-ung, -ident_key))
    else:
        raise ValueError(f"Unknown sort mode: {sort_mode}")

    codes2 = codes[order,:]
    ids2 = [ids[i] for i in order]
    metrics = {
        "order": order,
        "original_row": order + 1,  # 1-based original row
        "ungapped_length": ung[order],
        "identity_to_consensus": ident[order],
    }
    return codes2, ids2, metrics


def downsample_mode_2d(codes, out_h, out_w, n_codes=7):
    H, W = codes.shape
    out_h = min(out_h, H)
    out_w = min(out_w, W)
    # --- Downsample rows ---
    row_factor = int(np.ceil(H / out_h))
    H_pad = out_h * row_factor
    if H_pad != H:
        codes_r = np.pad(
            codes,
            ((0, H_pad - H), (0, 0)),
            mode="constant",
            constant_values=GAP,  # <-- key change
        )
    else:
        codes_r = codes
    codes_r = codes_r.reshape(out_h, row_factor, W)
    counts = np.zeros((out_h, W, n_codes), dtype=np.uint16)
    for k in range(n_codes):
        counts[:,:, k] = (codes_r == k).sum(axis=1)
    ds_rows = counts.argmax(axis=2).astype(np.uint8)
    # --- Downsample cols ---
    col_factor = int(np.ceil(W / out_w))
    W_pad = out_w * col_factor
    if W_pad != W:
        ds_rows = np.pad(
            ds_rows,
            ((0, 0), (0, W_pad - W)),
            mode="constant",
            constant_values=GAP,  # <-- key change
        )
    ds_rows = ds_rows.reshape(out_h, out_w, col_factor)
    counts2 = np.zeros((out_h, out_w, n_codes), dtype=np.uint16)
    for k in range(n_codes):
        counts2[:,:, k] = (ds_rows == k).sum(axis=2)
    return counts2.argmax(axis=2).astype(np.uint8)


def load_font(font_path, font_size):
    if font_path:
        try:
            return ImageFont.truetype(font_path, font_size)

        except Exception as e:
            raise RuntimeError(f"Could not load font '{font_path}': {e}")

    return ImageFont.load_default()


def nice_step(n_cols):
    """
    Choose a human-friendly tick step given alignment length.
    Produces ~5–10 ticks typically.
    """
    if n_cols <= 0:
        return 1

    target_ticks = 8
    raw = max(1, int(round(n_cols / target_ticks)))
    # round raw to 1,2,5 * 10^k
    pow10 = 10 ** int(np.floor(np.log10(raw)))
    for m in (1, 2, 5, 10):
        step = m * pow10
        if step >= raw:
            return step

    return 10 * pow10


def add_row_labels(rgb_img, n_seqs, label_margin, label_step, font):
    if label_margin <= 0:
        return rgb_img

    w, h = rgb_img.size
    canvas = Image.new("RGB", (w + label_margin, h), (255, 255, 255))
    canvas.paste(rgb_img, (label_margin, 0))
    draw = ImageDraw.Draw(canvas)
    for i in range(1, n_seqs + 1, max(1, label_step)):
        y = int(round((i - 1) * (h / n_seqs)))
        draw.line(
            [(label_margin - 8, y), (label_margin - 1, y)], fill=(0, 0, 0), width=1
        )
        draw.text((2, max(0, y - 7)), str(i), fill=(0, 0, 0), font=font)
    draw.line([(label_margin - 1, 0), (label_margin - 1, h)], fill=(0, 0, 0), width=1)
    return canvas


def add_top_scale(img, n_cols, top_margin, label_margin, font, x_step=None):
    """
    Add a top margin and draw a base-position scale.
    img: current image (may already include left margin)
    n_cols: original alignment columns
    label_margin: left margin width (so scale starts at alignment area)
    """
    if top_margin <= 0:
        return img

    w, h = img.size
    canvas = Image.new("RGB", (w, h + top_margin), (255, 255, 255))
    canvas.paste(img, (0, top_margin))
    draw = ImageDraw.Draw(canvas)
    # Alignment area spans x from label_margin to w
    x0 = label_margin
    x1 = w - 1
    axis_y = top_margin - 15  # axis line a bit above the alignment
    # Axis line
    draw.line([(x0, axis_y), (x1, axis_y)], fill=(0, 0, 0), width=1)
    if x_step is None:
        x_step = nice_step(n_cols)
    # ticks at positions 1..n_cols
    # map base position p (1-based) to pixel:
    # px = x0 + (p-1) * (alignment_px_width / n_cols)
    align_w = max(1, (x1 - x0))

    def pos_to_x(p):
        return int(round(x0 + (p - 1) * (align_w / max(1, n_cols))))

    # Always draw first and last
    positions = list(range(1, n_cols + 1, x_step))
    if positions[-1] != n_cols:
        positions.append(n_cols)
    for p in positions:
        x = pos_to_x(p)
        draw.line([(x, axis_y), (x, axis_y + 6)], fill=(0, 0, 0), width=1)
        label = str(p)
        # simple label placement; avoid drawing outside
        tx = max(0, min(w - 1, x - 10))
        ty = max(0, axis_y - 14)
        draw.text((tx, ty), label, fill=(0, 0, 0), font=font)
    # Title-ish label
    draw.text((x0, 2), "Position (bp)", fill=(0, 0, 0), font=font)
    return canvas


def write_rowmap(path, ids, metrics):
    ung = metrics["ungapped_length"]
    ident = metrics["identity_to_consensus"]
    orig_row = metrics["original_row"]
    with open(path, "w", encoding="utf-8") as f:
        f.write(
            "sorted_row\toriginal_row\tsequence_id\tungapped_length\tidentity_to_consensus\n"
        )
        for i, sid in enumerate(ids, start=1):
            idv = ident[i - 1]
            idv_str = "" if np.isnan(idv) else f"{idv:.6f}"
            f.write(f"{i}\t{orig_row[i - 1]}\t{sid}\t{int(ung[i - 1])}\t{idv_str}\n")


def save_image(im, out_path, fmt):
    fmt = fmt.lower()
    if fmt == "png":
        im.save(out_path, format="PNG")
    elif fmt == "pdf":
        im.convert("RGB").save(out_path, format="PDF")
    else:
        raise ValueError(f"Unknown format: {fmt}")


def main():
    args = parse_args()
    base = os.path.splitext(os.path.basename(args.alignment))[0]
    out = args.out or f"{base}.overview.{args.format}"
    label_map = args.label_map
    if label_map is None:
        label_map = f"{base}.rowmap.tsv"
    codes, ids = read_alignment_as_codes_and_ids(args.alignment, max_seqs=args.max_seqs)
    codes, ids, metrics = sort_codes_and_ids(
        codes,
        ids,
        sort_mode=args.sort,
        consensus_ignore_gaps=args.consensus_ignore_gaps,
        identity_ignore_gaps=args.identity_ignore_gaps,
    )
    n_seqs, n_cols = codes.shape
    ds = downsample_mode_2d(codes, out_h=args.height, out_w=args.width)
    rgb = PALETTE_RGB[ds]
    im = Image.fromarray(rgb, mode="RGB")
    font = load_font(args.font, args.font_size)
    # Left row ticks
    if args.label_margin > 0:
        im = add_row_labels(
            im,
            n_seqs=n_seqs,
            label_margin=args.label_margin,
            label_step=args.label_step,
            font=font,
        )
    # Top base scale (needs to know label_margin so axis starts at alignment)
    if args.top_margin > 0:
        im = add_top_scale(
            im,
            n_cols=n_cols,
            top_margin=args.top_margin,
            label_margin=max(0, args.label_margin),
            font=font,
            x_step=args.x_step,
        )
    if label_map:
        write_rowmap(label_map, ids, metrics)
    save_image(im, out, args.format)
    print(f"Wrote: {out}")
    print(f"Row map: {label_map}")
    print(f"Sort: {args.sort}")
    print(f"Input: {n_seqs} sequences × {n_cols} columns")
    print(f"Overview (alignment area): {args.width} × {args.height} px")
    print(
        f"Top scale: {'on' if args.top_margin > 0 else 'off'}  (x_step={args.x_step or nice_step(n_cols)})"
    )


if __name__ == "__main__":
    main()
