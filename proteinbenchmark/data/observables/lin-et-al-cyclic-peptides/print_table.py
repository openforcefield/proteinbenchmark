#!/usr/bin/env python3
"""Print NMR tables from JSON files for manual validation against SI sources."""
import dataclasses
import json
import glob
import sys

def main():
    files = sorted(glob.glob("*.json"))
    if not files:
        sys.exit("No JSON files found in current directory.")
    print_tables(files)

ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","GLU":"E","GLY":"G","ILE":"I",
       "LEU":"L","LYS":"K","PHE":"F","PRO":"P","SER":"S","THR":"T","VAL":"V"}
BB_H = {"HN","HA","HA2","HA3","HB","HB2","HB3","HB#"}
BB_C = {"CA","CB"}

def zip_strict(*args, **kwargs):
    return zip(*args, **kwargs, strict=True)

def label(r: dict[str, str], use_tlc: bool = False) -> str:
    tlc = r["three_letter_code"]
    olc = ONE[tlc]
    if use_tlc:
        return ("ᴅ-" if r["chirality"]=="D" else "  ") + tlc
    else:
        return olc.lower() if r["chirality"]=="D" else olc.upper()

def fmt_h(names: set[str], H: dict[str, list[float]]) -> list[str]:
    vals = [(n, shift) for n in sorted(names) if n in H for shift in H[n]]
    if not vals:
        return []
    if len(vals) == 1:
        return [f"{vals[0][1]:.3f}"]
    return [f"{n}:{v:.3f}" for n,v in vals]

@dataclasses.dataclass
class LinesBySource:
    """Assists in sorting tables"""
    _d: dict[str, list[str]] = dataclasses.field(default_factory=dict[str, list[str]])
    _cursor: str | None = None

    def set_cursor(self, src: str) -> None:
        self._cursor = src

    def add_line(self, line: str) -> None:
        assert self._cursor is not None
        assert isinstance(self._cursor, str)
        self._d.setdefault(self._cursor, []).append(line)

    def multiline(self, cols, widths):
        """Zip columns (each a list of strings) into printed lines."""
        height = max(len(c) for c in cols)
        for row in range(height):
            self.add_line("  " + "  ".join(f"{(cols[c][row] if row < len(cols[c]) else ''):<{widths[c]}}" for c in range(len(cols))))


    def print(self, *args, **kwargs) -> None:
        for _, lines in sorted(self._d.items(), key=lambda kv: int(kv[0].rsplit("S", 1)[-1])):
            for line in lines:
                print(line, *args, **kwargs)
        self._cursor = None

def print_tables(files: list[str]):
    for fname in files:
        data = json.load(open(fname))
        if "peptides" not in data:
            continue

        print(f"{'='*90}")
        print(fname)

        lines_by_source = LinesBySource()
        for pep in data["peptides"]:
            seq  = pep["sequence"]
            cond = pep["experimental_conditions"]
            lbl  = {r["position"]: label(r) for r in seq}
            lbl3 = {r["position"]: label(r, use_tlc=True) for r in seq}
            def cv(f): return cond[f]["value"] if cond and cond.get(f) else "?"

            # ── Chemical shifts ───────────────────────────────────────────────────
            cs = pep.get("chemical_shifts")
            if cs:
                # dict[residx, dict[atomname, dict[list[float]]]]
                H: dict[int, dict[str, list[float]]]={}
                C: dict[int, dict[str, list[float]]]={}
                N: dict[int, dict[str, list[float]]]={}
                for blk, dst in [("proton_shifts",H),("carbon_shifts",C),("nitrogen_shifts",N)]:
                    if cs.get(blk):
                        for ri,an,sh in zip_strict(cs[blk]["residue_index"],cs[blk]["atom_name"],cs[blk]["shift"]):
                            dst.setdefault(ri,{}).setdefault(an,[]).append(sh["value"])
                # dict[residx, float]
                dHN = {}
                if cs.get("amide_proton_temperature_coefficients"):
                    for ri,co in zip_strict(cs["amide_proton_temperature_coefficients"]["residue_index"],
                                    cs["amide_proton_temperature_coefficients"]["coefficient"]):
                        dHN[ri] = co["value"]

                jc = pep.get("j_couplings_and_torsional_restraints")
                J = {}
                if jc:
                    for ri,an,jv,pc,pr in zip_strict(jc["residue_index"],jc["ha_atom_name"],
                                            jc["j_coupling"],jc["phi_center"],jc["phi_half_range"]):
                        phi = f" {pc['value']:+.0f}±{pr['value']:.0f}°" if pc and pr else ""
                        J.setdefault(ri,[]).append(f"{an}:{jv['value']:.1f}{phi}" if jv else f"{an}:—")

                W = [6, 8, 8, 12, 7, 12, 7, 16, 18, 12, 20, 8]
                hdr = ["Res","HN","15N","Ha","13Ca","Hb","13Cb","Other H","Other C","Other N","J (Hz)","dHN/dT"]
                src = f"{cs['source'].get('table','')}"
                lines_by_source.set_cursor(src)
                lines_by_source.add_line(f"\n  {pep['name']}  Chemical shifts [{src} {fname}]")
                lines_by_source.add_line(f"  {cond['solvent'] if cond else '?'}  pH {cv('ph')}  {cv('temperature')} °C  {cv('spectrometer_frequency')} MHz")
                lines_by_source.add_line("  " + "  ".join(f"{h:<{w}}" for h,w in zip_strict(hdr,W)))
                lines_by_source.add_line("  " + "  ".join("-"*w for w in W))
                for pos in sorted(r["position"] for r in seq):
                    h, c, n = (H.get(pos,{}), C.get(pos,{}), N.get(pos,{}))
                    oh = [f"{an}:{shift:.3f}" for an in sorted(h) if an not in BB_H for shift in h[an]]
                    oc = [f"{an}:{shift:.2f}"  for an in sorted(c) if an not in BB_C for shift in c[an]]
                    on = [f"{an}:{shift:.2f}"  for an in sorted(n) if an != "N" for shift in n[an]]
                    cols = [
                        [f"{lbl3[pos]}{pos}"],
                        [f"{shift:.3f}" for shift in h['HN']] if "HN" in h else [],
                        [f"{shift:.2f}" for shift in n['N']]  if "N"  in n else [],
                        fmt_h({"HA","HA2","HA3"}, h),
                        [f"{shift:.2f}"for shift in c['CA']] if "CA" in c else [],
                        fmt_h({"HB","HB2","HB3","HB#"}, h),
                        [f"{shift:.2f}" for shift in c['CB']] if "CB" in c else [],
                        oh, oc, on,
                        J.get(pos,[]),
                        [f"{dHN[pos]:.1f}"] if pos in dHN else [],
                    ]
                    lines_by_source.multiline(cols, W)
                    lines_by_source.add_line("")  # blank line between residues for readability

            # ── Observed NOEs ─────────────────────────────────────────────────────
            onoe = pep.get("observed_noes")
            if onoe:
                n = len(onoe["residue_index_1"])
                has_strength = any(s is not None for s in onoe["strength"])
                src = f"{onoe['source'].get('table','')}"
                lines_by_source.set_cursor(src)
                lines_by_source.add_line(f"\n  {pep['name']}  Observed NOEs ({n} total)  [{src} {fname}]")
                lines_by_source.add_line(f"  {cond['solvent'] if cond else '?'}  pH {cv('ph')}  {cv('temperature')} °C  {cv('spectrometer_frequency')} MHz")
                w1 , w2, w3, w4 = (16, 8, 16, 8)
                hdr2 = f"  {'ppm':>{w2}} {'atom 1':<{w1}} {'ppm':>{w4}} {'atom 2':<{w3}}"
                if has_strength:
                    hdr2 += "  strength"
                lines_by_source.add_line(hdr2)
                lines_by_source.add_line(f"  {'-'*w2} {'-'*w1} {'-'*w4} {'-'*w3}" + ("  --------" if has_strength else ""))
                for i in range(n):
                    a1  = f"{onoe['residue_index_1'][i]}{lbl[onoe['residue_index_1'][i]]}-{onoe['atom_name_1'][i]}"
                    a2  = f"{onoe['residue_index_2'][i]}{lbl[onoe['residue_index_2'][i]]}-{onoe['atom_name_2'][i]}"
                    s1  = onoe["chemical_shift_1"][i]
                    s2  = onoe["chemical_shift_2"][i]
                    p1  = f"{s1['value']:.2f}" if s1 else ""
                    p2  = f"{s2['value']:.2f}" if s2 else ""
                    row = f"  {p1:>{w2}} {a1:<{w1}} {p2:>{w4}} {a2:<{w3}}"
                    if has_strength:
                        row += f"  {onoe['strength'][i] or ''}"
                    lines_by_source.add_line(row)


            # ── J-couplings / torsional restraints ───────────────────────────────
            jc = pep.get("j_couplings_and_torsional_restraints")
            if jc:
                n = len(jc["residue_index"])
                src = f"{jc['source'].get('table','')}"
                lines_by_source.set_cursor(src)
                lines_by_source.add_line(f"\n  {pep['name']}  J-couplings  [{src} {fname}]")
                lines_by_source.add_line(f"  {cond['solvent'] if cond else '?'}  pH {cv('ph')}  {cv('temperature')} °C  {cv('spectrometer_frequency')} MHz")
                lines_by_source.add_line(f"  {'Res':<6} {'atom':<6} {'J (Hz)':>7}  phi restraint")
                lines_by_source.add_line(f"  {'-'*6} {'-'*6} {'-'*7}  {'-'*20}")
                for i in range(n):
                    ri = jc["residue_index"][i]
                    an = jc["ha_atom_name"][i]
                    jv = jc["j_coupling"][i]
                    pc = jc["phi_center"][i]
                    pr = jc["phi_half_range"][i]
                    j_str   = f"{jv['value']:.1f}" if jv else "—"
                    phi_str = f"{pc['value']:+.0f} ± {pr['value']:.0f}°" if pc and pr else "—"
                    lines_by_source.add_line(f"  {lbl[ri]}{ri:<5} {an:<6} {j_str:>7}  {phi_str}")
                lines_by_source.add_line("")

            # ── NOE distance restraints ───────────────────────────────────────────
            noe = pep.get("noe_derived_distance_restraints")
            if noe:
                n = len(noe["residue_index_1"])
                src = f"{noe['source'].get('table','')}"
                lines_by_source.set_cursor(src)
                lines_by_source.add_line(f"  {pep['name']}  NOE distance restraints ({n} total)  [{src} {fname}]")
                lines_by_source.add_line(f"  {cond['solvent'] if cond else '?'}  pH {cv('ph')}  {cv('temperature')} °C  {cv('spectrometer_frequency')} MHz")
                lines_by_source.add_line(f"  {'#':<4} {'atom 1':<16} {'atom 2':<16} {'target':>8} {'lo':>7} {'up':>7}")
                lines_by_source.add_line(f"  {'-'*4} {'-'*16} {'-'*16} {'-'*8} {'-'*7} {'-'*7}")
                for i in range(n):
                    a1 = f"{noe['residue_index_1'][i]}{lbl[noe['residue_index_1'][i]]}-{noe['atom_name_1'][i]}"
                    a2 = f"{noe['residue_index_2'][i]}{lbl[noe['residue_index_2'][i]]}-{noe['atom_name_2'][i]}"
                    d,lo,up = noe["target_distance"][i]["value"],noe["lower_bound_delta"][i]["value"],noe["upper_bound_delta"][i]["value"]
                    lines_by_source.add_line(f"  {i+1:<4} {a1:<16} {a2:<16} {d:>8.3f} {lo:>7.3f} {up:>7.3f}")

        lines_by_source.print()
        print()

if __name__ == "__main__":
    main()
