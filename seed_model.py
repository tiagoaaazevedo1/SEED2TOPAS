"""A seed model is one plain TOPAS file: materials (Ma/), geometry (Ge/) and one source (So/), with the seed at the origin
and its axis along z. Lines such as "# sk_per_hist_Gy_cm2 = 4.03e-14" give what TOPAS does not need (see README).
block() repeats the seed for the plan: every Ge/ and So/ name gets the suffix _<i> and the seed is moved to its position."""
import os
import re

KEYS = ("name", "isotope", "half_life_d", "aliases", "sk_per_hist_Gy_cm2", "lambda_wafac", "status", "note", "references")
PARAM = re.compile(r"^\s*\w+:(Ge|So|Ma)/([^/\s]+)/([^\s=]+)\s*=\s*(.*)$")
OTHER = re.compile(r"^\s*\w+:(\w+)/")
SCORED = re.compile(r'^\s*\w+:Sc/[^/]+/Component\s*=\s*"([^"]+)"', re.I)


def _value(rhs):
    return rhs.split("#")[0].strip().strip('"')


class SeedModel:
    def __init__(self, path):
        self.file, self.path = os.path.basename(path), path
        self.meta, self.notes, ma, ge, so, scored, ignored = {}, [], [], {}, {}, set(), set()
        for raw in open(path, encoding="utf-8"):
            ln = raw.rstrip("\n")
            m = re.match(r"^\s*#\s*(\w+)\s*=\s*(.*)$", ln)
            if m and m.group(1) in KEYS:
                self.meta[m.group(1)] = m.group(2).strip()
            if not ln.strip() or ln.lstrip().startswith("#"):
                continue
            s = SCORED.match(ln)
            if s:
                scored.add(s.group(1))
            p = PARAM.match(ln)
            if not p:
                o = OTHER.match(ln)
                if o:
                    ignored.add(o.group(1))
                continue
            kind, name = p.group(1), p.group(2)
            if kind == "Ma":
                ma.append(ln)
            elif kind == "Ge":
                ge.setdefault(name, []).append(ln)
            else:
                so.setdefault(name, []).append(ln)
        ge.pop("World", None)
        self.materials = "\n".join(ma)

        field = lambda lines, key: next((_value(PARAM.match(x).group(4)) for x in lines if PARAM.match(x).group(3).lower() == key), None)
        ge = {k: v for k, v in ge.items() if field(v, "type")}
        parent = {k: field(v, "parent") for k, v in ge.items()}
        # components that only exist to score (a detector, a scoring grid), and everything inside them
        drop, grew = {c for c in scored if c in ge}, True
        while grew:
            more = {k for k, p in parent.items() if p in drop} - drop
            drop |= more
            grew = bool(more)
        for c in drop:
            ge.pop(c)
            parent.pop(c)
        if drop:
            self.notes.append(f"ignored the scoring components in the seed file: {', '.join(sorted(drop))}")
        if ignored:
            self.notes.append("ignored the parameters in the seed file other than Ma/, Ge/ and So/: " + ", ".join(sorted(ignored)))

        if not ge:
            raise ValueError(f"{self.file}: no geometry (Ge/) found")
        if len(so) != 1:
            raise ValueError(f"{self.file}: needs exactly one source (So/), found {len(so)}")
        self.source, self.so_lines = next(iter(so.items()))
        comp = field(self.so_lines, "component")
        if comp not in ge:
            raise ValueError(f'{self.file}: the source is in component "{comp}", which is not in the file')
        self.ge, self.parent = ge, parent
        self.tops = [k for k, p in parent.items() if p not in ge]
        one_group = len(self.tops) == 1 and field(ge[self.tops[0]], "type").lower() == "group"
        self.root = self.tops[0] if one_group else None
        # solids to draw: a group draws nothing, so its direct children
        self.visible = [k for k, p in parent.items() if p == self.root] if one_group else list(self.tops)
        self.names = "|".join(re.escape(n) for n in sorted(list(ge) + [self.source], key=len, reverse=True))

    @property
    def name(self):
        return self.meta.get("name") or os.path.splitext(self.file)[0]

    @property
    def status(self):
        return self.meta.get("status", "not validated")

    @property
    def sk_per_hist(self):
        try:
            return float(self.meta.get("sk_per_hist_Gy_cm2", ""))
        except ValueError:
            return None

    @property
    def half_life_d(self):
        try:
            return float(self.meta.get("half_life_d", ""))
        except ValueError:
            return None

    @property
    def aliases(self):
        return [a.strip().lower() for a in self.meta.get("aliases", "").split(",") if a.strip()]

    def matches(self, source_name):
        """True / False against the source name in the RT Plan; None if the file has no aliases to compare."""
        if not self.aliases:
            return None
        n = str(source_name).lower().replace(".", "")
        # "s17 plus" in the plan must not match a model that only has the alias "s17"
        return any(a in n for a in self.aliases) and not ("plus" in n and not any("plus" in a for a in self.aliases))

    def _rename(self, ln, i, parent_of_top):
        p = PARAM.match(ln)
        key = p.group(3).lower()
        ln = re.sub(rf"\b(Ge|So)/({self.names})/", lambda m: f"{m.group(1)}/{m.group(2)}_{i}/", ln)
        if p.group(1) == "Ge" and key == "parent" and p.group(2) in self.tops:
            return re.sub(r'=.*$', f'= "{parent_of_top}"', ln)
        if key in ("parent", "component"):
            ln = re.sub(r'"([^"]+)"', lambda m: f'"{m.group(1)}_{i}"' if m.group(1) in self.ge else m.group(0), ln)
        return ln

    def block(self, i, x, y, z, n, parent, parallel):
        trans = [f"d:Ge/{{}}/Trans{a} = {v:.4f} mm" for a, v in zip("XYZ", (x, y, z))]
        out = [f"# ---- seed {i}"]
        if self.root:
            top = parent
        else:                   # several solids at the top of the file: a group places them
            top = f"SeedPos_{i}"
            out += [f's:Ge/{top}/Type = "Group"', f's:Ge/{top}/Parent = "{parent}"'] + [t.format(top) for t in trans]
            if parallel:
                out.append(f'b:Ge/{top}/IsParallel = "True"')
        for name, lines in self.ge.items():
            for ln in lines:
                key = PARAM.match(ln).group(3).lower()
                if key == "isparallel" or (name == self.root and key in ("transx", "transy", "transz")):
                    continue
                out.append(self._rename(ln, i, top))
                if key == "type" and parallel:
                    out.append(f'b:Ge/{name}_{i}/IsParallel = "True"')
                if key == "parent" and name == self.root:
                    out += [t.format(f"{name}_{i}") for t in trans]
        for ln in self.so_lines:
            if PARAM.match(ln).group(3).lower() != "numberofhistoriesinrun":
                out.append(self._rename(ln, i, top))
        out.append(f"ic:So/{self.source}_{i}/NumberOfHistoriesInRun = {n}")
        return "\n".join(out)
