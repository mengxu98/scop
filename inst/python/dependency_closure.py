"""Read the active dependency closure without importing analysis packages."""

import sys
from importlib import metadata

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name


def dependency_closure(requirements):
    """Follow Requires-Dist, including selected extras and environment markers.

    A conflict belongs to the distribution imposing the requirement. Reverse
    dependents outside this closure may impose incompatible constraints on a
    shared dependency without invalidating the requested workflow's constraints.
    """
    pending = [Requirement(item) for item in requirements]
    selected = {}
    while pending:
        request = pending.pop()
        name = canonicalize_name(request.name)
        extras = set(request.extras)
        if name in selected and extras <= selected[name]:
            continue
        selected.setdefault(name, set()).update(extras)
        try:
            distribution = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            continue
        for raw in distribution.requires or ():
            dependency = Requirement(raw)
            if dependency.marker is None or any(
                dependency.marker.evaluate({"extra": extra})
                for extra in selected[name] | {""}
            ):
                pending.append(dependency)
    return sorted(selected)


if __name__ == "__main__":
    print("\n".join(dependency_closure(sys.argv[1:])))
